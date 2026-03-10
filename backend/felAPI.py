from fastapi import FastAPI, Body, HTTPException, APIRouter
from pydantic import BaseModel
from typing import Any, Callable, Dict, List, Optional
from ebeam import beam
from beamline import *
from schematic import *
from fastapi.middleware.cors import CORSMiddleware
import asyncio
import inspect
import importlib
import io
import base64
import logging
import tempfile
import threading
import traceback
import types
from excelElements import ExcelElements
import uvicorn
import os
from dotenv import load_dotenv
import copy
import math
import numpy as np
import matplotlib.pyplot as plt

load_dotenv('../.env')  # Only during dev testing when not using Dockerfile...
FRONTEND_PORT = os.getenv('FRONTEND_PORT', '5173')

ORIGINS = [f'http://localhost:{FRONTEND_PORT}', f"localhost:{FRONTEND_PORT}"]
#ORIGINS = ["http://localhost:5173", "localhost:5173"]
moduleName = 'beamline'

ALLOWED_BEAMLINE_CLASSES = {'driftLattice', 'qpfLattice', 'qpdLattice', 'dipole', 'dipole_wedge', 'lattice'}


def _instantiate_beamline_class(module, class_name: str, parameters: dict):
    if class_name not in ALLOWED_BEAMLINE_CLASSES:
        raise HTTPException(status_code=400, detail=f"Unknown beamline class: '{class_name}'")
    cls = getattr(module, class_name)
    return cls(**parameters)

class AxisTwiss(BaseModel):
    alpha: float
    beta: float
    phi: float
    epsilon: float

class TwissParameters(BaseModel):
    x: AxisTwiss
    y: AxisTwiss
    z: AxisTwiss

class BeamlineInfo(BaseModel):
    #__root__: Dict[str, Dict[str, Any]]
    segmentName: str
    parameters: Dict[str, Any]

MAX_PARTICLES = 100_000

class PlottingParameters(BaseModel):
    beamlineData: list[BeamlineInfo]
    beamType: str = 'electron'
    num_particles: int
    kineticE: float = 45.0
    interval: float = 1
    defineLim: bool = True
    saveData: bool = False
    matchScaling: bool = True
    scatter: bool = True
    beam_setup: str = 'twiss'
    twiss: TwissParameters = None
    #  I THINK WE NEED SAVE FIG AND SHAPE

class LineAxObject(BaseModel):
    axis: str # temporary placeholder axes
    twiss: str
    x_axis: list[float]
    beamsegment: list

class AxesPNGData(BaseModel):
    images: Dict[float, Any]
    line_graph: LineAxObject

class GraphParameters(BaseModel):
    beam_index: int
    target_parameter: str
    target_s_pos: float
    beamline_data: list[BeamlineInfo]
    min: int | float = 0
    max: int | float = 10
    custom_step: int | float = 1

class GraphPlotPointResponse(BaseModel):
    x: float | None
    y: float | None
    z: float | None
    twiss_parameter: str

class GraphPlotData(BaseModel):
    parameter_value: float
    data: List[GraphPlotPointResponse]


# ---------------------------------------------------------------------------
# Glyfada HttpEvaluator endpoint
# ---------------------------------------------------------------------------

logger = logging.getLogger(__name__)

glyfada_router = APIRouter(prefix="/glyfada", tags=["glyfada"])


class GlyfadaSetupRequest(BaseModel):
    """Configuration for the evaluation context.

    The objective function is built from beamOptimizer internals: given a
    beamline, beam distribution, segment variables, objectives, and start
    point, we construct the MSE closure that maps parameter values to a
    scalar cost.
    """
    beamline_data: list[BeamlineInfo]
    beam_type: str = "electron"
    kinetic_energy: float = 45.0
    num_particles: int = 1000
    segment_var: Dict[str, Any]
    objectives: Dict[str, Any]
    start_point: Dict[str, Any]
    variable_names: list[str]
    n_objectives: int = 1
    n_constraints: int = 0
    emit_constraints: bool = False


class GlyfadaSetupFromPickleRequest(BaseModel):
    """Lightweight setup: point at an existing pickled objective."""
    pickle_path: str
    n_objectives: int = 1
    n_constraints: int = 0
    emit_constraints: bool = False


class GlyfadaEvalRequest(BaseModel):
    """Request body from HttpEvaluator: params dict + optional metadata."""
    params: Dict[str, float]
    timeout: Optional[float] = None
    fidelity: Optional[float] = None


class _GlyfadaState:
    """Mutable evaluation state stored on app.state."""
    def __init__(self):
        self.evaluate_fn: Optional[Callable[[Dict[str, float]], Dict[str, Any]]] = None
        self.variable_names: list[str] = []
        self.n_objectives: int = 1
        self.n_constraints: int = 0
        self.eval_count: int = 0
        self.lock = asyncio.Lock()          # guards state setup operations
        self.eval_lock = threading.Lock()    # guards CPU-bound evaluate_fn calls


@glyfada_router.get("/health")
async def glyfada_health():
    """Health check for HttpEvaluator's circuit breaker."""
    return {"status": "ok"}


@glyfada_router.get("/evaluate/health")
async def glyfada_evaluate_health():
    """Alias: HttpEvaluator appends /health to the configured URL, so if it
    targets /glyfada/evaluate the probe hits /glyfada/evaluate/health."""
    return {"status": "ok"}


@glyfada_router.post("/setup")
async def glyfada_setup(req: GlyfadaSetupRequest):
    """Configure the evaluation context from beamline specification.

    Builds the objective function closure used by beamOptimizer, wraps it
    in the DH response format (objective_1, constraint_1, ...).
    """
    state: _GlyfadaState = glyfada_router.state

    # Currently the closure produces a single scalar MSE objective.
    # Extend here when multi-objective support is added.
    if req.n_objectives != 1:
        raise HTTPException(
            status_code=400,
            detail=f"n_objectives must be 1 (got {req.n_objectives}). Multi-objective not yet supported.",
        )

    try:
        beamline_mod = importlib.import_module("beamline")
        beamlist = []
        for seg in req.beamline_data:
            beamlist.append(_instantiate_beamline_class(beamline_mod, seg.segmentName, seg.parameters))

        lat = lattice(1)
        beamlist = lat.changeBeamType(req.beam_type, req.kinetic_energy, beamlist)

        eb = beam()
        particles = eb.gen_6d_gaussian(0, [1, 1, 1, 1, 0.1, 100], req.num_particles)

        from beamOptimizer import beamOptimizer
        opt = beamOptimizer(beamlist, particles)

        # Convert string keys in segment_var/objectives to int (JSON keys are strings)
        seg_var = {int(k): v for k, v in req.segment_var.items()}
        objectives = {int(k): v for k, v in req.objectives.items()}
        start_point = req.start_point

        # Wire up the optimizer internals (mirrors calc() setup without running scipy)
        opt.segmentVar = seg_var
        seen = {}
        for idx in opt.segmentVar:
            name = opt.segmentVar[idx][0]
            if name not in seen:
                seen[name] = True
        opt.variablesToOptimize = list(seen.keys())

        opt.objectives = copy.deepcopy(objectives)
        for key, value in opt.objectives.items():
            for goal in value:
                if goal["measure"][1] in opt.OBJECTIVEMETHODS:
                    goal["measure"][1] = opt.OBJECTIVEMETHODS[goal["measure"][1]]

        opt.trackGoals = {}
        for key, value in opt.objectives.items():
            for goal in value:
                opt.trackGoals[f"indice {key}: {goal['measure'][0]} {goal['measure'][1].__name__}"] = []

        opt.variablesValues = [1.0] * len(opt.variablesToOptimize)
        opt.bounds = [(None, None)] * len(opt.variablesToOptimize)
        for var in start_point:
            idx = opt.variablesToOptimize.index(var)
            if "start" in start_point[var]:
                opt.variablesValues[idx] = start_point[var]["start"]
            if "bounds" in start_point[var]:
                opt.bounds[idx] = tuple(start_point[var]["bounds"])

        opt.plotMSE = []
        opt.plotIterate = []
        opt.trackVariables = []
        opt.iterationTrack = 0

        variable_names = req.variable_names
        emit_constraints = req.emit_constraints

        def evaluate_fn(params: Dict[str, float]) -> Dict[str, Any]:
            vals = [float(params[name]) for name in variable_names]
            mse = float(opt._optiSpeed(vals))
            # Prevent unbounded memory growth from tracking lists
            opt.trackGoals = {k: [] for k in opt.trackGoals}
            opt.trackVariables = []
            opt.plotMSE = []
            opt.plotIterate = []
            if not math.isfinite(mse):
                mse = 1e6
            result = {"objective_1": -mse}
            if emit_constraints:
                result["constraint_1"] = 0.0 if (math.isfinite(mse) and mse < 1e4) else 1.0
            return result

        async with state.lock:
            state.evaluate_fn = evaluate_fn
            state.variable_names = variable_names
            state.n_objectives = req.n_objectives
            state.n_constraints = req.n_constraints
            state.eval_count = 0

        logger.info("Glyfada evaluation context configured: %d variables, %d objectives",
                     len(variable_names), req.n_objectives)
        return {
            "status": "configured",
            "variable_names": variable_names,
            "n_objectives": req.n_objectives,
            "n_constraints": req.n_constraints,
        }

    except Exception as e:
        logger.error("Glyfada setup failed: %s", traceback.format_exc())
        raise HTTPException(status_code=400, detail=str(e))


_PICKLE_ALLOWED_DIR = os.environ.get("GLYFADA_PICKLE_DIR", tempfile.gettempdir())


@glyfada_router.post("/setup-pickle")
async def glyfada_setup_pickle(req: GlyfadaSetupFromPickleRequest):
    """Configure evaluation from an existing pickled objective (glyfada_eval.py format).

    SECURITY WARNING: pickle.load executes arbitrary code. Access is restricted
    to files under GLYFADA_PICKLE_DIR (default: system temp directory) to limit
    path traversal, but only trusted pickles should ever be placed there.
    """
    state: _GlyfadaState = glyfada_router.state

    if req.n_objectives != 1:
        raise HTTPException(
            status_code=400,
            detail=f"n_objectives must be 1 (got {req.n_objectives}). Multi-objective not yet supported.",
        )

    try:
        resolved = os.path.realpath(req.pickle_path)
        allowed = os.path.realpath(_PICKLE_ALLOWED_DIR)
        if not resolved.startswith(allowed + os.sep) and resolved != allowed:
            raise HTTPException(
                status_code=400,
                detail=f"Pickle path must be under the allowed directory ({_PICKLE_ALLOWED_DIR}).",
            )

        try:
            import cloudpickle as pickle
        except ImportError:
            import pickle

        with open(resolved, "rb") as f:
            pkl_state = pickle.load(f)

        obj_func = pkl_state["objective_func"]
        variable_names = pkl_state["variable_names"]
        emit_constraints = req.emit_constraints

        def evaluate_fn(params: Dict[str, float]) -> Dict[str, Any]:
            vals = [float(params[name]) for name in variable_names]
            try:
                mse = float(obj_func(vals))
            except Exception:
                mse = 1e6
            if not math.isfinite(mse):
                mse = 1e6
            result = {"objective_1": -mse}
            if emit_constraints:
                result["constraint_1"] = 0.0 if (math.isfinite(mse) and mse < 1e4) else 1.0
            return result

        async with state.lock:
            state.evaluate_fn = evaluate_fn
            state.variable_names = variable_names
            state.n_objectives = req.n_objectives
            state.n_constraints = req.n_constraints
            state.eval_count = 0

        logger.info("Glyfada evaluation context loaded from pickle: %s (%d variables)",
                     req.pickle_path, len(variable_names))
        return {
            "status": "configured",
            "variable_names": variable_names,
            "n_objectives": req.n_objectives,
            "n_constraints": req.n_constraints,
        }

    except Exception as e:
        logger.error("Glyfada setup-pickle failed: %s", traceback.format_exc())
        raise HTTPException(status_code=400, detail=str(e))


@glyfada_router.post("/evaluate")
async def glyfada_evaluate(req: GlyfadaEvalRequest):
    """Evaluate a parameter set. Called by glyfada's HttpEvaluator."""
    state: _GlyfadaState = glyfada_router.state

    # Snapshot state atomically to avoid race with concurrent /setup calls
    async with state.lock:
        evaluate_fn = state.evaluate_fn
        variable_names = state.variable_names

    if evaluate_fn is None:
        raise HTTPException(
            status_code=400,
            detail="Evaluation context not configured. Call /glyfada/setup first."
        )

    missing = [v for v in variable_names if v not in req.params]
    if missing:
        raise HTTPException(
            status_code=400,
            detail=f"Missing required parameters: {missing}",
        )

    try:
        def _run():
            with state.eval_lock:
                state.eval_count += 1
                return state.evaluate_fn(req.params)

        loop = asyncio.get_running_loop()
        result = await loop.run_in_executor(None, _run)
        return result

    except Exception as e:
        logger.error("Glyfada evaluation failed: %s", traceback.format_exc())
        raise HTTPException(status_code=500, detail="Evaluation failed. Check server logs for details.")


@glyfada_router.get("/status")
async def glyfada_status():
    """Current state of the evaluation endpoint."""
    state: _GlyfadaState = glyfada_router.state
    return {
        "configured": state.evaluate_fn is not None,
        "variable_names": state.variable_names,
        "n_objectives": state.n_objectives,
        "n_constraints": state.n_constraints,
        "eval_count": state.eval_count,
    }


# ---------------------------------------------------------------------------

app = FastAPI()
ebeam = beam()  # Shared instance is safe: gen_6d_gaussian/gen_6d_from_twiss are stateless

# Attach glyfada state and register router
glyfada_router.state = _GlyfadaState()
app.include_router(glyfada_router)
# Allow requests from your frontend (CORS!)
app.add_middleware(
    CORSMiddleware,
    allow_origins=ORIGINS,  # In production, use your frontend's exact origin
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

def getPngObjFromBeamList(beamlist, plotParams: PlottingParameters):
    beam_dist = None
    print(plotParams.beam_setup)
    if plotParams.beam_setup == 'twiss': beam_dist = ebeam.gen_6d_from_twiss(plotParams.twiss.model_dump(), plotParams.num_particles)
    else: beam_dist = ebeam.gen_6d_gaussian(0,[1,1,1,1,0.1,100], plotParams.num_particles)
    schem = draw_beamline()
    axList, lineAxObj = schem.plotBeamPositionTransform(beam_dist, beamlist, plot=False, apiCall=True, scatter=True, interval=plotParams.interval)
    fig = lineAxObj['axis'].figure
    buf = io.BytesIO()
    fig.savefig(buf, format="png",bbox_inches="tight")
    buf.seek(0)
    lineAx_img = base64.b64encode(buf.read()).decode("utf-8")
    buf.close()

    images = {}
    for index, axes in axList.items():

        fig = axes.figure
        buf = io.BytesIO()
        fig.savefig(buf, format="png",bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")
        buf.close()

        images.update({index: img_base64})
     
    lineAxObj['axis'] = lineAx_img
    print(lineAxObj['twiss'])
    lineAxObj['twiss'] = lineAxObj['twiss'].to_json()
    beamsegmentJson = []
    #for segment in lineAxObj['beamsegment']:
    #    beamsegmentJson.append(segment.__dict__)
    lineAxObj['beamsegment'] = beamsegmentJson
    #print(beamsegmentJson)

    lineAxObj = LineAxObject(**lineAxObj)
    pngObject = AxesPNGData(**{'images': images, 'line_graph': lineAxObj})
    plt.close('all')
    return pngObject

@app.get("/")
def root():
    return {"Hello" : "World!"}

@app.post("/get-dist")
def gen_beam(particle_num : int):
    particle_num = min(particle_num, 100_000)
    beam_dist = ebeam.gen_6d_gaussian(0,[1,1,1,1,0.1,100], particle_num).tolist()
    return beam_dist

@app.post("/excel-to-beamline")
def excelToBeamline(excelJson: list[Dict[str, Any]]) -> list[dict[str, dict[str, Any]]]:
    try:
        excelHandler = ExcelElements(excelJson)
        beamlist = excelHandler.create_beamline()

        jsonBeamlist = []

        for segment in beamlist:
            clas = segment.__class__
            className = clas.__name__
            classSig= inspect.signature(clas.__init__)

            paramsDict = {}
            for name, param in classSig.parameters.items():
                if name == "self":
                    continue
                paramVal = getattr(segment, name, None)
                paramsDict.update({name: paramVal})
                    
            jsonBeamlist.append({className: paramsDict})

        return jsonBeamlist
    except Exception as e:
        print(e)
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/axes")
def loadAxes(plotParams: PlottingParameters) -> AxesPNGData:
    try:
        plotParams.num_particles = min(plotParams.num_particles, MAX_PARTICLES)
        beamline = importlib.import_module("beamline")
        beamlist = []
        beamlineData = plotParams.beamlineData
        for segment in beamlineData:
            beamlist.append(_instantiate_beamline_class(beamline, segment.segmentName, segment.parameters))

        latObj = lattice(1)
        beamlist = latObj.changeBeamType(plotParams.beamType, plotParams.kineticE, beamlist)

        pngObject = getPngObjFromBeamList(beamlist, plotParams)
        return pngObject
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.get("/beamsegmentinfo")
def getBeamSegmentInfo():
    module = importlib.import_module(moduleName)
    classes = inspect.getmembers(module, inspect.isclass)
    classes_in_module = [cls for name, cls in classes if cls.__module__ == moduleName and cls.__name__ not in ["Beamline", "lattice"]]
    beamSegInfo = {}

    for cls in classes_in_module:
        sig = inspect.signature(cls.__init__)
        params_info = {}
    
        for name, param in sig.parameters.items():
            if name == "self":
                continue
            default = (
                param.default
                if param.default is not inspect.Parameter.empty
                else 1  # or some other marker
            )
            params_info[name] = default

        params_info['color'] = cls.color  # Manually add class info about beam's color
    
        beamSegInfo[cls.__name__] = params_info

    return beamSegInfo

@app.post("/plot-parameters")
def plot_parameters(graphParams: GraphParameters) -> List[GraphPlotData]:
    LABELMAPPING = {
        r'$\epsilon$ ($\pi$.mm.mrad)': 'emittance',
        r'$\alpha$': 'alpha',
        r'$\beta$ (m)': 'beta',
        r'$\gamma$ (rad/m)': 'gamma',
        r'$D$ (m)': 'dispersion',
        r'$D^{\prime}$': 'dispersion_prime',
        r'$\phi$ (deg)': 'angle',
        r'Envelope $E$ (mm)': 'envelope'
    }
    if graphParams.custom_step <= 0:
        raise HTTPException(status_code=400, detail="custom_step must be positive")
    try:
        beamline = importlib.import_module(moduleName)
        beamlist = []
        beamlineData = graphParams.beamline_data
        for segment in beamlineData:
            beamlist.append(_instantiate_beamline_class(beamline, segment.segmentName, segment.parameters))

        cleanedBeamlist = beamlist[:graphParams.beam_index]

        schem = draw_beamline()
        ebeam = beam()
        beam_dist = ebeam.gen_6d_gaussian(0,[1,1,1,1,0.1,100], 1000)

        # print("Plotting initial beamline up to segment", cleanedBeamlist)

        #  100 chosen as a large number to speed up initial calculation
        schem.plotBeamPositionTransform(beam_dist, cleanedBeamlist, plot=False, interval=100, rendering=False)
        beam_dist = schem.matrixVariables

        beamObj = Beamline(beamlist)
        indexOfSSegment = beamObj.findSegmentAtPos(graphParams.target_s_pos)

        newSegment = copy.deepcopy(beamObj.beamline[indexOfSSegment])
        newSegment.length = graphParams.target_s_pos - beamObj.beamline[indexOfSSegment - 1].endPos
        optimized_beamlist = beamObj.beamline[graphParams.beam_index:indexOfSSegment]
        optimized_beamlist.append(newSegment)

        # for i in optimized_beamlist:
        #     print("Printing segment:", i) 

        plotInfo = [] 
        domain_range = np.arange(graphParams.min, graphParams.max, graphParams.custom_step).tolist()
        if graphParams.max not in domain_range: domain_range.append(graphParams.max)

        for i in domain_range:
            setattr(optimized_beamlist[0], graphParams.target_parameter, i)
            twiss = schem.plotBeamPositionTransform(beam_dist, optimized_beamlist, plot=False, interval=100, rendering=False)
            # col = LABELMAPPING.get(graphParams.twiss_target, graphParams.twiss_target)
            plotDict = {f'parameter_value': i, 
                         'data': [
                                    {
                                        **{name: None if axis[-1] is None or math.isnan(axis[-1]) or math.isinf(axis[-1]) else axis[-1]
                                        for name, axis in twiss[col].items()},
                                        'twiss_parameter': LABELMAPPING.get(col, col)
                                    }
                                    for col in twiss.columns
                                 ]
                        }
            # print(plotDict)
            plotInfo.append(plotDict)

        # RETURN ALL TWISS, MAKE USER SELECT WHICH ONE
        return plotInfo
    
    except Exception as e:
        print(e)
        raise HTTPException(status_code=400, detail=str(e))

# Don't use, doesn't check for changes and server reloads
#if __name__ == "__main__":
    #uvicorn.run(app, host="127.0.0.1", port=8000)
