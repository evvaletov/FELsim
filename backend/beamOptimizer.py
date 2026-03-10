#Author: Christian Komo

'''
helpful resources
https://www.youtube.com/watch?v=G0yP_TM-oag
https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
https://docs.scipy.org/doc/scipy/reference/optimize.html
gpt"can you xplain the parameters of a constraint paramter in scipy.minimize"
'''
#  NOTE: nelder-mead method doesn't work if starting search point is zero (if variablesValues = [0,0,0...]
#  NOTE: After testing, COBYLA and some other methods may try a test value of 0 for current, length, etc, that may throw back a difference of NAN
#        (because beamline object divides by 0). Have to figure out how to bound x variable values automatically so computer doesn't use weird values like 0, 
#        negative numbers (YOU MUST USE bounds for now so program doesn't use negative/zero numbers
#  NOTE: "measure" in self.objectives has to be a function call that returns a single value with a parameter of a
#         2d list of particles and each indice can only appear once as a key

import copy
import scipy.optimize as spo
from beamline import *
from ebeam import beam
import numpy as np
from schematic import *
import timeit
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import time

# NOTE: EACH BEAMOPTIMIZER OBJECT AFTER INSTANTIATION SHOULD ONLY BE USED TO 
# RUN A CALC() FUNCTION ONE TIME, CODE HAS NOT BEEN MODIFIED YET BEYOND ONE CALC() FUNCTION CALL

# For each indice of the beam segment in parameter, no proper error handling yet for invalid values, repeating indices with same objective, have to test...

#  Currently can only optimize one x variable for each segment indice, should we implement more than one variable in the future?

#   Do we begin plotting iterations at 0 or 1???

class beamOptimizer():
    def __init__(self, beamline, matrixVariables):
        '''
        Constructor for beamline and particle values to optimize over for given y objectives and x variables

        Parameters
        ----------
        beamline: list[beamline]
            list of beamline objects representing accelerator beam
        matrixVariables: np.array(list[float][float])
            2D numPy list of particle elements
        '''
        ebeam = beam()
        #  Methods included in class to return staistical information, add to dictionary if additional methods are created
        self.OBJECTIVEMETHODS = {"std": ebeam.std,"epsilon":ebeam.epsilon,"alpha":ebeam.alpha,
                                 "beta":ebeam.beta,"gamma":ebeam.gamma,"phi":ebeam.phi,
                                 "envelope": ebeam.envelope, "dispersion": ebeam.disper}
        
        self.matrixVariables = matrixVariables
        self.beamline = beamline
        self.use_apertures = False
        self.min_surviving_fraction = 0.1  # Return penalty if <10% survive
        self.transmission_weight = 0.0
    
    def _optiSpeed(self, variableVals):
        '''
        Simulates particle movement through a beamline, calculates y objective accuracy statistic for 
        particles after simulating with new x variables, and returns how accurate algorithm got by 
        using mean squared error stat

        Parameters
        ----------
        variableVals: list[float]
            test x value(s) to optimize, indice of each value corresponds to indice of variable in 
            variablesToOptimize

        Returns
        -------
        difference: float
            mean squared error statistic of how accurate all test objective values are to their goal
        '''
        particles = self.matrixVariables.copy()
        segments = self.beamline
        mse = []
        numGoals = 0
        
        n_initial = len(particles)

        #  Loop through beamline indices
        for i in range(len(segments)):
            #  Check if indice is in segmentVar (x variable dictionary)
            if i in self.segmentVar:
                try:
                    #  Adjust the segment's x variable value according to its mathematical relationship
                    yFunc = self.segmentVar.get(i)[2]
                    varIndex = self.variablesToOptimize.index(self.segmentVar.get(i)[0]) #  Get the index of the x variable to use with
                    newValue = yFunc(variableVals[varIndex])
                    param = self.segmentVar.get(i)[1]
                    particles = segments[i].useMatrice(particles, **{param: newValue})
                except TypeError as e:
                    raise ValueError(f"segment {i} has no parameter {param}")
            else:
                particles = segments[i].useMatrice(particles)
            if self.use_apertures:
                particles = segments[i].apply_aperture(particles)
                if len(particles) < 2:
                    # Can't compute Twiss with <2 particles — steep penalty
                    penalty = 1e6 * (1 - len(particles) / n_initial)
                    self.trackVariables.append(variableVals)
                    self.plotMSE.append(penalty)
                    self.plotIterate.append((self.iterationTrack) + 1)
                    self.iterationTrack += 1
                    return penalty
            #  Check if indice in objective dictionary
            if i in self.objectives:
                for goalDict in self.objectives[i]:
                    # add sum piece to calculate statistical accuracy with MSE
                    stat = (goalDict["measure"][1](particles, goalDict["measure"][0]))
                    goalDict["measured"] = stat
                    mse.append(((stat-goalDict["goal"])**2)*goalDict["weight"])
                    numGoals = numGoals+1
                    stringForm = "indice " + str(i) + ": " + goalDict["measure"][0] + " " + goalDict["measure"][1].__name__
                    self.trackGoals[stringForm].append(stat) #  for plotting in calc()

        # Transmission penalty: w_T × (1-T)²
        if self.use_apertures and self.transmission_weight > 0:
            T = len(particles) / n_initial
            mse.append(self.transmission_weight * (1 - T) ** 2)
            numGoals += 1
            self.trackGoals.setdefault("transmission", []).append(T)

        #  Calculate MSE
        if numGoals == 0:
            return np.inf
        difference = (np.sum(mse))/numGoals
        if not np.isfinite(difference):
            difference = 1e6

        #  For plotting purposes in calc()
        self.trackVariables.append(variableVals)
        self.plotMSE.append(difference)
        self.plotIterate.append((self.iterationTrack) + 1)
        self.iterationTrack = self.iterationTrack + 1

        return difference
    
    def calc(self, method, segmentVar, startPoint, objectives, plotProgress = False, plotBeam = False, printResults = False, use_apertures = False, transmission_weight = 0.0, **kwargs):
        '''
        optimizes beamline segment attribute values so y values are close to objective values as possible.
        Post optimization plotting supported

        Parameters
        ----------
        method: str
            name of minimization method/algorithm
        segmentVar: dict
            dictionary, each key is an indice corresponding to its value of a list of 
            x variable parameters 
        objectives: dict
            dictionary, each key is an indice corresponding to its value of a list of
            y objectives. In each list are dictionaries corresponding to the parameters of that 
            y objective 
        startPoint: dict
            dictionary, each key is an x variable corresponding to another dictionary of
            that variables' bounds, starting search point, and other parameters
        plotProgress: bool
            plot x variable and y objective values as a function of iterations
        plotBeam: bool
            plot beamline simulation with new x variables post-optimization
        printResults: bool
            output data in terminal

        Returns
        -------
        result: OptimizeResult
            Object containing resulting information about optimization process and results
        '''
        self.use_apertures = use_apertures
        self.transmission_weight = transmission_weight

        #  Variables for plotting purposes later
        self.plotMSE = []
        self.plotIterate = []
        self.trackVariables = []
        self.iterationTrack = 0
        self.trackGoals = {}

        #  Initialize ordered list of unique x variables to optimize
        self.segmentVar = segmentVar
        seen = {}
        for indice in self.segmentVar:
            if (indice < 0 or indice >= len(self.beamline)):
                raise IndexError(str(indice) + " is out of bounds for segmentVar dictionary")
            varName = self.segmentVar[indice][0]
            if varName not in seen:
                seen[varName] = True
        self.variablesToOptimize = list(seen.keys())

        #  Initialize objectives dictionary with measurement methods
        self.objectives = copy.deepcopy(objectives)
        for key, value in self.objectives.items():
            if (key not in range(len(self.beamline))):
                raise TypeError("Invalid indice: indice " + str(key) + " in objectives dict is out of bounds" )
            for goal in value:
                if goal["measure"][1] in self.OBJECTIVEMETHODS:
                    goal["measure"][1] = self.OBJECTIVEMETHODS[goal["measure"][1]]
                elif isinstance(goal["measure"][1], str):
                    raise TypeError("Invalid method name: No such method name exists in OBJECTIVESMETHOD dict")
                #  Used to keep track of data plotting through optimization
                #  Very rudementary, since looking for plotting of an objective relies on finding the same string name. Will have to improve in future
                self.trackGoals.update({"indice " + str(key) + ": " + goal["measure"][0] + " "  + goal["measure"][1].__name__: []})

        #  Create x variables' bounds and  start point list. 
        #  Order corresponding to x variable order in variablesToOptimize
        self.variablesValues = [] 
        self.bounds = []
        for i in self.variablesToOptimize:
            self.variablesValues.append(1) 
            self.bounds.append((None, None))
        for var in startPoint:
            index = self.variablesToOptimize.index(var)
            if "start" in startPoint.get(var): self.variablesValues[index] = startPoint.get(var).get("start")
            if "bounds" in startPoint.get(var): self.bounds[index] = startPoint.get(var).get("bounds")

        # Time speed to minimize difference of objective function
        startTime = time.perf_counter()
        if method == 'glyfada':
            from glyfadaAdapter import GlyfadaOptimizer

            CORE_KEYS = (
                'pop_size', 'max_gen', 'sigma', 'n_processes',
                'timeout_minutes', 'algorithm', 'debug',
                'n_objectives', 'constraints', 'seed_from_results', 'callback',
            )
            glyfada_kwargs = {k: kwargs[k] for k in CORE_KEYS if k in kwargs}

            # Everything else goes to extra_config
            extra_config = {k: v for k, v in kwargs.items()
                           if k not in CORE_KEYS and k not in
                           ('plotProgress', 'plotBeam', 'printResults')}
            if extra_config:
                glyfada_kwargs['extra_config'] = extra_config

            optimizer = GlyfadaOptimizer(
                objective_func=self._optiSpeed,
                variable_names=self.variablesToOptimize,
                bounds=self.bounds,
                default_values=self.variablesValues,
                **glyfada_kwargs
            )
            result = optimizer.optimize()
        else:
            result = spo.minimize(self._optiSpeed, self.variablesValues, method=method, bounds=self.bounds)
        endTime = time.perf_counter()

        if result.x is not None:
            # print out new values for each beam segment's attribute
            output = "\nx variables:"
            for indice in self.segmentVar:
                    variable = self.segmentVar.get(indice)[0]
                    index = self.variablesToOptimize.index(variable)
                    yFunc = self.segmentVar.get(indice)[2]
                    newVal = yFunc(result.x[index])
                    segAttr = self.segmentVar.get(indice)[1]
                    setattr(self.beamline[indice], segAttr, newVal)
                    if printResults:
                        output += "\nindice " + str(indice) + " new " + segAttr + " value: " + str(newVal)
            if printResults:
                output += "\n\ny objectives:\n"
                for indice, value in self.objectives.items():
                    for obj in value:
                        output += "indice " + str(indice) + ": " + obj["measure"][0] + " "  + obj["measure"][1].__name__ + " value of " + str(obj["measured"]) + "\n"
                if self.use_apertures and "transmission" in self.trackGoals and self.trackGoals["transmission"]:
                    output += f"Transmission: {self.trackGoals['transmission'][-1]*100:.1f}%\n"
                output += "Final difference: " + str(result.fun) + "\n"
                output += "\nTotal time: " + str(endTime-startTime) + " s\n"
                output +="Total iterations: " + str(self.iterationTrack) + "\n"
                print(output)
        else:
            # Multi-objective: Pareto front available in result.pareto_front
            if printResults:
                print("Multi-objective optimization complete. Pareto front available in result.")
                print(f"Total time: {endTime-startTime} s")
                print(f"Total iterations: {self.iterationTrack}")

        # Plot the progress of y objectives and x variables as a function of iterations
        if plotProgress:
            fig, ax = plt.subplots(2,1)
            handles = []

            # plot MSE line
            mseLine, =ax[1].plot(self.plotIterate, self.plotMSE, label = 'Mean Squared Error', color = 'black')
            ax[1].set_xlabel('Iterations')
            ax[1].set_yscale('log')
            ax[1].set_ylabel('Mean Squared Error')
            ax[1].set_title("MSE and objectives vs Iterations")
            ax[1].tick_params(axis='y')
            handles.append(mseLine)

            # Plot y goals
            ax2 = ax[1].twinx()
            mini = 0
            for i, key in enumerate(self.trackGoals):
                valLine, = ax2.plot(self.plotIterate, self.trackGoals[key], label = key)
                handles.append(valLine)
                tempMin = abs(min(self.trackGoals[key]))
                if i == 0 or mini>tempMin:
                    mini = tempMin
            ax2.set_yscale('symlog', linthresh=10**(np.ceil(np.log10(mini))))
            ax2.set_ylabel('Objective functions')
            ax2.legend(handles = handles, loc = 'upper right')

            # Plot x variables + sec/iteration
            tempTrackVari = np.array(self.trackVariables)
            handles = []
            for i in range(len(tempTrackVari[0])):
                varLine, = ax[0].plot(self.plotIterate, tempTrackVari[:,i], label = self.variablesToOptimize[i])
                handles.append(varLine)
            ax[0].set_xlabel('Iterations')
            ax[0].set_ylabel('Variable values')
            ax[0].set_title('Variable Values through each Iteration')
            timeLine = mlines.Line2D([], [], color = 'white', label=f'{round((endTime-startTime)/max(self.iterationTrack, 1),4)} s/iteration')
            handles.append(timeLine)
            ax[0].legend(handles = handles, loc = 'upper right')

           
            plt.tight_layout()
            plt.show()

        # Plot beam simulation with new values
        if plotBeam:
            schem = draw_beamline()
            tempPart = self.matrixVariables.copy()
            schem.plotBeamPositionTransform(tempPart, self.beamline)

        return result
   
