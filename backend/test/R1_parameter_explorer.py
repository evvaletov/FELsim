"""R1: Interactive Parameter Explorer

Plotly/Dash dashboard for exploring S4 (1D) and S5 (2D) parameter scan
results. Three tabs: 1D scans, 2D feasibility maps, S4/S5 comparison.

Usage:
    python R1_parameter_explorer.py [--port 8050]

Then open http://localhost:8050 in a browser.

Author: Eremey Valetov
"""

import sys
import csv
import math
import argparse
from pathlib import Path

import numpy as np

try:
    import dash
    from dash import dcc, html, Input, Output, State
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError:
    print("R1 requires: pip install dash plotly")
    sys.exit(1)

backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from UHM_beamline_opt_05ps_params import MSE_THRESHOLDS, QUAD_INDICES

S4_DIR = Path(__file__).resolve().parent / 'results' / 'params_05ps'
S5_DIR = Path(__file__).resolve().parent / 'results' / 'params_05ps_2d'

SCAN_1D = {
    'energy_spread': {
        'label': 'Energy Spread σ_E (%)',
        'file': 'scan_energy_spread.csv',
    },
    'chirp': {
        'label': 'Chirp h (×10⁹/s)',
        'file': 'scan_chirp.csv',
        'display_scale': 1e-9,
    },
    'emittance': {
        'label': 'Emittance ε_n (π·mm·mrad)',
        'file': 'scan_emittance.csv',
    },
}

SCAN_2D = {
    's5a': {
        'param1': 'energy_std_percent', 'param2': 'h',
        'p1_label': 'σ_E (%)', 'p2_label': 'h (×10⁹/s)',
        'p2_display_scale': 1e-9,
        'title': 'S5a: σ_E × h (ε_n=8)',
    },
    's5b': {
        'param1': 'energy_std_percent', 'param2': 'epsilon_n',
        'p1_label': 'σ_E (%)', 'p2_label': 'ε_n (π·mm·mrad)',
        'p2_display_scale': 1.0,
        'title': 'S5b: σ_E × ε_n (h=5e9)',
    },
    's5c': {
        'param1': 'h', 'param2': 'epsilon_n',
        'p1_label': 'h (×10⁹/s)', 'p2_label': 'ε_n (π·mm·mrad)',
        'p1_display_scale': 1e-9, 'p2_display_scale': 1.0,
        'title': 'S5c: h × ε_n (σ_E=0.5%)',
    },
}


def read_csv(filepath):
    if not filepath.exists():
        return []
    with open(filepath) as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def classify(mse):
    if math.isnan(mse):
        return 'Failed'
    for label in ['Excellent', 'Acceptable', 'Marginal']:
        if mse < MSE_THRESHOLDS[label]:
            return label
    return 'Failed'


QUALITY_COLORS = {
    'Excellent': '#2ecc71', 'Acceptable': '#f1c40f',
    'Marginal': '#e67e22', 'Failed': '#e74c3c',
}

# Display RMS = sqrt(MSE) in plots and hover text
RMS_THRESHOLDS = {k: math.sqrt(v) for k, v in MSE_THRESHOLDS.items()}


# ── Robust y-axis limits ─────────────────────────────────────────────────────

def _robust_range(values, pad_frac=0.15):
    """Shortest Covering Interval for linear-scale y-axis limits.

    Finds the shortest interval covering ≥80% of the data, extends by
    50% of its span on each side. Consistent with the SCI approach used
    in _robust_log_range and the static matplotlib plots.
    """
    clean = sorted(v for v in values if np.isfinite(v))
    if not clean:
        return [-1, 1]
    if len(clean) < 3:
        lo, hi = min(clean), max(clean)
        margin = max((hi - lo) * 0.3, 0.5)
        return [lo - margin, hi + margin]
    n = len(clean)
    k = max(math.ceil(0.8 * n), min(n, 4))
    best_span, best_lo, best_hi = np.inf, clean[0], clean[k - 1]
    for i in range(n - k + 1):
        span = clean[i + k - 1] - clean[i]
        if span < best_span:
            best_span = span
            best_lo, best_hi = clean[i], clean[i + k - 1]
    ext = max(best_span * 0.5, 0.1)
    lo, hi = best_lo - ext, best_hi + ext
    if hi - lo < 1e-10:
        mid = (lo + hi) / 2
        lo, hi = mid - 0.5, mid + 0.5
    pad = (hi - lo) * pad_frac
    return [lo - pad, hi + pad]


def _robust_log_range(values, pad_decades=0.5):
    """Determine log-scale y-axis range via Shortest Covering Interval.

    Finds the shortest interval in log10 space that covers ≥80% of the
    data, then extends by 50% of its span on each side.  This naturally
    identifies the densest region and excludes extreme stragglers without
    relying on IQR (which fails for N<15 with heavy tails).
    """
    clean = sorted(v for v in values if np.isfinite(v) and v > 0)
    if not clean:
        return [-10, 0]
    if len(clean) < 3:
        lo = np.log10(min(clean))
        hi = np.log10(max(clean))
        margin = max((hi - lo) * 0.3, 1.0)
        return [lo - margin, hi + margin]

    lv = np.log10(clean)
    n = len(lv)
    k = max(math.ceil(0.8 * n), min(n, 4))

    # Sliding window: find shortest interval covering k points
    best_span = np.inf
    best_lo = lv[0]
    best_hi = lv[k - 1]
    for i in range(n - k + 1):
        span = lv[i + k - 1] - lv[i]
        if span < best_span:
            best_span = span
            best_lo = lv[i]
            best_hi = lv[i + k - 1]

    # Extend core by 50% of its span on each side
    ext = max(best_span * 0.5, 0.5)  # minimum 0.5 decade extension
    lo = best_lo - ext
    hi = best_hi + ext

    # Enforce minimum 1-decade total span
    if hi - lo < 1.0:
        mid = (lo + hi) / 2
        lo, hi = mid - 0.5, mid + 0.5

    return [lo - pad_decades, hi + pad_decades]


# ── 1D scan figure ──────────────────────────────────────────────────────────

def make_1d_figure(scan_name):
    cfg = SCAN_1D[scan_name]
    rows = read_csv(S4_DIR / cfg['file'])
    if not rows:
        return go.Figure().add_annotation(text="No data", showarrow=False)

    scale = cfg.get('display_scale', 1.0)
    x = [r['param_value'] * scale for r in rows]
    rms = [math.sqrt(r['mse']) if r['mse'] >= 0 else float('nan') for r in rows]
    colors = [QUALITY_COLORS[classify(r['mse'])] for r in rows]

    # Hover text with quad currents
    hover = []
    for r in rows:
        rms_val = math.sqrt(r['mse']) if r['mse'] >= 0 else float('nan')
        lines = [f"RMS: {rms_val:.2e} ({classify(r['mse'])})"]
        lines.append(f"β_x: {r['beta_x']:.4f} m, β_y: {r['beta_y']:.4f} m")
        lines.append(f"α_x: {r['alpha_x']:.4f}, α_y: {r['alpha_y']:.4f}")
        lines.append("─ Stage 11 quads ─")
        for qi in [87, 93, 95, 97]:
            lines.append(f"  q{qi}: {r[f'quad_{qi}']:.4f} A")
        hover.append("<br>".join(lines))

    fig = make_subplots(rows=2, cols=1, row_heights=[0.6, 0.4],
                        subplot_titles=["RMS Twiss Mismatch", "Twiss Parameters"],
                        vertical_spacing=0.15)

    fig.add_trace(go.Scatter(
        x=x, y=rms, mode='markers+lines',
        marker=dict(color=colors, size=8, line=dict(width=1, color='black')),
        hovertext=hover, hoverinfo='text',
        name='RMS',
    ), row=1, col=1)

    for label, thresh in RMS_THRESHOLDS.items():
        fig.add_hline(y=thresh, line_dash='dot', line_color='gray',
                      annotation_text=label, row=1, col=1)

    rms_range = _robust_log_range(rms)
    fig.update_yaxes(type='log', title_text='RMS Twiss Mismatch', row=1, col=1,
                     range=rms_range)

    # Twiss panel — use robust y-limits excluding outliers
    all_twiss = []
    for key, name, color in [
        ('beta_x', 'β_x', '#3498db'), ('beta_y', 'β_y', '#e74c3c'),
        ('alpha_x', 'α_x', '#2ecc71'), ('alpha_y', 'α_y', '#f39c12'),
    ]:
        vals = [r[key] for r in rows]
        all_twiss.extend(vals)
        fig.add_trace(go.Scatter(
            x=x, y=vals, mode='markers+lines',
            name=name, marker=dict(size=5), line=dict(color=color),
        ), row=2, col=1)

    twiss_lo, twiss_hi = _robust_range(all_twiss)
    fig.update_xaxes(title_text=cfg['label'], row=2, col=1)
    fig.update_yaxes(title_text='Value', row=2, col=1,
                     range=[twiss_lo, twiss_hi])
    fig.update_layout(height=700, showlegend=True,
                      title=f"S4: {scan_name.replace('_', ' ').title()} Sweep")
    return fig


# ── 2D scan figure ──────────────────────────────────────────────────────────

def make_2d_figure(scan_name):
    cfg = SCAN_2D[scan_name]
    rows = read_csv(S5_DIR / f'scan_{scan_name}.csv')
    if not rows:
        return go.Figure().add_annotation(text="No data", showarrow=False)

    p1_key, p2_key = cfg['param1'], cfg['param2']
    p1_scale = cfg.get('p1_display_scale', 1.0)
    p2_scale = cfg.get('p2_display_scale', 1.0)

    p1_vals = sorted(set(r[p1_key] for r in rows))
    p2_vals = sorted(set(r[p2_key] for r in rows))

    mse_grid = np.full((len(p2_vals), len(p1_vals)), np.nan)
    qual_grid = np.full((len(p2_vals), len(p1_vals)), '', dtype=object)
    hover_grid = np.full((len(p2_vals), len(p1_vals)), '', dtype=object)

    p1_map = {v: i for i, v in enumerate(p1_vals)}
    p2_map = {v: i for i, v in enumerate(p2_vals)}

    for r in rows:
        i = p1_map.get(r[p1_key])
        j = p2_map.get(r[p2_key])
        if i is not None and j is not None:
            mse_grid[j, i] = r['mse']
            qual_grid[j, i] = classify(r['mse'])
            rms_val = math.sqrt(r['mse']) if r['mse'] >= 0 else float('nan')
            hover_grid[j, i] = (
                f"{cfg['p1_label']}: {r[p1_key] * p1_scale:.3g}<br>"
                f"{cfg['p2_label']}: {r[p2_key] * p2_scale:.3g}<br>"
                f"RMS: {rms_val:.2e}<br>"
                f"Quality: {classify(r['mse'])}"
            )

    p1_d = [v * p1_scale for v in p1_vals]
    p2_d = [v * p2_scale for v in p2_vals]

    rms_grid = np.where(np.isnan(mse_grid), 10.0,
                        np.sqrt(np.clip(mse_grid, 1e-20, 1e5)))

    fig = go.Figure(data=go.Heatmap(
        z=np.log10(rms_grid),
        x=p1_d, y=p2_d,
        text=hover_grid, hoverinfo='text',
        colorscale='Viridis', reversescale=False,
        colorbar=dict(title='log₁₀(RMS)'),
        zmin=-5, zmax=1,
    ))

    # Threshold contour lines (in RMS)
    for label, color in [
        ('Excellent', 'lime'), ('Acceptable', 'yellow'), ('Marginal', 'orange'),
    ]:
        rms_thresh = RMS_THRESHOLDS[label]
        fig.add_trace(go.Contour(
            z=np.log10(rms_grid), x=p1_d, y=p2_d,
            contours=dict(start=np.log10(rms_thresh),
                          end=np.log10(rms_thresh),
                          size=0, coloring='none'),
            line=dict(color=color, width=2, dash='dash'),
            showscale=False, name=label, hoverinfo='skip',
        ))

    fig.update_layout(
        title=cfg['title'],
        xaxis_title=cfg['p1_label'],
        yaxis_title=cfg['p2_label'],
        height=550,
    )
    return fig


# ── Dash app ────────────────────────────────────────────────────────────────

app = dash.Dash(__name__)
app.title = "FELsim Parameter Explorer"

app.layout = html.Div([
    html.H2("FELsim Parameter Explorer (R1)",
            style={'textAlign': 'center', 'marginBottom': '10px'}),

    dcc.Tabs(id='tabs', value='tab-1d', children=[
        dcc.Tab(label='S4: 1D Scans', value='tab-1d'),
        dcc.Tab(label='S5: 2D Maps', value='tab-2d'),
    ]),

    html.Div(id='tab-content', style={'padding': '20px'}),
])


@app.callback(Output('tab-content', 'children'), Input('tabs', 'value'))
def render_tab(tab):
    if tab == 'tab-1d':
        return html.Div([
            html.Div([
                html.Label("Select sweep:"),
                dcc.Dropdown(
                    id='scan-1d-selector',
                    options=[{'label': v['label'], 'value': k}
                             for k, v in SCAN_1D.items()],
                    value='emittance',
                    style={'width': '300px'},
                ),
            ], style={'marginBottom': '15px'}),
            dcc.Graph(id='graph-1d'),
        ])

    elif tab == 'tab-2d':
        return html.Div([
            html.Div([
                html.Label("Select 2D scan:"),
                dcc.Dropdown(
                    id='scan-2d-selector',
                    options=[{'label': v['title'], 'value': k}
                             for k, v in SCAN_2D.items()],
                    value='s5b',
                    style={'width': '400px'},
                ),
            ], style={'marginBottom': '15px'}),
            dcc.Graph(id='graph-2d'),
        ])


# Separate callbacks for each tab's graph to avoid missing component errors
@app.callback(
    Output('graph-1d', 'figure'),
    Input('scan-1d-selector', 'value'),
    prevent_initial_call=False,
)
def update_1d(scan_name):
    return make_1d_figure(scan_name or 'emittance')


@app.callback(
    Output('graph-2d', 'figure'),
    Input('scan-2d-selector', 'value'),
    prevent_initial_call=False,
)
def update_2d(scan_name):
    return make_2d_figure(scan_name or 's5b')


def main():
    parser = argparse.ArgumentParser(description='FELsim Parameter Explorer')
    parser.add_argument('--port', type=int, default=8050)
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()

    print(f"Starting R1 Parameter Explorer on http://localhost:{args.port}")
    print(f"  S4 data: {S4_DIR}")
    print(f"  S5 data: {S5_DIR}")
    app.run(host='0.0.0.0', port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()
