"""
Interactive dashboard (skeleton).

To keep dependencies minimal, this module provides a stub interface. For a full
Dash-based app, install dash and implement run_dashboard() accordingly.
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd


def run_dashboard(results_csv: Path, host: str = "127.0.0.1", port: int = 8050):
    try:
        import dash  # type: ignore
        from dash import html, dcc
        import plotly.express as px
    except Exception as e:
        raise RuntimeError("Dash is not installed. Install 'dash' to enable the dashboard.") from e

    df = pd.read_csv(results_csv)
    app = dash.Dash(__name__)
    # Controls
    targets = sorted([t for t in df.get("target_uniprot", pd.Series(dtype=str)).dropna().unique()])
    compounds = sorted([c for c in df.get("compound_id", pd.Series(dtype=str)).dropna().unique()])

    app.layout = html.Div([
        html.H2("MCP Drug Discovery Dashboard"),
        html.Div([
            html.Div([
                html.Label("Targets"),
                dcc.Dropdown(options=[{"label": t, "value": t} for t in targets], multi=True, id="dd-targets")
            ], style={"width": "48%", "display": "inline-block", "verticalAlign": "top"}),
            html.Div([
                html.Label("Compounds"),
                dcc.Dropdown(options=[{"label": c, "value": c} for c in compounds], multi=True, id="dd-compounds")
            ], style={"width": "48%", "display": "inline-block", "marginLeft": "4%", "verticalAlign": "top"}),
        ]),
        html.Div([
            dcc.Graph(id="hist-fused"),
            dcc.Graph(id="scatter-sim-score"),
            dcc.Graph(id="hist-evidence"),
        ])
    ])

    @app.callback(
        dash.dependencies.Output("hist-fused", "figure"),
        dash.dependencies.Output("scatter-sim-score", "figure"),
        dash.dependencies.Output("hist-evidence", "figure"),
        dash.dependencies.Input("dd-targets", "value"),
        dash.dependencies.Input("dd-compounds", "value"),
    )
    def update_figs(sel_targets, sel_compounds):
        dff = df.copy()
        if sel_targets:
            dff = dff[dff["target_uniprot"].isin(sel_targets)]
        if sel_compounds:
            dff = dff[dff["compound_id"].isin(sel_compounds)]
        # Hist fused
        if "fused_score" in dff.columns:
            fig1 = px.histogram(dff, x="fused_score", nbins=30, title="Fused Score Distribution")
        else:
            fig1 = px.scatter(x=[0], y=[0], title="No fused_score column")
        # Scatter similarity vs fused
        if "cosine_max" in dff.columns and "fused_score" in dff.columns:
            fig2 = px.scatter(dff, x="cosine_max", y="fused_score", title="Similarity vs Fused Score")
        else:
            fig2 = px.scatter(x=[0], y=[0], title="Missing columns: cosine_max/fused_score")
        # Evidence strength hist
        if "evidence_strength" in dff.columns:
            fig3 = px.histogram(dff, x="evidence_strength", nbins=20, title="Evidence Strength")
        else:
            fig3 = px.scatter(x=[0], y=[0], title="No evidence_strength column")
        return fig1, fig2, fig3

    app.run_server(debug=False, host=host, port=port)
