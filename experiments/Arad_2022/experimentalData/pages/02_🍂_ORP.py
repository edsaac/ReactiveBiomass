import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots

"# 🍂 Oxidation-reduction potential (ORP)"

## Globals
colors = st.session_state.colors
z = np.array([25.0, 65.0, 105.0, 145.0])/100.
probes = [*"ABCD"]
add_flood_and_dry_periods = st.session_state["add_flood_and_dry_periods"]

@st.cache_data
def read_data():
    orp = pd.read_excel(
        ".hiddendata/ORP.xlsx", 
        usecols=["Time (min)"] + [f"ORP_{i}" for i in probes]
    )
    orp["Time (d)"] = orp["Time (min)"]/(60*24)
    return orp

orp = read_data()

## Dataframe render
with st.expander("🗃️ **Dataframes**"):
    st.dataframe(orp, height=200, use_container_width=True)

## Timeseries render
with st.expander("📈 **Timeseries**"):

    fig = go.Figure()
    add_flood_and_dry_periods(fig)
    for probe, color in zip(probes, colors):

        # Raw data
        fig.add_trace(
            go.Scatter(
                x=orp["Time (d)"],
                y=orp[f"ORP_{probe}"],
                line=dict(color=color, width=3),
                name=f"Probe {probe}"
            )
        )

    fig.update_xaxes(title_text="Time (d)")
    fig.update_yaxes(title_text="ORP [mV]")

    fig.update_layout(
        title="Oxidation Reduction Potential (ORP)",
        legend=dict(
            orientation="h",
            yanchor="bottom", y=0.99,
            xanchor="center", x=0.50
        )
    )

    st.plotly_chart(fig, use_container_width=True)

## Heatmap evolution

with st.expander("🌽 **Spatial distribution**", expanded=True):
    t = orp["Time (d)"].to_numpy()
    orpnp = orp[[f"ORP_{i}" for i in probes]].to_numpy()
    
    fig = go.Figure()

    fig.add_trace(
        go.Heatmap(
            x=t, y=z, z=orpnp.T,
            hovertemplate="<b>ORP = %{z:.1f} mV</b> <br> z = %{y:.2f} m <br> t = %{x:.1f} d",
            colorscale="Bluered",
            name="ORP",
            colorbar=dict(
                title_text="mV",
                **st.session_state.colorbar_kwargs
            )
        ),
    )
    
    fig.update_xaxes(
        title_text="Time [d]", 
        showgrid=False,
        **st.session_state.axis_setup)

    fig.update_yaxes(
        title_text="Depth [m]", 
        range=[1.5, 0.0], 
        showgrid=False,
        **st.session_state.axis_setup)
    
    fig.update_layout(
        title="Oxidation Reduction Potential (ORP)",
    )

    st.plotly_chart(fig, use_container_width=True)