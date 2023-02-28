import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots

"# 🌵 Water content (θ)"

## Globals
colors = st.session_state.colors
z = np.array([25.0, 65.0, 105.0, 145.0])/100.
probes = [*"ABCD"]
add_flood_and_dry_periods = st.session_state["add_flood_and_dry_periods"]

@st.cache_data
def read_data():
    vwc = pd.read_excel(
        ".hiddendata/ORP.xlsx", 
        usecols=["Time (min)"] + [f"VWC_{i}" for i in probes]
    )
    vwc["Time (d)"] = vwc["Time (min)"]/(60*24)
    return vwc

vwc = read_data()

## Dataframe render
with st.expander("🗃️ **Dataframes**"):
    st.dataframe(vwc, height=200, use_container_width=True)

## Timeseries render
with st.expander("📈 **Timeseries**"):

    fig = go.Figure()
    add_flood_and_dry_periods(fig)

    for probe, color in zip(probes, colors):

        # Raw data
        fig.add_trace(
            go.Scatter(
                x=vwc["Time (d)"],
                y=vwc[f"VWC_{probe}"],
                line=dict(color=color, width=3),
                name=f"Probe {probe}"
            )
        )

    fig.update_xaxes(title_text="Time (d)")
    fig.update_yaxes(title_text="Vol. water content (θ) [%]", range=[5,44])

    fig.update_layout(
        title="Volumetric water content (θ)",
        legend=dict(
            orientation="h",
            yanchor="bottom", y=0.99,
            xanchor="center", x=0.50
        )
    )

    st.plotly_chart(fig, use_container_width=True)

## Heatmap evolution

with st.expander("🌽 **Spatial distribution**", expanded=True):
    t = vwc["Time (d)"].to_numpy()
    orpnp = vwc[[f"VWC_{i}" for i in probes]].to_numpy()
    
    fig = go.Figure()

    fig.add_trace(
        go.Heatmap(
            x=t, y=z, z=orpnp.T,
            hovertemplate="<b>θ = %{z:.1f}%</b> <br> z = %{y:.2f} m <br> t = %{x:.1f} d",
            colorscale="Blues",
            name="θ",
            colorbar=dict(
                title_text="%",
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
        title="Volumetric water content (θ)",
    )

    st.plotly_chart(fig, use_container_width=True)