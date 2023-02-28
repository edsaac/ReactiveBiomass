import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots

"# 🌡️ Temperature [°C]"

## Globals
colors = st.session_state.colors
z = np.array([25.0, 65.0, 105.0, 145.0])/100.
probes = [*"ABCD"]
add_flood_and_dry_periods = st.session_state["add_flood_and_dry_periods"]

@st.cache_data
def read_data():
    temp = pd.read_excel(
        ".hiddendata/ORP.xlsx", 
        usecols=["Time (min)"] + [f"Temp_{i}" for i in probes]
    )
    temp["Time (d)"] = temp["Time (min)"]/(60*24)
    return temp

temp = read_data()

## Dataframe render
with st.expander("🗃️ **Dataframes**"):
    st.dataframe(temp, height=200, use_container_width=True)

## Timeseries render
with st.expander("📈 **Timeseries**"):

    fig = go.Figure()
    add_flood_and_dry_periods(fig)
    for probe, color in zip(probes, colors):

        # Raw data
        fig.add_trace(
            go.Scatter(
                x=temp["Time (d)"],
                y=temp[f"Temp_{probe}"],
                line=dict(color=color, width=3),
                name=f"Probe {probe}"
            )
        )

    fig.update_xaxes(title_text="Time (d)")
    fig.update_yaxes(title_text="Vol. water content (θ) [%]")

    fig.update_layout(
        title="Temperature [°C]",
        legend=dict(
            orientation="h",
            yanchor="bottom", y=0.99,
            xanchor="center", x=0.50
        )
    )

    st.plotly_chart(fig, use_container_width=True)

## Heatmap evolution

with st.expander("🌽 **Spatial distribution**", expanded=True):
    t = temp["Time (d)"].to_numpy()
    tempnp = temp[[f"Temp_{i}" for i in probes]].to_numpy()
    
    fig = go.Figure()

    fig.add_trace(
        go.Heatmap(
            x=t, y=z, z=tempnp.T,
            zmin=18, zmax=22,
            hovertemplate="<b>T = %{z:.1f} °C</b> <br> z = %{y:.2f} m <br> t = %{x:.1f} d",
            colorscale="Reds",
            name="T",
            colorbar=dict(
                title_text="C",
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
        title="Temperature [°C]",
    )

    st.plotly_chart(fig, use_container_width=True)