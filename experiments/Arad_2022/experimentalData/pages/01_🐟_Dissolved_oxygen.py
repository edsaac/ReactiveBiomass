import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots

"# 🐟 Dissolved oxygen (DO)"

## Globals
colors = st.session_state.colors
z = np.array([25.0, 65.0, 105.0, 145.0])/100.
probes = [*"ABCD"]
add_flood_and_dry_periods = st.session_state["add_flood_and_dry_periods"]

@st.cache_data
def read_data():
    dox = pd.read_excel(".hiddendata/DissolvedOxygen.xlsx")
    dox["Time (s)"] = (dox["DateTime"] - dox["DateTime"].iloc[0]).apply(lambda x: x.total_seconds())
    return dox

dox = read_data()

## Dataframe render
with st.expander("🗃️ **Dataframes**"):
    st.dataframe(dox, height=200, use_container_width=True)

## Timeseries render
with st.expander("📈 **Timeseries**"):

    fig = go.Figure()
    add_flood_and_dry_periods(fig)

    for probe, color in zip(probes, colors):

        # Raw data
        fig.add_trace(
            go.Scatter(
                x=dox["Time (s)"]/(60*60*24), 
                y=dox[f"DO_{probe}"],
                line=dict(color=color, width=3),
                name=f"Probe {probe}",
                hovertemplate="<b>DO = %{y:.1f} mg/L</b><br> t = %{x:.2f} d",
            )
        )

    fig.update_xaxes(title_text="Time (d)")
    fig.update_yaxes(title_text="Dissolved Oxygen [mg/L]", range=[0, 9.0])

    fig.update_layout(
        title="Dissolved Oxygen",
        legend=dict(
            orientation="h",
            yanchor="bottom", y=1.03,
            xanchor="center", x=0.50
        )
    )

    st.plotly_chart(fig, use_container_width=True)

## Heatmap evolution

with st.expander("🌽 **Spatial distribution**", expanded=True):
    t = dox["Time (s)"].to_numpy()/86400.
    doxnp = dox[[f"DO_{i}" for i in probes]].to_numpy()
    
    fig = go.Figure()

    fig.add_trace(
        go.Heatmap(
            x=t, y=z, z=doxnp.T,
            zmin=0.0, zmax=9.0,
            hovertemplate="<b>DO = %{z:.1f} mg/L</b> <br> z = %{y:.2f} m <br> t = %{x:.1f} d",
            colorscale="Blues",
            name="DO",
            colorbar=dict(
                title_text="mg/L",
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
        title="Dissolved Oxygen"
    )

    st.plotly_chart(fig, use_container_width=True)