import numpy as np
import pandas as pd
import xarray as xr
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots

species = ["DOC", "NO3", "NH4"]

@st.cache_data()
def read_data(var:str):
    raw = pd.read_excel(".hiddendata/SampledData.xlsx")
    times = raw["Time (hr)"].to_numpy()/24.0
    results = raw[[f"{var}_{i}" for i in [*"IABCD"]]].to_numpy().T
    data = xr.DataArray(
        results, 
        dims=("z","t"), 
        coords={
            "z": np.array([0, 0.25, 0.65, 1.05, 1.45]), 
            "t": times})
    
    return data

all_data = {var:read_data(var) for var in species}

tabs = st.tabs([f"**{i}**" for i in species])

for tab, (specie, data) in zip(tabs, all_data.items()):
    with tab:
        fig = go.Figure()

        fig.add_trace(
            go.Heatmap(
                x=data.t, y=data.z, z=data,
                #zmin=18, zmax=22,
                hovertemplate=f"<b>{specie}"+" = %{z:.1f} mg/L</b> <br> z = %{y:.2f} m <br> t = %{x:.1f} d",
                colorscale="Reds",
                name="T",
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
            title=f"{specie} [mg/L]",
        )

        st.plotly_chart(fig, use_container_width=True)