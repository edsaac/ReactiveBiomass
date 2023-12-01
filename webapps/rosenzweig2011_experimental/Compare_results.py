import streamlit as st
from pathlib import Path
from dataclasses import dataclass
from myusefultools.parser import vtk_to_xarray

import numpy as np
import xarray as xr
import pandas as pd
import streamlit as st

import plotly.graph_objects as go
from plotly.graph_objs.layout import XAxis

import git
import json

st.set_page_config(initial_sidebar_state="expanded")

st.session_state["css_style"] = """
    <style>
        h1{
            font-size: 5rem;
            text-align:center;
            text-opacity: 50%;
        }

        [data-testid="stSidebar"] ul li:first-child{
            font-size: 1.2rem;
            font-style: italic;
        }
    </style>
"""
st.markdown(st.session_state.css_style, unsafe_allow_html=True)

st.session_state["axis_setup"] = dict(
    mirror=False,
    ticks="outside",
    color="#bbb",
    showline=True)

casesList = ["inoculated"]
casename = st.radio("Select simulation:", casesList)
st.session_state.casepath = Path("fit_ravid", casename)

@dataclass
class ModelResult:
    """Class for organizing data"""
    var_name: str
    equation_repr: str = ""
    long_name: str = ""
    description: str = ""
    units: str = ""
    data: xr.DataArray|None = None
    is_vtk: bool = True
    is_logscale: bool = True
    zmin: float|None = None
    zmax: float|None = None
    colormap: str = "cividis"
    figure: go.Figure|None = None

    def __post_init__(self):
        if self.is_vtk:
            self.data = vtk_to_xarray(field=self.var_name, casepath=st.session_state.casepath)

    def plot_heatmap(self):
        
        if self.is_logscale:
            scalar = np.log10(self.data)
            # hovertemplate = f"log({self.equation_repr})" + " = %{z:.2e} "
            colorbar_title_text = f"<b>log({self.equation_repr})</b><br><i>log[{self.units}]</i>"       
        
        else: 
            scalar = self.data
            colorbar_title_text = f"<b>{self.equation_repr}</b><br><i>[{self.units}]</i>"

        hovertemplate = "<b>" + self.equation_repr + " = %{customdata:.2e} " + self.units + "</b><br> t = %{x:.2f} d <br> z = %{y:.1f} m" 

        fig = go.Figure()

        fig.add_trace(
            go.Heatmap(
                x=self.data.t/(60*60*24), y=self.data.z, z=scalar,
                customdata=self.data.to_numpy(),
                hovertemplate=hovertemplate,
                zmin=self.zmin, zmax=self.zmax,
                name=self.long_name,
                colorscale=self.colormap,
                colorbar=dict(
                    len=0.7,
                    lenmode="fraction",
                    thickness=0.02,
                    thicknessmode="fraction",
                    title_text=colorbar_title_text,
                    titlefont_size=18,
                    outlinewidth=0,
                    y=0.50,
                    x=1.05
                    )
            )
        )

        fig.update_xaxes(
            title=dict(
                text="Time [d]",
                font_color="#444"
                ), 
            showgrid=False,
            range=[0, 22.9],
            **st.session_state.axis_setup)

        fig.update_yaxes(
            title=dict(
                text="Depth [m]",
                font_color="#444"),
            range=[0.00, 0.65], 
            showgrid=True,
            **st.session_state.axis_setup)

        fig.update_layout(
            margin=dict(t=10, pad=5)
        )
        self.figure = fig
        return None

st.session_state.ModelResult = ModelResult

@st.cache_data
def get_results_blob(group:str):
    GIT_REPO_PATH = git.Repo('.', search_parent_directories=True).working_tree_dir
    with open(f"{GIT_REPO_PATH}/misc/data.json") as f:
        details = json.load(f)
    details = details[group]

    return { k:ModelResult(k, **v) for k,v in details.items() }

st.session_state.get_results_blob = get_results_blob

if st.button("Clear All"):
    st.cache_data.clear()