import pyvista as pv
import plotly.graph_objects as go
import numpy as np
from myusefultools.parser import getVTKList
import os
import re

import streamlit as st

print("----------")
col1,col2 = st.columns([1,3])

REFERENCE_PERMEABILITY = 9.99e-11
PATH_TO_VTK = "./sandbox/VTK"
all_vtk_paths = [os.path.join(PATH_TO_VTK,f) for f in getVTKList(PATH_TO_VTK)]
pattern = re.compile(r'_(\d+).vtk')

chemOptions = [
  "perm",
  "BAP",
  "DOC",
  "EPS",
  "NH4",
  "NO3",
  "O2",
  "POCr",
  "XAR",
  "XDN",
  "XI",
  "XN",
  "h",
  "katt_BAP",
  "katt_POCr",
  "kdet_EPS",
  "kdet_XI",
  "n",
  "rDN",
  "rH",
  "rN",
  "U"
]

def convertPathToTimeLabel(f):
    t = pattern.search(f).group(1)
    return f"{float(t)/4320/2:.1f} day"

with col1:
    st.markdown("### Time:")
    selected_time = st.select_slider(" ",all_vtk_paths[1:],format_func=convertPathToTimeLabel)
    st.markdown("----")
    st.markdown("### Variable:")
    selected_chem = st.selectbox(" ",chemOptions)

mesh = pv.read(selected_time)
slice = mesh.slice(normal='y')
data = slice[selected_chem]

fig = go.Figure()
fig.update_layout(
    autosize=False,
    height=800,
    title={
        'text': selected_chem,
        'y': 0.9,
        'x': 0.5,
        'xanchor': 'center',
        'yanchor': 'top'},
    xaxis={
        'title':"X (m)",
        'exponentformat' : "power"},
    yaxis={
        'title':"Z (m)",
        'exponentformat' : "power",
        'scaleanchor':'x'},
    legend={
        'orientation':"h",
        'yanchor':"top",
        'y':-0.1,
        'xanchor':"left",
        'x':0.01,
        'title':None,
        'groupclick':"toggleitem"},
    font={
        'size': 14,
        'color': "darkblue"}
    )

if selected_chem == "perm": 
    with col1: log_scale = st.checkbox("Normalized log-scale?",False)
    if log_scale: data = np.log10(data/REFERENCE_PERMEABILITY)

if selected_chem == "U": 
    with col1: u_component = st.selectbox("Component",["Mag",0,1,2])
    if u_component == "Mag":
        data = np.linalg.norm(data,axis=1)  #Shows only the Z-component
    else:
        data = data.T[u_component]

# From blockMeshDict
lenght, width, depth = 0.441, 0.300, 0.010
deltaLenght, deltaWidth = 0.005, 0.005
LCells, WCells = round(lenght/deltaLenght), round(width/deltaWidth) + 1

fig.add_trace(go.Heatmap(
        dx = deltaWidth, dy = deltaLenght, z = data.reshape([LCells,WCells]),
        name=str(convertPathToTimeLabel(selected_time)) + " d",
        colorbar={'exponentformat':"power"}))

## Show in webpage
with col2: st.plotly_chart(fig)
