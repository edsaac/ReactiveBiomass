import pyvista as pv
from matplotlib.cm import get_cmap
from matplotlib.colors import to_hex

import plotly.graph_objects as go
import numpy as np
from myusefultools.parser import getVTKList
import os
from itertools import cycle

import streamlit as st

PATH_TO_VTK = "./column1/VTK"
all_vtk_paths = [os.path.join(PATH_TO_VTK,f) for f in getVTKList(PATH_TO_VTK)]

chemOptions = ["DOC","BAP","EPS","XAR","XN","XDN","XI","O2","NH4","NO3","POCr","n","perm"]


col1,col2 = st.columns([1,5])
chem = col1.selectbox("Select the variable:",options=chemOptions,index=1)

if chem:

    cmap = get_cmap('copper')
    deltaTime = 0.5
    totalFiles = len(all_vtk_paths)

    fig = go.Figure()
    fig.update_layout(
        autosize=False,
        height=600,
        title={
            'text': chem,
            'y': 0.9,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        xaxis_title="Concentration",
        yaxis_title="Depth",
        legend_title="Time",
        font={
            'size': 14,
            'color': "RebeccaPurple"}
        )
    
    if chem == 'perm': 
        fig.update_xaxes(type="log")

    for i,vtk_file in enumerate(all_vtk_paths):
        timestep = i * deltaTime
        if timestep % 1 == 0:
            grid = pv.read(vtk_file)
            color = to_hex(cmap(i/totalFiles))
            dataLine = grid.sample_over_line((0, 0, 0.52),(0, 0, 0),resolution=100)
            depthArray = np.where(dataLine['vtkValidPointMask'],dataLine['Distance'],np.NaN)
        
            fig.add_trace(go.Scatter(
                    x = dataLine[chem],y = depthArray,
                    name=str(timestep),
                    mode='lines',
                    line=dict(color=color,width=5)))
      
    col2.plotly_chart(fig,use_container_width=True)