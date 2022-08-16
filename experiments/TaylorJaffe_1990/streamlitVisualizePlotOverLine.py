import pyvista as pv
from matplotlib.cm import get_cmap
from matplotlib.colors import to_hex

import plotly.graph_objects as go
import numpy as np
from myusefultools.parser import getVTKList
import os
import pandas as pd
from natsort import natsorted
import re

import streamlit as st
st.set_page_config(layout="wide",initial_sidebar_state="collapsed")
st.markdown(
    '''
    <style>
    .stDataFrame {
        width: 100%;
    }
    </style>
    ''',unsafe_allow_html=True)


PATH_TO_VTK = "./column1/VTK"
all_vtk_paths = [os.path.join(PATH_TO_VTK,f) for f in getVTKList(PATH_TO_VTK)]

PATH_TO_CSV = "./experimentalData"
all_csv_paths = [os.path.join(PATH_TO_CSV,f) for f in os.listdir(PATH_TO_CSV) if "Fig1" in f]
all_csv_paths = natsorted(all_csv_paths)

REFERENCE_PERMEABILITY = 2.57e-10
COLUMN_LENGHT = 0.520
experimentalResults = dict()
pattern = re.compile(r't(\d+)d')

for f in all_csv_paths:
    
    experiment_df = pd.read_csv(f)
    # st.dataframe(experiment_df.style.format('{:.2E}'),height=250)
    experiment_df["Column Distance (m)"] = COLUMN_LENGHT - experiment_df["Column Distance (cm)"]/100
    experiment_df["Permeability (m²)"] = REFERENCE_PERMEABILITY * experiment_df["Permeability Reduction (-)"]
       
    timefromfile = pattern.search(f).group(1)
    experimentalResults[timefromfile] = experiment_df


chemOptions = ["perm","DOC","BAP","EPS","XAR","XN","XDN","XI","O2","NH4","NO3","POCr","n"]

chem = st.sidebar.selectbox("Select the variable:",options=chemOptions,index=0)

if chem:
    cmap = get_cmap('cool')
    deltaTime = 1.0
    totalFiles = len(all_vtk_paths)

    fig = go.Figure()
    
    fig.update_layout(
        autosize=False,
        height=600,
        # width=800,
        title={
            'text': chem,
            'y': 0.9,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        xaxis_title="Value",
        yaxis_title="Depth",
        legend={
            'yanchor':"top",
            'y':-0.05,
            'xanchor':"left",
            'x':0.01,
            'title':None,
            'groupclick':"toggleitem"},
        font={
            'size': 14,
            'color': "darkblue"}
        )
    
    if chem == 'perm': 
        fig.update_xaxes(type="log")

        ## Experimental data
        for k,df in experimentalResults.items():
            fig.add_trace(
                go.Scatter(
                    x = df['Permeability Reduction (-)'],y = df['Column Distance (m)'],
                    name = k + ' d',
                    mode = 'lines+markers',
                    legendgroup = "experiments",
                    legendgrouptitle_text = "Experimental data",
                    marker={'size' : 8,
                            'color': to_hex(cmap(float(k)/totalFiles)),
                            'symbol': "x" 
                            }
                    ))
        

    for i,vtk_file in enumerate(all_vtk_paths):
        timestep = i * deltaTime
        
        if timestep % 7 == 0 and timestep > 0:
            grid = pv.read(vtk_file)
            color = to_hex(cmap(i/totalFiles))
            dataLine = grid.sample_over_line((0, 0, 0.52),(0, 0, 0),resolution=99)
            data = dataLine[chem]

            if chem == "perm": data = data/REFERENCE_PERMEABILITY

            depthArray = np.where(dataLine['vtkValidPointMask'],dataLine['Distance'],np.NaN)
        
            fig.add_trace(go.Scatter(
                    x = data, y = depthArray,
                    name=str(timestep) + " d",
                    mode='lines',
                    legendgroup="model",
                    legendgrouptitle_text="Numerical model",
                    line={
                        'color':color,
                        'width':5}))
      
    st.plotly_chart(fig,use_container_width=True)