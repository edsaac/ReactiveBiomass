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

st.header("Taylor & Jaffé (1990)",anchor="main_title")
PATH_TO_VTK = f"./column1/VTK"

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

selected_chem = st.selectbox("Select the variable:",options=chemOptions)

if selected_chem:
    cmap = get_cmap('cool')
    deltaTime = 1.0
    totalFiles = len(all_vtk_paths)

    fig = go.Figure()
    
    fig.update_layout(
        autosize=False,
        hovermode='closest',
        height=800,
        title={
            'text': selected_chem,
            'y': 0.9,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        xaxis={
            'title':"Value",
            'exponentformat' : "power",
            'gridwidth':1,
            'gridcolor':"CadetBlue"},
        yaxis={
            'title':"Depth (m)",
            'exponentformat' : "power",
            'gridwidth':0.2,
            'gridcolor':"CadetBlue"},
        legend={
            'orientation':"h",
            'yanchor':"top",
            'y':-0.1,
            'xanchor':"left",
            'x':0.01,
            'title':None,
            'groupclick':"toggleitem"},
        font={
            'size': 14}
        )
    
    if selected_chem == 'perm': 
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
        
        if timestep % 14 == 0 and timestep > 0:
            grid = pv.read(vtk_file)
            color = to_hex(cmap(i/totalFiles))
            dataLine = grid.sample_over_line((0, 0, 0.52),(0, 0, 0),resolution=99)
            data = dataLine[selected_chem]

            if selected_chem == "perm": data = data/REFERENCE_PERMEABILITY
            if selected_chem == "U": data = data.T[-1]  #Shows only the Z-component

            depthArray = np.where(dataLine['vtkValidPointMask'],dataLine['Distance'],np.NaN)
        
            fig.add_trace(go.Scatter(
                    x = data, y = depthArray,
                    name=str(timestep) + " d",
                    mode='lines',
                    legendgroup="model",
                    legendgrouptitle_text="Numerical model",
                    line={
                        'color':color,
                        'width':5},
                    hovertemplate = '<b>%{x:.2E}</b>'))
    
    st.plotly_chart(fig,use_container_width=True)
    st.write(dataLine)

st.markdown("""<style>
h2 {
  position: sticky;
  top: 10px;
  padding: 5px;
  background-color: #cae8ca;
  border: 2px solid #4CAF50;
  text-align:left;
}
</style>""",unsafe_allow_html=True)