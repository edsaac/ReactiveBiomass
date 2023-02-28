import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go

st.set_page_config(layout="wide")

st.session_state["axis_setup"] = dict(
    mirror=False,
    ticks="outside",
    color="#bbb",
    showline=True)

st.session_state["colors"] = ['#e41a1c','#377eb8','#4daf4a','#ff7f00','#ffff33']

st.session_state["colorbar_kwargs"] = dict(
    len=0.5,
    lenmode="fraction",
    thickness=0.04,
    thicknessmode="fraction",
    titlefont_size=15,
    outlinewidth=0,
    x=1.10,
    y=0.50,
    xanchor="center"
    )

def add_flood_and_dry_periods(fig:go.Figure):
    #FP/DP
    fig.add_vrect(
        x0=0, x1=1, line_width=0, 
        fillcolor="#0000ff", opacity=0.1, 
        annotation_text="Flood")
    
    fig.add_vrect(
        x0=1, x1=2, line_width=0, 
        fillcolor="#ffefd5", opacity=0.1, 
        annotation_text="Dry")

    fig.add_vrect(
        x0=2, x1=3, line_width=0, 
        fillcolor="#0000ff", opacity=0.1, 
        annotation_text="Flood")

st.session_state["add_flood_and_dry_periods"] = add_flood_and_dry_periods


"# Improving soil aquifer treatment efficiency using air injection into the subsurface"
"> Experiments by Ido Arad"

st.image(".hiddendata/image1.jpg")

cols = st.columns(2)

with cols[0]:
    r"""
    - Start: 20/12/21 @ 11:10
    - Duration: 72 hours 
        - (24 hr. flooding, 24 hr. drying, 24 hr. flooding)
    - No active air injection
    - Measuring ports and their depth: 
        - IN 2 (0 cm)
        - A (25 cm)
        - B (65 cm)
        - C (105 cm)
        - D (145 cm)
    - 160 cm tall column
    """

with cols[1]:
    r"""
    Available data:
    - Dissolved oxygen [mg/L]
    - Temperature
    - Oxidation-reduction potential [mV]
    - Volumetric water content [%]
    """

r"""
## Inflow composition
Basic parameters of the synthetic effluent in the three main experiments $\pm$ mean ± SD

|Species          |         mg/L|
|-----------------|-------------|
|$\mathsf{NH}_4^+$|$2.62\pm0.98$|
|$\mathsf{TKN}$   |$8.75\pm0.56$|
|$\mathsf{NO}_3^-$|$0.85\pm0.66$|
|$\mathsf{DOC}$   |$41.2\pm1.36$|

"""