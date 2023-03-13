import pandas as pd
import numpy as np
import streamlit as st
import xarray as xr
import plotly.graph_objects as go

st.set_page_config(initial_sidebar_state="collapsed")
st.markdown(st.session_state.css_style, unsafe_allow_html=True)
st.title("🌊")
# st.warning("Add Sw to XTotal to get the experimental volumetric water content?")

casepath = st.session_state.casepath
ModelResult = st.session_state.ModelResult
get_results_blob = st.session_state.get_results_blob

####### Simulation data #########################
"## **💻 Numerical simulation**"
hydraulics = get_results_blob("hydraulics")

selection = st.selectbox(
    "Select variable:", 
    hydraulics.keys(), 
    format_func=lambda x:hydraulics[x].long_name
    )

hydraulics[selection].plot_heatmap()
st.plotly_chart(hydraulics[selection].figure, use_container_width=True)

####### Experimental data #######################
"## **🧪 Experimental data**"
is_control = st.checkbox("Show control experiment?")
sheet_name = ("Inoculated", "Control")[is_control]

tabs = st.tabs(["**Water content** (θ)", "**Head** (h)"])
z = np.array([0.57, 0.52, 0.47, 0.42, 0.37, 0.32])

with tabs[0]:
    ### Format experimental head data as a ModelResult object
    θ = pd.read_excel(
        "./experimentalData/.hiddendata/TDR-theta.xlsx", 
        sheet_name=sheet_name)

    t = θ["Time (hr)"].to_numpy()*(60*60)  # convert to seconds
    θ = θ[[f"θ{i}" for i in range(1,7)]].to_numpy()

    data = xr.DataArray(
        θ.T, 
        dims=("z","t"), 
        coords={
            "z": z, 
            "t": t})

    exp_θ = ModelResult(
        var_name="theta", data=data,
        equation_repr = "θ",
        long_name="Vol. water content", 
        is_vtk=False, is_logscale=False,
        zmin=0.15, zmax=0.22, colormap="blues")
    
    exp_θ.plot_heatmap()
    st.plotly_chart(exp_θ.figure, use_container_width=True)

with tabs[1]:
    ### Format experimental head data as a ModelResult object
    head = pd.read_excel(
        "./experimentalData/.hiddendata/matricHead.xlsx", 
        sheet_name=sheet_name)

    t = head["Time (min)"].to_numpy()*60 # convert to seconds
    h = head[[i for i in head.columns if "h_" in i]].to_numpy()/100.


    data = xr.DataArray(
        h.T, 
        dims=("z","t"), 
        coords={
            "z": z, 
            "t": t})

    exp_head = ModelResult(
        var_name="head", data=data,
        long_name="Hydraulic head", 
        is_vtk=False, is_logscale=False,
        zmin=-0.45, zmax=-0.10, colormap="blues")

    exp_head.plot_heatmap()
    st.plotly_chart(exp_head.figure, use_container_width=True)

"## 💻 Post-processing"

r"""
For the simulation, the water saturation is calculated as

$$
\begin{equation}
    S_w(z,t) = \dfrac{V_{\textrm{water}}(z,t)}{V_{\textrm{voids}}(z,t)} = \dfrac{V_{\textrm{water}}(z,t)}{V_{\textrm{voids,t=0}} - V_{\textrm{biomass}}(z,t)}
\end{equation}
$$

In the experiments, the water content is 

$$
\begin{equation}
    S'_{w,b}(z,t) = \dfrac{V_{\textrm{water}}(z,t) + f_{\textrm{w/b}}V_{\textrm{biomass}}(z,t)}{V_{\textrm{voids,t=0}}}
\end{equation}
$$

This can be calculated form the simulation as

$$
S_{w,b}(z,t) = \dfrac{n(z,t)\,S_w(z,t) + f_{\textrm{w/b}}\,\dfrac{\sum{X_i}(z,t)}{\rho_X}} {n(z,t) + f_{\textrm{w/b}}\,\dfrac{\sum{X_i}(z,t)}{\rho_X}}
$$

"""


post_hydr = get_results_blob("postprocessed_hydraulics")
biomass = get_results_blob("biomass")
post_biom = get_results_blob("postprocessed_biomass")

post_biom["Xactive"].data = biomass["XAR"].data + biomass["XN"].data + biomass["XDN"].data
post_biom["Xinactive"].data = biomass["XI"].data + biomass["EPS"].data
post_biom["Xtotal"].data = post_biom["Xactive"].data + post_biom["Xinactive"].data

RHOX = 0.5
post_hydr["Swb"].data = ( (hydraulics["porosity"].data * hydraulics["Sw"].data) + (post_biom["Xtotal"].data / RHOX) ) / ( hydraulics["porosity"].data + (post_biom["Xtotal"].data / RHOX) )

post_hydr["Swb"].plot_heatmap()
st.plotly_chart(post_hydr["Swb"].figure, use_container_width=True)