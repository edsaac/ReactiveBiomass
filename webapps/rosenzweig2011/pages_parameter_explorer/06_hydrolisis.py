import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from pathlib import Path
from typing import Literal
from myusefultools.pyopenfoam import OpenFOAM
import multiprocessing as mp
import colorcet as cc
import pandas as pd
import numpy as np

st.set_page_config(page_title=None, page_icon=None, layout="wide", initial_sidebar_state="auto", menu_items=None)

cmaps = st.session_state.cmaps

def main():

    RHO_X_LIST = st.session_state.rho_x_list
    KAPPA_LIST = st.session_state.d_growth_list
    HYDROLYSIS_DICT = st.session_state.hydrolysis_dict

    with st.sidebar.form("Select case:"):
        doc_val = 50
        rho_x = st.selectbox(r"$\rho_X$", RHO_X_LIST, index=2)
        d_growth = st.selectbox(r"$\kappa$", KAPPA_LIST, index=2)
        z_probe = st.slider("z", 0.00, 0.60, 0.58, 0.01)
        submit_btn = st.form_submit_button("Submit", use_container_width=True)
    
    if submit_btn:

        groups = dict(
                biomass = ["XAR", "EPS", "XI"],
                dissolved = ["DOC", "O2", "NH4"],
                tracers = ["NO3", "tracer"],
                hydraulics = ["h", "Sw", "porosity"]
            )
        
        for hydrolysis, folder_path in HYDROLYSIS_DICT.items():
            
            st.subheader(f"Hydrolysis Rate `{hydrolysis}`")

            ## Check CASES folder exists
            cases_folder = Path(f"{folder_path}/CASES_{doc_val}mgDOC")
            template_folder = Path("template")
        
            if cases_folder.exists():
                
                identifier = cases_folder / f"rhox_{rho_x}__dgrowth_{d_growth}"

                of_high = OpenFOAM(
                    path_case=identifier, 
                    write_to_log = False, 
                    path_template=template_folder
                )
            

                with mp.Pool() as pool:
                    data = pool.starmap(
                        OpenFOAM.read_as_xarray_dataset,
                        [(of_high, group) for group in groups.values()],
                    )
                
                with st.expander("Raw data as an `xarray`"):
                    for i,v in enumerate(data):
                        f"### {i}"
                        st.components.v1.html(v._repr_html_(), height=250, scrolling=True)
            
            else:  
                st.toast("Path does not exist!")

if __name__ == "__main__":
    main()