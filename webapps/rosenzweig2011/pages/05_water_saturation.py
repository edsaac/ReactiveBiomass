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

### Format experimental head data as a ModelResult object
z_list = [0.57, 0.52, 0.47, 0.42, 0.37, 0.32]
θ = pd.read_excel(
    "../../experiments/Rosenzweig_2011/experimentalData/.hiddendata/TDR-theta.xlsx", 
    sheet_name="Inoculated")
t = θ["Time (hr)"].to_numpy()*(60*60)  # convert to seconds
POROSITY_0 = 0.385

def main():

    RHO_X_LIST = st.session_state.rho_x_list
    DOC_LIST = st.session_state.doc_list
    KAPPA_LIST = st.session_state.d_growth_list
    HYDROLYSIS_DICT = st.session_state.hydrolysis_dict
    
    with st.sidebar.form("Select case:"):
        doc_val = st.selectbox(r"$[DOC]$", DOC_LIST, index=2)
        rho_x = st.selectbox(r"$\rho_X$", RHO_X_LIST, index=2)
        d_growth = st.selectbox(r"$\kappa$", KAPPA_LIST, index=2)
        z_probe = st.slider("z", 0.00, 0.60, 0.58, 0.01)
        submit_btn = st.form_submit_button("Submit", use_container_width=True)
    
    if submit_btn:
        
        for hydrolysis, folder_path in HYDROLYSIS_DICT.items():
            
            st.subheader(f"Hydrolysis Rate `{hydrolysis}`")

            ## Check CASES folder exists
            cases_folder = Path(f"{folder_path}/CASES_{doc_val}mgDOC")
            template_folder = Path("template")
            
            if not cases_folder.exists():
                st.toast("Path does not exist!")
            
            else:
                identifier = cases_folder / f"rhox_{rho_x}__dgrowth_{d_growth}"

                of = OpenFOAM(
                        path_case=identifier, 
                        write_to_log = False, 
                        path_template=template_folder
                    )
                
                groups = dict(
                    biomass = ["XAR", "EPS", "XI"],
                    hydraulics = ["h", "Sw", "porosity"]
                )

                with mp.Pool() as pool:
                    data = pool.starmap(
                        OpenFOAM.read_as_xarray_dataset,
                        [(of, group) for group in groups.values()],
                    )
                
                with st.expander("Raw data as an `xarray`"):
                    for i,v in enumerate(data):
                        f"### {i}"
                        st.components.v1.html(v._repr_html_(), height=250, scrolling=True)
                
                numeric = data[1].sel(z=z_probe, method="nearest")
                sw = numeric["Sw"]
                n = numeric["porosity"]
                time = numeric.t/86_400

                # tabs = st.tabs([f"θ{i}" for i in range(1,7)])

                # for i, (tab,z) in enumerate(zip(tabs, z_list), start=1):
                #     cols = tab.columns(3)
                #     numeric = data[1].sel(z=z, method="nearest")
                #     sw = numeric["Sw"]
                #     n = numeric["porosity"]
                #     time = numeric.t/86_400

                cols = st.columns(3)
                fig, ax = plt.subplots()
                ax.plot(time, n)
                ax.set_xlabel("Time [d]")
                ax.set_ylabel("Porosity $n$ [Vv/Vt]")
                # with tab:
                cols[0].pyplot(fig)
                
                fig, ax = plt.subplots()
                ax.plot(time, sw)
                ax.set_xlabel("Time [d]")
                ax.set_ylabel("Water saturation $Sw$ [Vw/Vv]")
                ax.set_ylim(0,1)
                ax.set_title(z_probe)
                # with tab:
                cols[1].pyplot(fig)

                eff_sw = ((sw*n) + 0.99 * (POROSITY_0 - n))/POROSITY_0

                fig, ax = plt.subplots()
                ax.plot(time, eff_sw, label="Numerical")
                
                try: 
                    i = z_list.index(z_probe) + 1
                    ax.plot(t/86_400, θ[f"θ{i}"]/POROSITY_0, label="Lab")
                except ValueError:
                    ...
                ax.set_xlabel("Time [d]")
                ax.set_ylabel("Water saturation with biomass $S'w$ [(Vw + Vb)/Vv]")
                ax.set_ylim(0,1)
                ax.set_title(z_probe)
                # with tab:
                cols[2].pyplot(fig)

        ...

if __name__ == "__main__":
    main()