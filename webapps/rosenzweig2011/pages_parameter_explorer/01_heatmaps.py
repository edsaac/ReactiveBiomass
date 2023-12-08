import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from pathlib import Path
from typing import Literal
from myusefultools.pyopenfoam import OpenFOAM
import multiprocessing as mp
import colorcet as cc

st.set_page_config(page_title=None, page_icon=None, layout="wide", initial_sidebar_state="auto", menu_items=None)

cmaps = st.session_state.cmaps

def draw_heatmaps(data):

    is_biomass = "XAR" in list(data.data_vars)
    is_dissolved = "O2" in list(data.data_vars)
    
    ## Heatmaps
    igt = 0
    fig, (caxs, axs) = plt.subplots(
        2, len(data), figsize=[15, 5], dpi=120,
        gridspec_kw={"height_ratios": [0.2, 5]}, sharey="row"
    )

    for ax,cax,cmap, (name,scalar) in zip(axs, caxs, cmaps, data.items()):
        
        img = ax.pcolormesh(scalar.t[igt:] / 86400, scalar.z, scalar[:, igt:])
        
        if is_biomass: 
            img.set_norm(colors.LogNorm(5e-5, 1e-0))
            img.set_cmap(cmap)
        else:
            img.set_cmap("cet_gray_r")
        
        if is_dissolved:
            img.set_clim(vmin=0)

        ax.spines.right.set_visible(False)
        ax.set_xlabel("Time $t$ [d]")
        plt.colorbar(img, cax=cax, orientation="horizontal")
        cax.set_title(rf"{name} [g/L]")

    axs[0].set_ylabel("Depth $z$ [m]")
    fig.tight_layout()
    return fig

def main():

    RHO_X_LIST = st.session_state.rho_x_list
    DOC_LIST = st.session_state.doc_list
    KAPPA_LIST = st.session_state.d_growth_list
    HYDROLYSIS_DICT = st.session_state.hydrolysis_dict

    with st.sidebar.form("Select case:"):
        doc_val = st.selectbox(r"$[DOC]$", DOC_LIST, index=2)
        rho_x = st.selectbox(r"$\rho_X$", RHO_X_LIST, index=2)
        d_growth = st.selectbox(r"$\kappa$", KAPPA_LIST, index=2)

        submit_btn = st.form_submit_button("Submit", use_container_width=True)
    
    if submit_btn:
        
        for hydrolysis, folder_path in HYDROLYSIS_DICT.items():
            
            f"### Hydrolysis Rate `{hydrolysis}`"
            
            ## Check CASES folder exists
            cases_folder = Path(f"{folder_path}/CASES_{doc_val}mgDOC")
            template_folder = Path(f"{folder_path}/template")
            
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
                    dissolved = ["DOC", "O2", "NH4"],
                    tracers = ["NO3", "tracer"],
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
                
                with mp.Pool() as pool:
                    figures = pool.map(draw_heatmaps, data)
                
                for fig in figures:
                    st.pyplot(fig)
            
if __name__ == "__main__":
    main()