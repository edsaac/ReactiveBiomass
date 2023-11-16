import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from pathlib import Path
from typing import Literal
from myusefultools.pyopenfoam import OpenFOAM
import multiprocessing as mp

cmaps = st.session_state.cmaps
linecolors = st.session_state.linecolors

def show_nice_option(i):
    return [r"ρ_X", "[DOC]", r"κ"][i]


def main():
    RHO_X_LIST = st.session_state.rho_x_list
    DOC_LIST = st.session_state.doc_list
    KAPPA_LIST = st.session_state.d_growth_list
    LISTS = [RHO_X_LIST, DOC_LIST, KAPPA_LIST]

    dynamic_var = st.sidebar.selectbox("Select a variable", range(3), format_func=show_nice_option)

    looped_list = LISTS[dynamic_var]

    with st.sidebar.form("Select case:"):
        if dynamic_var != 0:
            rho_x = st.selectbox(r"$\rho_X$", RHO_X_LIST, index=2)
        if dynamic_var != 1:
            doc_val = st.selectbox(r"$[\mathsf{DOC}]$", DOC_LIST, index=2)
        if dynamic_var != 2:
            d_growth = st.selectbox(r"$\kappa$", KAPPA_LIST, index=2)
        
        submit_btn = st.form_submit_button("Submit", use_container_width=True)
    
    if submit_btn:
        ## Check CASES folder exists
        ofs = list()

        for val in looped_list:
            # identifier = Path(f"../../experiments/Rosenzweig_2011/fit_ravid/CASES_{doc_val}mgDOC") / f"rhox_{rho_x}__dgrowth_{d_growth}"
            if dynamic_var == 0:
                identifier = Path(f"../../experiments/Rosenzweig_2011/fit_ravid/CASES_{doc_val}mgDOC") / f"rhox_{val}__dgrowth_{d_growth}"
            elif dynamic_var == 1:
                identifier = Path(f"../../experiments/Rosenzweig_2011/fit_ravid/CASES_{val}mgDOC") / f"rhox_{rho_x}__dgrowth_{d_growth}"
            elif dynamic_var == 2:
                identifier = Path(f"../../experiments/Rosenzweig_2011/fit_ravid/CASES_{doc_val}mgDOC") / f"rhox_{rho_x}__dgrowth_{val}"

            template_folder = Path("template")

            if not identifier.exists(): 
                st.toast(f"Some case was not found: \n`{str(identifier)}`")

            ofs.append(OpenFOAM(identifier, template_folder, False))
        
        print(ofs)

        with mp.Pool() as pool:
            
            data_list = pool.starmap(
                OpenFOAM.read_as_xarray_dataset,
                [(of, ["XAR", "EPS", "XI"]) for of in ofs],
            )

        totalsxar = [data.XAR.sum(dim="z") for data in data_list]
        totalseps = [data.EPS.sum(dim="z") for data in data_list]
        totalsxi = [data.XI.sum(dim="z") for data in data_list]

        totalsum = [xar+eps+xi for xar,eps,xi in zip(totalsxar, totalseps, totalsxi)]

        tabs = st.tabs([str(i) for i in looped_list])

        for k,tab,data in zip(looped_list, tabs, data_list):
        
            fig, ax = plt.subplots()

            totalxar = data.XAR.sum(dim="z")
            totaleps = data.EPS.sum(dim="z")
            totalxi = data.XI.sum(dim="z")

            ax.stackplot(data.XAR.t/3600, totalxar, totaleps, totalxi, colors=linecolors)
            ax.ticklabel_format(useMathText=True)
            ax.spines.right.set_visible(False)
            ax.set_xlabel("Time $t$ [d]")
            ax.set_ylabel("Total biomass $\int{X}dV$ [g/L * V]")

            ax.set_title(k)
            ax.set_ylim(0, max([float(s.max()) for s in totalsum]))

            with tab:
                st.pyplot(fig)

if __name__ == "__main__":
    main()
