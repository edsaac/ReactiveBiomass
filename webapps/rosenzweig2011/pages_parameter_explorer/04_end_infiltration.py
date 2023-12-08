import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patheffects as pe
from pathlib import Path
from typing import Literal
from myusefultools.pyopenfoam import OpenFOAM
import multiprocessing as mp
from itertools import product, cycle
import numpy as np


cmaps = st.session_state.cmaps
linecolors = st.session_state.linecolors

def show_nice_option(i):
    return [r"ρ_X", "[DOC]", r"κ"][i]

def main():
    RHO_X_LIST = st.session_state.rho_x_list
    DOC_LIST = st.session_state.doc_list
    KAPPA_LIST = st.session_state.d_growth_list
    LISTS = [RHO_X_LIST, DOC_LIST, KAPPA_LIST]
    HYDROLYSIS_DICT = st.session_state.hydrolysis_dict

    static_var = st.sidebar.selectbox("Select a variable", range(3), format_func=show_nice_option)
    looped_list = LISTS[static_var]

    with st.sidebar.form("Select case:"):
        var_value = st.selectbox(show_nice_option(static_var), looped_list, index=2)
        submit_btn = st.form_submit_button("Submit", use_container_width=True)
    
    if submit_btn:
        for hydrolysis, folder_path in HYDROLYSIS_DICT.items():
            
            st.title(f"Hydrolysis Rate `{hydrolysis}`")
            ## Check CASES folder exists
            ofs = list()

            template_folder = Path(f"{folder_path}/template")

            if static_var == 0:  #<- rho_x
                rho_val = var_value
                xlabel, ylabel = show_nice_option(1), show_nice_option(2)
                for doc_val, d_growth in product(x_list:= DOC_LIST, y_list:= KAPPA_LIST):
                    identifier = Path(f"{folder_path}/CASES_{doc_val}mgDOC") / f"rhox_{rho_val}__dgrowth_{d_growth}"
                    
                    if identifier.exists(): 
                        ofs.append(OpenFOAM(identifier, template_folder, False))
                    else:
                        st.toast(f"Some case was not found: \n`{str(identifier)}`")
            
            elif static_var == 1: #<- DOC
                doc_val = var_value
                xlabel, ylabel = show_nice_option(0), show_nice_option(2)
                for rho_val, d_growth in product(x_list:=RHO_X_LIST, y_list:= KAPPA_LIST):
                    identifier = Path(f"{folder_path}/CASES_{doc_val}mgDOC") / f"rhox_{rho_val}__dgrowth_{d_growth}"
                    
                    if identifier.exists(): 
                        ofs.append(OpenFOAM(identifier, template_folder, False))
                    else:
                        st.toast(f"Some case was not found: \n`{str(identifier)}`")

            elif static_var == 2:  #<- kappa
                d_growth = var_value
                xlabel, ylabel = show_nice_option(0), show_nice_option(1)
                for rho_val, doc_val in product(x_list:=RHO_X_LIST, y_list:= DOC_LIST):
                    identifier = Path(f"{folder_path}/CASES_{doc_val}mgDOC") / f"rhox_{rho_val}__dgrowth_{d_growth}"
                    
                    if identifier.exists(): 
                        ofs.append(OpenFOAM(identifier, template_folder, False))
                    else:
                        st.toast(f"Some case was not found: \n`{str(identifier)}`")
            
            for of in ofs: 
                of.process_boundaryProbes()
            
            
            datas = [of.boundaryProbes[1].array_data["Uz"] for of in ofs]
            flow_end = np.array([float(data.isel(probe=1).isel(time=-1))/float(data.isel(probe=1).isel(time=0)) for data in datas]).reshape(len(x_list), len(y_list))

            # flow_end

            fig, ax = plt.subplots(dpi=120)
            im = ax.imshow(flow_end, cmap="Blues_r", vmin=0, vmax=1, aspect="auto")
            cb = fig.colorbar(im, shrink=0.5, aspect=10, pad=0.05)

            cb.ax.set_ylabel(r"$\dfrac{\mathsf{End\,flow\,rate}}{\mathsf{Initial\,flow\,rate}}$", fontsize=9, labelpad=10)
            cb.ax.tick_params(labelsize=9)
            cb.ax.set_yticks(np.linspace(0,1,3))

            # Show all ticks and label them with the respective list entries
            ax.set_xticks(np.arange(len(y_list)), labels = [str(y) for y in y_list])
            ax.set_yticks(np.arange(len(x_list)), labels = [str(x) for x in x_list])

            # labels=["$10^{{ {:.{}f} }}$".format(log10(float(i)), j) for i,j in zip(y_list, cycle([0,1]))])
            # labels=["0"] + [f"$10^{{{log10(i):.0f}}}$" for i in x_list[1:]])

            ax.tick_params(labelsize=9)
            ax.set_title(f"With {show_nice_option(static_var)} = {var_value}", fontsize=8)

            # Loop over data dimensions and create text annotations.
            for i in range(len(x_list)):
                for j in range(len(y_list)):
                    text = ax.text(
                        j, i, 
                        f"{flow_end[i, j]*100:.1f}%",
                        path_effects=[pe.withStroke(linewidth=2, foreground="white")],
                        ha="center", va="center", color="k", fontdict=dict(size=9)
                    )

            ax.set_ylabel(xlabel)
            ax.set_xlabel(ylabel)
            # ax.ticklabel_format(useMathText=True)

            ax.set_xticks(np.arange(len(y_list))-.5, minor=True)
            ax.set_yticks(np.arange(len(x_list))-.5, minor=True)
            ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
            ax.tick_params(which="minor", bottom=False, left=False)
            cols = st.columns([1,4,1])
            cols[1].pyplot(fig)

if __name__ == "__main__":
    main()