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
    HYDROLYSIS_DICT = st.session_state.hydrolysis_dict
    
    rho_x = 10
    doc_val = 100
    d_growth = 1e-11

    for hydrolysis, folder_path in HYDROLYSIS_DICT.items():
        
        f"### Hydrolysis Rate `{hydrolysis}`"

        ## Check CASES folder exists
        identifier = Path(f"{folder_path}/CASES_{doc_val}mgDOC") / f"rhox_{rho_x}__dgrowth_{d_growth}"
        template_folder = Path(f"../../exploration/time_flipper/ROSENZ_TEMPLATE")
        
        if not identifier.exists():
            st.toast(f"Some case was not found: \n`{str(identifier)}`")
        
        else:

            of = OpenFOAM(
                    path_case=identifier, 
                    write_to_log = False, 
                    path_template=template_folder
                )

            of.process_boundaryProbes()

            # of.boundaryProbes[0].array_data

            with st.expander("Raw data as an `xarray`"):
                for probe in of.boundaryProbes:
                    st.components.v1.html(probe.array_data._repr_html_(), height=250, scrolling=True)

            # st.components.v1.html(of.boundaryProbes[1].array_data["Ux"]._repr_html_(), height=250, scrolling=True)

            h_at_boundary = dict()
            u_at_boundary = dict()
            k_at_boundary = dict()
            doc_at_boundary = dict()


            for timestep in of.list_times:
                # h_at_boundary[timestep] = float(of.get_value_from_foamDictionary(f"{timestep}/h", "boundaryField.top.value" ))
                u_at_boundary[timestep] = float(of.get_value_from_foamDictionary(f"{timestep}/U", "boundaryField.top.value", is_vector=True)[-1])
                k_at_boundary[timestep] = float(of.get_value_from_foamDictionary(f"{timestep}/hydraulicCond", "boundaryField.top.value" ))
                doc_at_boundary[timestep] = float(of.get_value_from_foamDictionary(f"{timestep}/DOC", "boundaryField.top.value" ))

            # st.write(doc_at_boundary)
            # fig,ax = plt.subplots()
            # ax.plot([float(k)/86400 for k in sw_at_boundary.keys()], k_at_boundary.values())
            # st.pyplot(fig)

            # # Infiltration
            # cols = st.columns(2)
            # fig,ax = plt.subplots()
            # uv = of.boundaryProbes[1].array_data["Uz"]
            # # uv.isel(probe=0).plot.line(ax=ax, label="Out")
            # # ax.plot(uv.time/86400, -uv.isel(probe=1), label="Influx", lw=3)
            # ax.plot([float(k)/86400 for k in u_at_boundary.keys()], [-i for i in u_at_boundary.values()], lw=4)
            # # ax.set_ylim(bottom = 0)
            # ax.set_xlim(left=0)
            # ax.set_yscale("log") 
            # ax.set_title("Infiltration")
            # # ax.ticklabel_format(useMathText=True)
            # ax.set_xlabel("Time $t$ [d]")
            # ax.legend()
            # cols[0].pyplot(fig)

            ## Sw
            cols = st.columns([1,4,1])
            fig,ax = plt.subplots()
            uv = of.boundaryProbes[0].array_data["DOC"]
            # uv.isel(probe=0).plot.line(ax=ax, label="Out")
            # ax.plot(uv.time/86400, uv.isel(probe=0), label="Outflux", lw=3)
            ax.plot(uv.time/86400, uv.isel(probe=1), label="Influx", lw=3)
            ax.plot([float(k)/86400 for k in doc_at_boundary.keys()], doc_at_boundary.values(), lw=4)

            # ax.set_ylim(0, 0.32)
            ax.set_xlim(left=0)
            ax.set_title("h")
            ax.ticklabel_format(useMathText=True)
            ax.set_xlabel("Time $t$ [d]")
            ax.legend()
            cols[1].pyplot(fig)



if __name__ == "__main__":
    main()