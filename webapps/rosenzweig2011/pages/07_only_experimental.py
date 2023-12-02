import streamlit as st

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import LogLocator, MaxNLocator

from pathlib import Path
import multiprocessing as mp
from math import log10

from myusefultools.pyopenfoam import OpenFOAM
from pandas import read_excel

cmaps = st.session_state.cmaps
linecolors = st.session_state.linecolors
POROSITY_0 = 0.385

# 🏁 End of experiment results
EXPERIMENT_PATH = st.session_state.experiments_path

cfu = read_excel(EXPERIMENT_PATH / "Live-counts.xlsx", sheet_name="Inoculated")
protein = read_excel(EXPERIMENT_PATH / "Protein-bradfordMethod.xlsx", sheet_name="Inoculated")
theta = read_excel(EXPERIMENT_PATH / "theta-endexperiment.xlsx", sheet_name="Inoculated")
tdr = read_excel(EXPERIMENT_PATH / "TDR-theta.xlsx", sheet_name="Inoculated")

def main():
    
    with st.expander("Dataframes", expanded=False):
        which_df = st.multiselect("Which to show", ["CFU", "Protein", "θ(z)", "θ(t)"])
        if "CFU" in which_df:
            st.dataframe(cfu, use_container_width=True)
        if "Protein" in which_df:
            st.dataframe(protein, use_container_width=True)
        if "θ(z)" in which_df:
            st.dataframe(theta, use_container_width=True)
        if "θ(t)" in which_df:
            st.dataframe(tdr, use_container_width=True)

    tabs = st.tabs(["End of experiment", "Time series"])

    ## End of experiment
    with tabs[0]:
        fig,axs = plt.subplots(1,3, figsize=[10,5], sharey=True, gridspec_kw=dict(wspace=0.1))
        
        ## CFU
        ax=axs[0]
        ax.errorbar(cfu["CFU/mL"], cfu["z (m)"], xerr=cfu["stdev"],
            fmt="none",ecolor=linecolors[0], capsize=2)
        ax.plot(cfu["CFU/mL"], cfu["z (m)"],
            marker="o",markersize=5, c=linecolors[0], lw=1, ls="dashed",)
        ax.set_xscale("log")
        ax.set_ylim(bottom=0.2)
        ax.set_xlim(1e7, 1e10)
        ax.set_xlabel("Microbial counts\n[CFU/mL]")    
        ax.set_ylabel("Depth [m]", fontsize=12)
        ax.xaxis.set_major_locator(LogLocator())
        ax.xaxis.set_minor_locator(LogLocator(subs="all"))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_minor_locator(MaxNLocator(17))

        ## Proteins
        ax=axs[1]
        ax.plot(protein["mg protein/g dry sand"], protein["z (m)"],
            marker="o", markersize=5, c=linecolors[1], lw=1, ls="dashed",)
        ax.set_xscale("log")
        ax.set_xlim(2e-4, 9e-1)
        ax.set_xlabel("Protein content\n[mg protein/g dry sand]")
    
        ## Water content at the end
        ax=axs[2]
        ax.plot(theta["θ"]/POROSITY_0, theta["z (m)"], 
            marker="o", markersize=5, c=linecolors[2], lw=2, ls="dashed",)
        ax.set_xlabel("Water saturation\n$S_w$ [Vw/Vv]",)
        ax.set_xlim(0.1, 1.0)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.xaxis.set_minor_locator(MaxNLocator(10))

        st.pyplot(fig)
    
    with tabs[1]:
        
        labels = ["0.57", "0.52", "0.47"]  # z-coords
        marker_style = dict(linestyle=':', color='0.8', markersize=10,
                            markerfacecolor="tab:blue", markeredgecolor="tab:blue")

        fig, ax = plt.subplots(figsize=[6,3])
        for i, label in enumerate(labels, start=1):
            ax.plot(
                tdr["Time (hr)"], 
                tdr[f"θ{i}"]/POROSITY_0,
                marker=".", mec=linecolors[i], ms=3, mfc=linecolors[i], 
                lw=0.5, ls=":",
            )
            ax.plot(
                tdr["Time (hr)"], 
                tdr[f"θ{i}"].rolling(6).mean()/POROSITY_0, 
                label=f"{label} cm", 
                c=linecolors[i], lw=3
            )
        
        
        
        ax.set_ylim(0.2, 0.6)
        ax.set_xlim(left=0)
        ax.set_xlabel("Time [hr]",)
        ax.set_ylabel("Water saturation\n$S_w$" + R"[$\mathrm{V_w}$/$\mathrm{V_{v0}}$]",)

        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_minor_locator(MaxNLocator(17))
        ax.legend(
            title=R"$z$",
            title_fontproperties=dict(size=10, weight=100),
            loc="lower right",
        )
        st.pyplot(fig)
    

if __name__ == "__main__":
    main()