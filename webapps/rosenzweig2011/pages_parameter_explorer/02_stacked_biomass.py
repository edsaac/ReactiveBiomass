import streamlit as st

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import LogLocator

from pathlib import Path
import multiprocessing as mp
from math import log10

from myusefultools.pyopenfoam import OpenFOAM
from pandas import read_excel

cmaps = st.session_state.cmaps
linecolors = st.session_state.linecolors

# 🏁 End of experiment results
EXPERIMENT_PATH = st.session_state.experiments_path
cfu = read_excel(EXPERIMENT_PATH / "Live-counts.xlsx", sheet_name="Inoculated")
protein = read_excel(
    EXPERIMENT_PATH / "Protein-bradfordMethod.xlsx", sheet_name="Inoculated"
)


def show_nice_option(i):
    return [r"ρ_X", "[DOC]", r"κ"][i]


def plot_stacked(k, data):
    ## Some data processing
    totalxar = data.XAR.sum(dim="z")
    totaleps = data.EPS.sum(dim="z")
    totalxi = data.XI.sum(dim="z")

    fig_stacked, ax = plt.subplots(figsize=[5, 3])

    ## Stacked plot over time
    ax.stackplot(
        data.XAR.t / 86400,
        totalxar,
        totaleps,
        totalxi,
        colors=linecolors,
        labels=[R"$X_{\mathsf{aerobic\;resp.}}$", R"$X_{\mathsf{labile}}$", R"$X_{\mathsf{recalcitrant}}$"],
        # baseline='zero',
    )

    ax.ticklabel_format(useMathText=True)
    ax.spines.right.set_visible(False)
    ax.set_xlabel("Time $t$ [d]")
    ax.set_ylabel("Total biomass\n" + R"$\int{X_j} \; dV$")
    # ax.set_title(f"With ${show_nice_option(dynamic_var)}$ = {k}")
    # ax.set_ylim(0, max([float(s.max()) for s in totalsum]))
    # ax.set_ylim(0, 10.0)
    ax.legend(loc="upper left")

    ##############################

    fig_profiles, axs = plt.subplots(1,2, figsize=[8,4], sharey=True, gridspec_kw=dict(wspace=0.1))

    ## Final active biomass over depth
    ax = axs[0]
    l, = ax.plot(data.XAR.isel(t=-1), data.XAR.z,
        label="Simulation",
        c=linecolors[0], alpha=0.9, lw=3, ls="dashed")
    ax.set_xscale("log")
    ax.set_xlim(2e-6, 2e-0)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("Depth [m]")
    ax.set_xlabel(R"Active biomass [g/L]")
    ax.xaxis.set_major_locator(LogLocator())
    ax.xaxis.set_minor_locator(LogLocator(subs="auto"))

    ax2 = ax.twiny()
    ax2.errorbar(cfu["CFU/mL"], cfu["z (m)"], xerr=cfu["stdev"],
        fmt="none",ecolor="gray", capsize=2
        )
    l2, = ax2.plot(cfu["CFU/mL"], cfu["z (m)"],
        label=R"Experimental",
        marker="o",markersize=5, mfc = linecolors[0],
        c="gray", lw=1, ls="dashed", 
        )
    ax2.set_xscale("log")
    ax2.set_ylim(ax.get_ylim())
    ax2.spines.top.set_visible(True)
    ax2.set_xlim(5e6, 5e10)
    ax2.set_xlabel("Microbial counts\n[CFU/mL]")
    ax2.legend(
        [l, l2],
        [l.get_label(), l2.get_label()],
        title=R"$X_{\mathsf{aerobic\;resp.}}$",
        title_fontproperties=dict(size=10, weight=100),
        loc="lower right",
    )
    ax2.spines.top.set_position(('axes', 1.05))

    ## Final inactive biomass over depth
    ax = axs[1]
    (l,) = ax.plot(
        data.EPS.isel(t=-1) + data.XI.isel(t=-1),
        data.EPS.z,
        label="Simulation",
        c=linecolors[1],
        alpha=0.9,
        lw=3,
        ls="dashed",
    )

    ax.set_xscale("log")
    ax.set_xlim(2e-6, 10e-0)
    ax.set_ylim(bottom=0)
    # ax.set_ylabel("Depth [m]")
    ax.set_xlabel(R"Inactive biomass  [g/L]")
    ax.xaxis.set_major_locator(LogLocator())
    ax.xaxis.set_minor_locator(LogLocator(subs="all"))
    # ax.tick_params(axis="both", which="minor", width=1, length=10, color="k")

    DENSITY_SAND = 1_530 # g/L_REV
    MG_TO_G = 1000 
    PROTEIN_PERCENT_IN_EPS = 0.678 # See https://doi.org/10.1016/j.biortech.2009.01.053
    # ax2 = ax.twiny()

    (l2,) = ax.plot(
        protein["mg protein/g dry sand"] * DENSITY_SAND / MG_TO_G / PROTEIN_PERCENT_IN_EPS,  ##<- convert to g/L_REV
        protein["z (m)"],
        label=R"Experimental",
        marker="o",
        markersize=5,
        c="gray",
        lw=1,
        ls="dashed",
        mfc=linecolors[1]
    )

    # ax2.set_xscale("log")
    # ax2.set_ylim(ax.get_ylim())
    # ax2.spines.top.set_visible(True)
    # ax2.set_xlim(8e-6, 2e-0) #(1e-4, 1e0)
    # ax2.set_xlabel("Derived from Protein content\n[g/L]")

    ax.legend(
        [l, l2],
        [l.get_label(), l2.get_label()],
        title=R"$X_{\mathsf{labile}} + X_{\mathsf{recalcitrant}}$",
        title_fontproperties=dict(size=10, weight=100),
        loc="upper left",
    )

    return (fig_stacked, fig_profiles)


def main():
    cols = st.columns(2)
    with cols[0]:
        "### Colony formation units"
        st.dataframe(cfu)
    with cols[1]:
        "### Protein content"
        st.dataframe(protein)

    RHO_X_LIST = st.session_state.rho_x_list
    DOC_LIST = st.session_state.doc_list
    KAPPA_LIST = st.session_state.d_growth_list
    LISTS = [RHO_X_LIST, DOC_LIST, KAPPA_LIST]
    HYDROLYSIS_DICT = st.session_state.hydrolysis_dict

    global dynamic_var
    dynamic_var = st.sidebar.selectbox(
        "Select a variable", range(3), format_func=show_nice_option
    )

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
        for hydrolysis, folder_path in HYDROLYSIS_DICT.items():
            f"### Hydrolysis Rate `{hydrolysis}`"
            
            ## Check CASES folder exists
            ofs = list()

            for val in looped_list:
                # identifier = Path(f"../../experiments/Rosenzweig_2011/fit_ravid/CASES_{doc_val}mgDOC") / f"rhox_{rho_x}__dgrowth_{d_growth}"
                if dynamic_var == 0:
                    identifier = (
                        Path(f"{folder_path}/CASES_{doc_val}mgDOC")
                        / f"rhox_{val}__dgrowth_{d_growth}"
                    )
                elif dynamic_var == 1:
                    identifier = (
                        Path(f"{folder_path}/CASES_{val}mgDOC")
                        / f"rhox_{rho_x}__dgrowth_{d_growth}"
                    )
                elif dynamic_var == 2:
                    identifier = (
                        Path(f"{folder_path}/CASES_{doc_val}mgDOC")
                        / f"rhox_{rho_x}__dgrowth_{val}"
                    )

                template_folder = Path(f"{folder_path}/template")

                if not identifier.exists():
                    st.toast(f"Some case was not found: \n`{str(identifier)}`")

                ofs.append(OpenFOAM(identifier, template_folder, False))

            print(ofs)

            with mp.Pool() as pool:
                data_list = pool.starmap(
                    OpenFOAM.read_as_xarray_dataset,
                    [(of, ["XAR", "EPS", "XI"]) for of in ofs],
                )

            with mp.Pool() as pool:
                figs = pool.starmap(plot_stacked, list(zip(looped_list, data_list)))

            tabs = st.tabs(
                [f"{show_nice_option(dynamic_var)} = {i}" for i in looped_list]
            )
            for tab, fig in zip(tabs, figs):
                with tab:
                    st.pyplot(fig[0])
                    st.pyplot(fig[1])


if __name__ == "__main__":
    main()
