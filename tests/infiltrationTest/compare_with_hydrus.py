from collections import namedtuple
from subprocess import check_output
import subprocess

import streamlit as st

import matplotlib.pyplot as plt
import matplotlib as mpl

import numpy as np
import pandas as pd

import pyvista as pv
from stpyvista import stpyvista

from myusefultools import parser
from phydrus.read import read_nod_inf

import plotly.graph_objects as go

if "first_run" not in st.session_state:
    st.session_state.first_run = True

    ## Check if xvfb is already running on the machine
    is_xvfb_running = subprocess.run(["pgrep", "Xvfb"], capture_output=True)

    if is_xvfb_running.returncode == 1:
        st.toast("[red](Xvfb was not running...)")
        pv.start_xvfb()
    else:
        st.toast(f"Xvfb is running! \n\n`PID: {is_xvfb_running.stdout.decode('utf-8')}`")


REPO_PATH = (
    check_output(["git", "rev-parse", "--show-toplevel"]).decode("utf-8").strip()
)
plt.style.use(f"{REPO_PATH}/misc/edwin.mplstyle")

"# Tests for unsaturated flow solver"

OPTIONS = ["LOAMYSAND", "THREELAYERS", "TOPCLOG", "HANGCOL"]
COLUMN_LENGTH_DICT = {k: v for k, v in zip(OPTIONS, [6, 6, 5, 1])}  # m
CASE = st.selectbox("Select test case", options=OPTIONS)
COLUMN_LENGTH = COLUMN_LENGTH_DICT[CASE]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Retrive Hydrus-1D results
# It assumes the units are L:meters and T:seconds
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HYDRUS_NOD_INF = f"{CASE}/HYDRUS_RESULTS/Nod_Inf.out"
hydrus_profiles = read_nod_inf(HYDRUS_NOD_INF)

# Remove time zero
del hydrus_profiles[0.0]

# Convert to hours for readability
hydrus_times_hours = [t / 3600 for t in hydrus_profiles.keys()]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sample colors for plotting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cmap = mpl.cm.Dark2
colors = cmap(np.linspace(0, 1, len(hydrus_times_hours)))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Retrive OpenFOAM results
# It assumes the units are L:meters and T:seconds
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OPENFOAM_FOLDER = f"{CASE}/{CASE}"
of_vtk_files = parser.getVTKList(f"{OPENFOAM_FOLDER}/VTK")
of_times_seconds = parser.getTimeList(OPENFOAM_FOLDER)
of_times_hours = np.array([float(t) / 3600 for t in of_times_seconds])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Visualize soil propertires distribution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
with st.sidebar:
    plotter = pv.Plotter()
    of_vtk_soil_properties = parser.getVTKList(f"{OPENFOAM_FOLDER}/VTK_soilProperties")

    mesh = pv.read(f"{OPENFOAM_FOLDER}/VTK_soilProperties/{of_vtk_soil_properties[0]}")
    soil_property = st.selectbox("Soil property", list(set(mesh.array_names)), index=3)
    plotter.add_mesh(mesh, scalars=soil_property, cmap="bwr")
    plotter.background_color = "#f0f2f6"
    plotter.camera.zoom("tight")
    plotter.view_xz()
    plotter.camera.azimuth = 45
    plotter.camera.elevation = 25
    plotter.window_size = [240, 450]
    stpyvista(plotter, horizontal_align="center")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check correspondance between times from both simulations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
with st.expander("VTK <-> Hydrus time step comparisons"):
    """
    ## From Hydrus-1D:
    The following times were found
    """
    st.code(hydrus_times_hours)

    """
    ## From OpenFOAM:
    The following times were found
    """
    st.code(of_times_hours)

    of_selected_vtk_files = dict()
    for vtk, time_s, time_hr in zip(of_vtk_files, of_times_seconds, of_times_hours):
        # if round(float(time_hr), 2) in hydrus_times_hours:
        if min(np.abs(np.array(hydrus_times_hours) - float(time_hr))) < 0.01:
            vtk, time_s, time_hr, "*"
            of_selected_vtk_files[round(float(time_hr), 2)] = vtk
        else:
            vtk, time_s, time_hr

if len(hydrus_profiles) == len(of_selected_vtk_files) == len(colors):
    st.toast("1:1 times found!")
else:
    st.toast("[red](**Something went wrong**)")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot using matplotlib
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"## Results comparison"
VARIABLE_LABELS = ["h(z)", "θ(h)", "K(h)", "C(h)"]
tabs = st.tabs([f"**{v}**" for v in VARIABLE_LABELS])
of_pandas_summary = dict()

PLOT_ENGINE = st.sidebar.radio("Plot with", ["matplotlib", "plotly"], index=1)

Plot = namedtuple("Plot", ["fig", "ax"])

if PLOT_ENGINE == "matplotlib":
    plot_kwargs = {
        "figure": dict(figsize=[5, 7]),
        "openfoam": dict(lw=5, alpha=0.5, zorder=1),
        "hydrus1d": dict(lw=1.5, ls="dashed", zorder=2),
    }

    plots = {k: Plot(*plt.subplots(**plot_kwargs["figure"])) for k in VARIABLE_LABELS}

    for (time_hydrus, hydrus_profile), (time_of, of_vtk), color in zip(
        hydrus_profiles.items(), of_selected_vtk_files.items(), colors
    ):
        mesh = pv.read(f"{OPENFOAM_FOLDER}/VTK/" + of_vtk)
        line = pv.Line(a := [0, 0, mesh.bounds[5]], b := [0, 0, mesh.bounds[2]])

        sample = mesh.sample_over_line(a, b)

        for var, plot in plots.items():
            ax = plot.ax

            if var == "h(z)":
                xdata_of = sample["h"]
                xdata_hy = hydrus_profile["Head"]
            elif var == "θ(h)":
                xdata_of = sample["Sw"] * sample["porosity"]
                xdata_hy = hydrus_profile["Moisture"]
            elif var == "K(h)":
                xdata_of = sample["hydraulicCond"]
                xdata_hy = hydrus_profile["K"]
                ax.set_xscale("log")
            elif var == "C(h)":
                xdata_of = sample["capillarity"] * sample["porosity"]
                xdata_hy = hydrus_profile["C"]

            ##  OpenFOAM results
            ax.plot(
                xdata_of,
                sample.points[:, 2],  # <- Depth (m)
                label=f"{time_of:.2f} hr",
                c=color,
                **plot_kwargs["openfoam"],
            )

            ##  Hydrus1D results
            ax.plot(
                xdata_hy,
                COLUMN_LENGTH + hydrus_profile["Depth"],  # <- Depth (m)
                c=color,
                **plot_kwargs["hydrus1d"],
            )

            ax.set_xlabel(rf"${var}$")
            ax.grid(True, ls="dashed", alpha=0.5)
            ax.legend(
                loc="center left", bbox_to_anchor=[1.01, 0.5], title="Time", fontsize=12
            )
            ax.set_ylabel("Depth [m]")

        # Summary in pandas format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        of_pandas_summary[time_of] = pd.DataFrame(
            {
                "z": sample.points[:, 2],
                "h": sample["h"],
                "theta": sample["Sw"] * sample["porosity"],
                "hydraulic_cond": sample["hydraulicCond"],
                "capillarity": sample["capillarity"] * sample["porosity"],
            }
        )


elif PLOT_ENGINE == "plotly":
    plots = {k: go.Figure() for k in VARIABLE_LABELS}

    for (time_hydrus, hydrus_profile), (time_of, of_vtk), color in zip(
        hydrus_profiles.items(), of_selected_vtk_files.items(), colors
    ):
        color = np.round(color * 255)
        color_rgb = f"rgb({color[0]}, {color[1]}, {color[2]})"
        color_rgba = f"rgba({color[0]}, {color[1]}, {color[2]}, 0.5)"

        mesh = pv.read(f"{OPENFOAM_FOLDER}/VTK/" + of_vtk)
        line = pv.Line(a := [0, 0, mesh.bounds[5]], b := [0, 0, mesh.bounds[2]])

        sample = mesh.sample_over_line(a, b)

        for var, fig in plots.items():
            if var == "h(z)":
                xdata_of = sample["h"]
                xdata_hy = hydrus_profile["Head"]
            elif var == "θ(h)":
                xdata_of = sample["Sw"] * sample["porosity"]
                xdata_hy = hydrus_profile["Moisture"]
            elif var == "K(h)":
                xdata_of = sample["hydraulicCond"]
                xdata_hy = hydrus_profile["K"]
                fig.update_xaxes(type="log")
            elif var == "C(h)":
                xdata_of = sample["capillarity"] * sample["porosity"]
                xdata_hy = hydrus_profile["C"]

            # Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ##  OpenFOAM results
            fig.add_trace(
                go.Scatter(
                    x=xdata_of,
                    y=sample.points[:, 2],
                    name="OpenFOAM",
                    legendgroup=f"{time_of:.2f} h",
                    legendgrouptitle_text=f"{time_of:.2f} h",
                    line=dict(color=color_rgba, width=10),
                )
            )
            ##  Hydrus1D results
            fig.add_trace(
                go.Scatter(
                    x=xdata_hy,
                    y=COLUMN_LENGTH + hydrus_profile["Depth"],
                    name="Hydrus-1D",
                    legendgroup=f"{time_of:.2f} h",
                    line=dict(color=color_rgb, width=3, dash="dash"),
                )
            )

            fig.update_yaxes(
                title=dict(text="Depth [m]", font_color="#444", font_size=30),
                showgrid=True,
                mirror=False,
                ticks="outside",
                color="#bbb",
                showline=True,
            )

            fig.update_xaxes(
                title=dict(text=var, font_color="#444", font_size=30),
                showgrid=True,
                mirror=False,
                ticks="outside",
                color="#bbb",
                showline=True,
            )

            fig.update_layout(
                margin=dict(t=10, pad=5),
                height=800,
            )

for plot, tab in zip(plots.values(), tabs):
    with tab:
        if PLOT_ENGINE == "matplotlib":
            st.pyplot(plot.fig, use_container_width=True)

        elif PLOT_ENGINE == "plotly":
            st.plotly_chart(plot, use_containter_width=True)
