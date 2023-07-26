import streamlit as st

import matplotlib.pyplot as plt
from subprocess import check_output
import matplotlib as mpl
import numpy as np
import pyvista as pv
from myusefultools import parser
from phydrus.read import read_nod_inf
import pandas as pd
from stpyvista import stpyvista
import pyvista as pv
import subprocess

## Check if xvfb is already running on the machine
is_xvfb_running = subprocess.run(["pgrep", "Xvfb"], capture_output=True)

if is_xvfb_running.returncode == 1:
    with st.sidebar:
        st.warning("Xvfb was not running...")
    pv.start_xvfb()
else:
    with st.sidebar:
        st.info(f"Xvfb is running! \n\n`PID: {is_xvfb_running.stdout.decode('utf-8')}`")


REPO_PATH = check_output(['git', 'rev-parse', '--show-toplevel']).decode('utf-8').strip()
plt.style.use(f'{REPO_PATH}/misc/edwin.mplstyle')

"# Tests for unsaturated flow solver"

OPTIONS = ["LOAMYSAND", "THREELAYERS", "TOPCLOG"]
COLUMN_LENGTH_DICT = {k:v for k,v in zip(OPTIONS, [6, 6, 5])} # m
CASE = st.selectbox("Select test case", options=OPTIONS)

COLUMN_LENGTH = COLUMN_LENGTH_DICT[CASE]

HYDRUS_NOD_INF = f"{CASE}/HYDRUS_RESULTS/Nod_Inf.out"
hydrus_profiles = read_nod_inf(HYDRUS_NOD_INF)

# Remove time zero
del hydrus_profiles[0.0]

hydrus_times_hours = hydrus_profiles.keys()

## Sample colors for plot
cmap = mpl.cm.tab20c
colors = cmap(np.linspace(0, 1, len(hydrus_times_hours)))

OPENFOAM_FOLDER = f"{CASE}/{CASE}"
of_vtk_files = parser.getVTKList(f"{OPENFOAM_FOLDER}/VTK")
of_times_seconds = parser.getTimeList(OPENFOAM_FOLDER)
of_times_hours = np.array([float(t)/3600 for t in of_times_seconds])

of_selected_vtk_files = dict()

plotter = pv.Plotter()
of_vtk_soil_properties = parser.getVTKList(f"{OPENFOAM_FOLDER}/VTK_soilProperties")

mesh = pv.read(f"{OPENFOAM_FOLDER}/VTK_soilProperties/{of_vtk_soil_properties[0]}")
## Send to streamlit
soil_property = st.selectbox("Soil property", list(set(mesh.array_names)), index=3)
plotter.add_mesh(mesh, scalars=soil_property, cmap="bwr")
plotter.background_color = "#bbbbbb"
plotter.camera.zoom("tight")
plotter.view_xz()
plotter.camera.azimuth = 45
plotter.camera.elevation = 25
plotter.window_size = [200, 500]
stpyvista(plotter, horizontal_align="center")

with st.expander("VTK <-> Hydrus time step comparisons"):
    f"""
    ## From Hydrus-1D:
    The following times were found
    """
    st.code(hydrus_times_hours)

    f"""
    ## From OpenFOAM:
    The following times were found
    """
    st.code(of_times_hours)

    
    for vtk, time_s, time_hr in zip(of_vtk_files, of_times_seconds, of_times_hours):
        if round(float(time_hr), 2) in hydrus_times_hours:
            vtk, time_s, time_hr, "*"
            of_selected_vtk_files[round(float(time_hr), 2)] = vtk
        else:
            vtk, time_s, time_hr
    
if len(hydrus_profiles) == len(of_selected_vtk_files) == len(colors):
    st.success("1:1")
else:
    st.error("Something went wrong")

############

fig_head ,ax_head = plt.subplots()
fig_theta ,ax_theta = plt.subplots()
fig_K ,ax_K = plt.subplots()
fig_capil ,ax_capil = plt.subplots()

of_pandas_summary = dict()

plot_kwargs = {
    "openfoam": dict(lw=5, alpha=0.5, zorder=1),
    "hydrus1d": dict(lw=1.5, ls="dashed", zorder=2)
}

for (time_hydrus, hydrus_profile), (time_of, of_vtk), color in zip(hydrus_profiles.items(), of_selected_vtk_files.items(), colors):

    mesh = pv.read(f"{OPENFOAM_FOLDER}/VTK/" + of_vtk)
    line = pv.Line(
    a:=[0, 0, mesh.bounds[5]],
    b:=[0, 0, mesh.bounds[2]])
    
    sample = mesh.sample_over_line(a,b)

    # Head ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##  OpenFOAM results    
    ax_head.plot(
        sample["h"],                #<- Head  (m)
        sample.points[:, 2],        #<- Depth (m)
        label=f"{time_of:.2f} hr",
        c=color, 
        **plot_kwargs["openfoam"]
    )

    ##  Hydrus1D results    
    ax_head.plot(
        hydrus_profile["Head"]/100,  #<- Head (m)
        COLUMN_LENGTH + hydrus_profile["Depth"]/100, #<- Depth (m)
        c=color,
        **plot_kwargs["hydrus1d"]
    ) 

    # Water content ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##  OpenFOAM results    
    ax_theta.plot(
        sample["Sw"] * sample["porosity"],      #<- theta = Sw * porosity
        sample.points[:, 2],          #<- Depth (m)
        label=f"{time_of:.2f} hr",
        c=color, 
        **plot_kwargs["openfoam"]
    )

    ##  Hydrus1D results    
    ax_theta.plot(
        hydrus_profile["Moisture"],  #<- Theta (m)
        COLUMN_LENGTH + hydrus_profile["Depth"]/100, #<- Depth (m)
        c=color,
        **plot_kwargs["hydrus1d"]
    )

    # Hydraulic conduct ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##  OpenFOAM results    
    ax_K.plot(
        sample["hydraulicCond"],      #<- K (m/s)
        sample.points[:, 2],          #<- Depth (m)
        label=f"{time_of:.2f} hr",
        c=color,
        **plot_kwargs["openfoam"]
    )

    ##  Hydrus1D results    
    ax_K.plot(
        hydrus_profile["K"]/(100*3600),  #<- K (m/s)
        COLUMN_LENGTH + hydrus_profile["Depth"]/100, #<- Depth (m)
        c=color,
        **plot_kwargs["hydrus1d"]
    )

    # Capillarity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##  OpenFOAM results    
    ax_capil.plot(
        sample["capillarity"] * sample["porosity"],   #<- C(h) (1/m)
        sample.points[:, 2],                #<- Depth (m)
        label=f"{time_of:.2f} hr",
        c=color,
        **plot_kwargs["openfoam"]
    )

    ##  Hydrus1D results    
    ax_capil.plot(
        hydrus_profile["C"] * 100,  #<- C(h) (m/s)
        COLUMN_LENGTH + hydrus_profile["Depth"]/100, #<- Depth (m)
        c=color,
        **plot_kwargs["hydrus1d"]
    )

    of_pandas_summary[time_of] = pd.DataFrame(
        {
            'z': sample.points[:, 2],
            'h': sample['h'],
            'theta': sample["Sw"] * sample["porosity"], 
            'hydraulic_cond': sample['hydraulicCond'],
            'capillarity': sample['capillarity'] * sample['porosity']
        }
    )

for ax in [ax_head, ax_theta, ax_K, ax_capil]:
    ax.grid(True, ls='dashed', alpha=0.5)
    ax.legend(loc='center left', bbox_to_anchor=[1.01, 0.5], title="Time")
    ax.set_ylabel("Depth [m]")

ax_head.set_title(r"Head profile $h$ [m]")
ax_head.set_xlabel(r"Head $h$ [m]")

ax_theta.set_title(r"Water content profile $\theta$ [-]")
ax_theta.set_xlabel(r"Water content $\theta$ [-]")

ax_K.set_title(r"Hydraulic conductivity profile $K$ [m/s]")
ax_K.set_xlabel(r"Hydraulic conductivity $K$ [m/s]")
ax_K.set_xscale('log')

ax_capil.set_title(r"Capillarity $C = \dfrac{d\theta}{dh}$ [1/m]")
ax_capil.set_xlabel(r"Capillarity $C$ [1/m]")

"## Results comparison"
tabs = st.tabs(["h(z)", "θ(h)", "K(h)", "C(h)"])
for fig, tab in zip([fig_head, fig_theta, fig_K, fig_capil], tabs):
    with tab: 
        st.pyplot(fig, use_container_width=True)
