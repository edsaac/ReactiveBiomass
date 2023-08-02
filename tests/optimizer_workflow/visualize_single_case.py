import numpy as np
import pandas as pd
import os, subprocess, shutil
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pyvista as pv
from myusefultools.parser import getVTKList
import xarray as xr
import streamlit as st

from stpyvista import stpyvista
import pyvista as pv

REPO_PATH = subprocess.check_output(['git', 'rev-parse', '--show-toplevel']).decode('utf-8').strip()
plt.style.use(f'{REPO_PATH}/misc/edwin.mplstyle')

caseName = st.selectbox("Case:", ["constantHead","highFlowRate","oxygenReplenish"], index=1)

PATH_TO_VTK = f"{caseName}/VTK"
LIST_OF_VTK = getVTKList(PATH_TO_VTK)

@st.cache_data
def get_probe(field) -> np.float64:

    point = np.array([[0.04, 0.04, 0.00]])

    ## Extract VTK results (this should be done with a probe but meh)
    all_vtk_paths = [os.path.join(PATH_TO_VTK, f) for f in getVTKList(PATH_TO_VTK)]
    nTimes = len(all_vtk_paths)

    ## Initialize array to store data
    results = np.zeros([nTimes])

    ## Extract field for each vtk field
    for t,vtk in enumerate(all_vtk_paths):
        mesh = pv.read(vtk)
        probing = mesh.probe(point)
        results[t] = probing[field][0]

    return results 

@st.cache_data
def get_mean_velocity() -> np.float64:

    ## Extract VTK results (this should be done with a probe but meh)
    all_vtk_paths = [os.path.join(PATH_TO_VTK, f) for f in getVTKList(PATH_TO_VTK)]
    nTimes = len(all_vtk_paths)

    ## Initialize array to store data
    results = np.zeros([nTimes])

    ## Extract field for each vtk field
    for t,vtk in enumerate(all_vtk_paths):
        mesh = pv.read(vtk)
        integrated_volume = mesh.integrate_data()
        volume = integrated_volume['Volume'][0]
        z_velocity = integrated_volume['U'][0,-1]
        results[t] = z_velocity/volume

    return results

@st.cache_data
def get_total_biomass(field:str) -> np.float64:

    ## Extract VTK results (this should be done with a probe but meh)
    all_vtk_paths = [os.path.join(PATH_TO_VTK, f) for f in getVTKList(PATH_TO_VTK)]
    nTimes = len(all_vtk_paths)

    ## Initialize array to store data
    results = np.zeros([nTimes])

    ## Extract field for each vtk field
    for t,vtk in enumerate(all_vtk_paths):
        mesh = pv.read(vtk)
        integrated_volume = mesh.integrate_data()
        total_field = integrated_volume[field][0]
        results[t] = total_field

    return results

@st.cache_data
def getAllMeshes(field:str):
    ## Extract VTK result (this should be done with a probe but meh)
    all_vtk_paths = [os.path.join(PATH_TO_VTK, f) for f in getVTKList(PATH_TO_VTK)]
    nTimes = len(all_vtk_paths)
    times = [float(t) for t in subprocess.check_output("foamListTimes", cwd=f"./{caseName}").decode("utf-8").splitlines() if t[0:2].isnumeric()]

    ## Use dimensions from the first VTK
    mesh = pv.read(all_vtk_paths[0])
    line = pv.Line(
        a:=[0, 0, mesh.bounds[5]],
        b:=[0, 0, mesh.bounds[2]])
    sample = mesh.sample_over_line(a,b)
    nPoints = len(sample[field])

    ## Initialize array to store data
    results = np.zeros([nPoints, nTimes])

    ## Extract field for each vtk field
    for t,vtk in enumerate(all_vtk_paths):
        mesh = pv.read(vtk)
        sample = mesh.sample_over_line(a,b)
        results[:,t] = sample[field]

    data = xr.DataArray(
        results, 
        dims=("z","t"), 
        coords={
            "z": sample.points[:,2], 
            "t": times})
 
    return data


## Careful with this gradient, might need a very fine mesh 
## to solve for such an abrupt change

st.header("Field over time")

scalarList = sorted("EPS NO3 rH rN XDNp XI BAP NH4 h porosity perm_clog perm_satu XDN kdet XNp tracer POCr XARp XAR rDN Sw XN katt O2 hydraulicCond DOC".split(), key=str.casefold)

cols = st.columns(3)
with cols[0]: 
    scalarName = st.selectbox("Select field", options=scalarList, index=0)

    plotter = pv.Plotter(window_size=[150,400])
    reader = pv.get_reader(f"{PATH_TO_VTK}/{LIST_OF_VTK[0]}")
    mesh = reader.read()
    plotter.add_mesh(mesh, scalars=scalarName, cmap='bwr')
    plotter.view_isometric()
    stpyvista(plotter, horizontal_align="center")

with cols[1]: 
    with_logscale = st.checkbox("Log-scale?", False)

with cols[2]: 
    label = st.text_input("Label", value=scalarName, max_chars=50)
    units = st.text_input("Units", value="mg/L", max_chars=15)

scalar = getAllMeshes(scalarName)
#fig, (ax,cax) = plt.subplots(1,2, figsize=[8,6], gridspec_kw={"width_ratios":[5,0.2]}, sharex=False)

cols = st.columns([2,1])
fig, (cax,ax) = plt.subplots(2,1, figsize=[5,7], gridspec_kw={"height_ratios":[0.2,5]}, sharex=False)
fig.set_facecolor("#ffffff00")
igt = 4  # Ignore the first timesteps for plor
if with_logscale:
    img = ax.pcolormesh(scalar.t[igt:]/86400, scalar.z, scalar[:,igt:], cmap="copper",
                        norm=colors.LogNorm(vmin=max(scalar.min(),1.0E-8), vmax=scalar.max())
                        #norm=colors.LogNorm(vmin=scalar.min(), vmax=2.067E-4)
                        #vmin=0.000, vmax=scalar.max()
                        )
else:
    img = ax.pcolormesh(scalar.t[igt:]/86400, scalar.z, scalar[:,igt:], cmap="copper")

ax.spines.right.set_visible(False)
ax.set_xlabel("Time $t$ [d]")
ax.set_ylabel("Depth $z$ [m]")
#plt.colorbar(img, cax=cax, label=f"{scalarName} {units}")
plt.colorbar(img, cax=cax, orientation="horizontal")
cax.set_title(fr"${label}$ {units}")

with cols[0]:
    st.pyplot(fig)

# Profile at the end of the experiment
fig, ax = plt.subplots(figsize=[2.5,8.0])
ax.plot(scalar[:,-1], scalar.z)
ax.set_ylabel("Depth $z$ [m]")
ax.set_xscale(("linear","log")[with_logscale])
ax.set_xlabel(f"${label}$ [{units}]")
with cols[1]:
    st.pyplot(fig)


######################################################
st.header("Flow over time")

K_0 = 2.067E-4
biomass_labels = ["XAR","EPS","XI"]
biomasses = np.array([get_total_biomass(i) for i in biomass_labels])
velocity = get_mean_velocity()


fig, ax  = plt.subplots(figsize=(5,7))
fig.set_facecolor("#ffffff00")
ax.plot(scalar.t/86400, -velocity)
ax.axhline(y=K_0, ls="dotted", c="gray", label="$K_0$")
ax.set_xlabel("Time $t$ [d]")
ax.set_ylabel("Specific discharge $q$ [m/s]")
ax.set_ylim(bottom=0.0)
ax.ticklabel_format(axis='y', scilimits=(0,0))
ax.legend()
ax.spines.right.set_visible(False)

st.pyplot(fig)

st.header("Stacked total biomass")

biomass_colors = ["tomato", "khaki", "darkkhaki"]
biomass_labels = [r"X_{\mathsf{AR}}",r"X_{\mathsf{Inert, l}}",r"X_{\mathsf{Inert, r}}"]

fig, (rax,ax)  = plt.subplots(2,1, figsize=(8,7), gridspec_kw={"height_ratios":[1,5]}, sharex=True)
fig.set_facecolor("#ffffff00")
ax.stackplot(scalar.t/86400, *biomasses, baseline="sym",
    labels = [ rf"${{{i}}}$" for i in biomass_labels ],
    colors = biomass_colors)

ax.text(13.0, -4.5E-4, r"$X_{\mathsf{Active}}$", 
    ha='center', va='center', fontsize=18, rotation=10)

ax.text(12.5, 2.0E-4, r"$X_{\mathsf{Inert},\,\mathsf{labile}}$", 
    ha='center', va='center', fontsize=18, rotation=-5)

ax.text(15.0, 5.2E-4, r"$X_{\mathsf{Inert},\,\mathsf{recalcitrant}}$", 
    ha='center', va='center', fontsize=18, rotation=-15)

#ax.legend(loc="lower left")
ax.set_xlabel("Time $t$ [d]")
ax.set_ylabel("Total immobile biomass $X$ [g/L]")
ax.spines.right.set_visible(False)
ax.ticklabel_format(axis='y', scilimits=(0,0))
ax.set_xlim(0,20)

rax.plot(scalar.t/86400, np.gradient(np.sum(biomasses, axis=0), scalar.t/86400),
    label="Biomass generation \n rate $dX/dt$ [g/L/d]")
rax.spines.right.set_visible(False)
rax.ticklabel_format(axis='y', scilimits=(0,0))
rax.axhline(y=0, ls='dotted', c='gray', lw=0.8)
rax.legend(fontsize=8)
rax.set_xlim(0,20)
st.pyplot(fig, facecolor='#ffffff00')

st.header("Nutrient transformations")

labels = ["DOC","NO3","NH4","O2"]
nutrient_colors = ["dimgray", "darkorchid", "darkturquoise", "red"]
initConc = [10,2,1,9]
probes = [get_probe(i) for i in labels]


labels = ["DOC","NO_3^-","NH_4^+","O_2"]
fig, ax = plt.subplots(figsize=(8,7))
fig.set_facecolor("#ffffff00")
for i, probe in enumerate(probes):
    ax.plot(scalar.t/86400, probe/(initConc[i]*1.0E-3), 
        label=rf"$C_{{\mathsf{{{labels[i]}}}}}$", color=nutrient_colors[i], lw=3)

ax.set_xlabel("Time $t$ [d]")
ax.set_ylabel("Nutrient tranformation \t $C_{\mathsf{out}}/C_{\mathsf{in}}$ [-]")
ax.spines.right.set_visible(False)
ax.legend(loc="lower left")
ax.axhline(y=1.0, ls='dotted', c='gray', lw=0.8)
st.pyplot(fig)