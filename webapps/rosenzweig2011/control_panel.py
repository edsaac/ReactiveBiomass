import streamlit as st

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from pathlib import Path
from typing import Literal

from myusefultools.pyopenfoam import OpenFOAM

st.session_state.cmaps = ["Greens", "Purples", "Greys", "Oranges"]
st.session_state.linecolors = ["forestgreen", "purple", "grey", "darkorange"]

st.session_state.doc_list = [1, 2, 5, 10, 20, "XX"]
st.session_state.rho_x_list = [1, 3, 10, 31, 100, 316, 1000]  #<- Values chosen so they are nice powers of 10
st.session_state.d_growth_list = [0, 1e-11, 1e-10, 1e-9, 1e-8]

# if "submitted" not in st.session_state:
#     st.session_state.submitted = False


    # if st.form_submit_button("Read file"):

# cols = st.columns([1,2])
# type_of_plot = cols[0].selectbox("Type of plot", ["Probe over time", "Profile over depth"])
# log_scale = cols[0].checkbox("log scale?")

# if type_of_plot == "Probe over time":
#     z_selected = cols[0].selectbox("Depth:", biomass.z.values)

#     fig,ax = plt.subplots(figsize=[6, 5], dpi=120)
#     for (k, var), linecolor in zip(biomass.items(), linecolors):
#         ax.plot(var.t/3600, var.sel(z=z_selected, method="nearest"), alpha=0.9, lw=3, label=k, c=linecolor)

#     if log_scale: ax.set_yscale("log")
#     ax.set_ylim(1e-5, 1e-0)
#     ax.set_title(f"z = {z_selected:.2f} m", fontproperties=dict(size=12, weight=100))
#     ax.set_ylim(bottom=0)
#     ax.set_ylabel(r"Biomass $X$ [g/L]")
#     ax.set_xlabel("Time [hr]")

#     ax.legend(
#         title="Immobile biomass\n[g/L]", 
#         title_fontproperties=dict(size=10, weight=100), 
#         loc="center left", bbox_to_anchor=(1.0, 0.5)
#     )

#     cols[1].pyplot(fig)

# elif type_of_plot == "Profile over depth":
#     t_selected = cols[0].selectbox("Time:", biomass.t.values, index=len(biomass.t.values)-1)
    
#     fig,ax = plt.subplots(figsize=[6, 5], dpi=120)
#     for (k, var), linecolor in zip(biomass.items(), linecolors):
#         ax.plot(var.sel(t=t_selected, method="nearest"), var.z, alpha=0.9, lw=3, label=k, c=linecolor)

#     if log_scale: ax.set_xscale("log")
#     ax.set_xlim(1e-5, 1e-0)
#     ax.set_ylim(0, 0.62)
#     ax.set_title(f"t = {t_selected/3600:.2f} hr", fontproperties=dict(size=12, weight=100))
#     ax.set_xlabel(r"Biomass $X$ [g/L]")
#     ax.set_ylabel("Depth [m]")

#     ax.legend(
#         title="Immobile biomass\n[g/L]", 
#         title_fontproperties=dict(size=10, weight=100), 
#         loc="center left", bbox_to_anchor=(1.0, 0.5)
#     )

#     cols[1].pyplot(fig)
