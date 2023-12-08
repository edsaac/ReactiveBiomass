import streamlit as st
from pathlib import Path

# import colorcet as cc
# import matplotlib.pyplot as plt
# import matplotlib.colors as colors
# import matplotlib as mpl

# from typing import Literal
# from myusefultools.pyopenfoam import OpenFOAM

st.session_state.cmaps = ["Greens", "Purples", "Greys", "Oranges"]
st.session_state.linecolors = ["forestgreen", "purple", "grey", "darkorange"]

st.session_state.doc_list = [2, 5, 10, 20, 50, 100, 1000]
st.session_state.rho_x_list = [
    1,
    3,
    10,
    31,
    100,
    316,
    1000,
]  # <- Values chosen so they are nice powers of 10
st.session_state.d_growth_list = [0, 1e-11, 1e-10, 1e-9, 1e-8]
# st.session_state.hydrolysis_dict = dict(
#     Low="../../experiments/Rosenzweig_2011/fit_ravid/cases_lowHydrol",
#     High="../../experiments/Rosenzweig_2011/fit_ravid/cases_highHydrol",
# )

st.session_state.hydrolysis_dict = dict(
    High = "../../exploration/time_flipper/CASES",
    Low = "../../experiments/Rosenzweig_2011/fit_ravid/cases_lowHydrol",
)

st.session_state.experiments_path = Path(
    "../../experiments/Rosenzweig_2011/experimentalData/.hiddendata"
)

R"""Two experiments were run: 

# Rosenzweig (2011) - The effect of biofilms on the hydraulic properties of unsaturated soils

- Low flowrate: 1mL/min
- High flowrate: 30mL/min

Data from the high flowrate experiment is not the best, we will only focus on the Low flowrate case.

$$
    Q = 1 \, \mathsf{mL/min} = 1.67 \times 10^{-8} \, \mathsf{m^3/s}
$$

Column dimensions are:

- Length: 60cm
- Diameter: 9cm (-1cm plexiglass edges) = 8cm
- Cross area: 5.026E-3 m²

Initial Darcy flux
$$
    q = Q/A = \dfrac{1.67 \times 10^{-8} \, \mathsf{m^3/s} }{5.026 \times 10^{-3} \, \mathsf{m^2}} = 3.316 \times 10^{-6}  \, \mathsf{m/s}
$$

### 8.2.1 Bacteria, inoculum, growth media and sand

|Item|Description|Notes|
|:---|:----------|:----|
|**Bacteria**|*Pseudomonas Putida F1*|Strict respiratory metabolism, cell size 2x0.5 um|
|**Growth media**| Luria-Bertani (LB) broth 1:75|
|**Sand**| Caesarea sand sterilized |

&nbsp;

#### Sand properties

Unsaturated parameters extracted from Table 4.1. Hydraulic conductivity from Chapter 8, p132.

|Property|Value|Notes|
|----:|:----------|:----|
|$\theta_s$| $0.385$ | Assume $n = \theta_s$ |
|$\theta_r$| $0.012$ | $S_r = 0.0312$ |
|$n$| $7.26$ | Van-Genuchten exponent |
|$\alpha$| $2.79 \, \mathsf{m^{-1}}$ |  |
|$K_s$| $2.067 \times 10^{-4} \, \mathsf{m/s}$ | $1.24 \, \mathsf{cm/s}$ |

&nbsp;

To achieve a flow rate of 1mL/min, $h = -0.432 m$ given free drainage (grad(h) = 0, so grad(h+z) = 1).

#### Carbon removal
> TOC analyses of the feeding solution and the effluent conducted at the end of the experiment showed 89% and
> 63% removal of the substrate in the low and high flow rate experiments, respectively.
"""

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
