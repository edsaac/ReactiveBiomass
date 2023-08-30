import streamlit as st

st.set_page_config(layout="wide")

st.session_state.css_style = """
    <style>
        [data-testid=stImage]{
            text-align: center;
            display: block;
            margin-left: auto;
            margin-right: auto;
            width: 90%;
        }

        [data-testid="stSidebar"] ul li:first-child {
            font-size: 1.2rem;
            font-style: italic;
        }
    </style>
    """

st.markdown(st.session_state.css_style, unsafe_allow_html=True)

st.session_state["axis_setup"] = dict(
    mirror=False,
    ticks="outside",
    color="#bbb",
    showline=True)

st.image(".hiddendata/column_schematic.png")

"# Rosenzweig (2011) - The effect of biofilms on the hydraulic properties of unsaturated soils"

r"""Two experiments were run: 

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

"""