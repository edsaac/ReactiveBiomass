import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.interpolate import InterpolatedUnivariateSpline
import streamlit as st
import subprocess

REPO_PATH = subprocess.check_output(['git', 'rev-parse', '--show-toplevel']).decode('utf-8').strip()
plt.style.use(f'{REPO_PATH}/misc/edwin.mplstyle')

with st.sidebar:
    if st.checkbox("Add CSS", True):
        # Add some styling with CSS selectors
        with open(f"{REPO_PATH}/.streamlit/style.css") as f:
            st.markdown(f"""
            <style>
                {f.read()}
            </style>
            """, unsafe_allow_html=True)

r"# Constitutive relations for unsaturated flow"

## Globals
# α = 2.79 #1/m
# θs = 0.385
# θr = 0.012
# n = 7.26
# m = 1 - 1/n
# Ks = 2.07E-4  #m/s
# qtarget = 3.32E-6 #m/s

α = 1.433 #1/m
θs = 0.3308
θr = 0.0
n = 1.506
m = 1 - 1/n
Ks = 6.944E-5  #m/s
qtarget = 3.32E-6 #m/s

with st.sidebar:
    with st.form("Parameters:"):
        h_lowerbound = st.number_input("Capillary head lower bound ", -50.0, -1.0, -15.0, format="%.3f")
        α = st.number_input("Capillary head scale parameter α [1/m]:", 0.0, 10.0, α, format="%.3f")
        θs = st.number_input("Saturated water content θs [m³/m³]:", 0.0, 1.0, θs, format="%.4f")
        θr = st.number_input("Residual water content θs [m³/m³]:", 0.0, 1.0, θr, format="%.4f")
        n = st.number_input("Van Genuchten fit exponent [-]:", 1.1, 10.0, n, format="%.3f")
        m = 1 - 1/n
        Ks = st.number_input("Saturated hydraulic conductivity [m/s]:", 1.0E-8, 1.0E-2, Ks, format="%.4e")
        qtarget = st.number_input("Target infiltration rate [m/s]:", 1.0E-10, 1.0E-2, qtarget, format="%.4e")
        st.form_submit_button("Calculate!", help="Calculate water content functions")

def vanGenuchten(h:np.array) -> np.array:
    θe = np.where(
        h >= 0,
        1.0,
        np.power(1 + np.power(α*h*np.sign(h), n), -m)
    )
    return θe

def waterSaturation(h:np.array) -> np.array:
    θe = vanGenuchten(h)
    θ = θe * (θs - θr) + θr
    return θ

def mualemPermeability(h:np.array) -> np.array:
    θe = vanGenuchten(h)
    kr = np.sqrt(θe) * np.power(1-np.power(1-np.power(θe, 1/m), m), 2)
    return kr

def capillarity(h) -> np.array:
    x = np.where(
        h >= 0,
        0.0,
        α*m*n * np.power(α*h*np.sign(h), n-1) * np.power(1 + np.power(α*h*np.sign(h), n), -m-1)
        )
    x *= (θs - θr)
    return x

h = np.concatenate([np.linspace(h_lowerbound, -1.0E-3, 500), np.linspace(0,0.5,10)])
θe = vanGenuchten(h)
θ = waterSaturation(h)
kr = mualemPermeability(h)
Ch = capillarity(h)

tabs = st.tabs(["θe(h)", "θ(h)", "kr(h)", "K(h)", "C(h)", "Calculators"])

with tabs[0]:
    r"## $\theta_e(h)$: Rel. water retention curve"
    
    cols = st.columns(2)
    with cols[0]: 
        "### Van Genuchten model"

    with cols[1]:
        st.latex(r"\theta_e(h) = \left( 1 + \left( \alpha h \right)^n \right)^{-m}")
        st.caption("with $m = 1 - 1/n$")


    fig, ax = plt.subplots(figsize=[5,4])
    ax.plot(h,θe, label=r"$\theta_e(h)$")
    ax.set_ylabel("Relative water content \t $θ_e$")
    ax.set_xlabel("Capillary head \t $h$ [m]")
    ax.set_ylim(-0.02, 1.02)
    ax.legend()
    st.pyplot(fig)


with tabs[1]:
    r"## $\theta(h)$: Water content curve"
    cols = st.columns(2)
    with cols[0]:
        st.caption(r"Bound the water content around residual ($\theta_r$) and saturation values ($\theta_s$)")
    with cols[1]:
        st.latex(r"\theta = \theta_e \left( \theta_s - \theta_r \right) + \theta_r")

    fig, ax = plt.subplots(figsize=[5,4])
    ax.plot(h,θ, label=r"$\theta(h)$")
    ax.axhline(y=θs, ls='dashed', lw=1, c='darkorange', label=fr"$\theta_s={θs}$")
    # ax.axhline(y=θfc, label=rf"$\theta_{{fc}}=${θfc:.3f}", ls=":", lw=1, c="purple")
    ax.axhline(y=θr, ls='dashed', lw=1, c='orange', label=fr"$\theta_r={θr}$")
    ax.set_ylabel("Water content (θ)")
    ax.set_xlabel("Capillary head \t $h$ [m]")
    ax.legend()
    st.pyplot(fig)


with tabs[2]:
    r"## $k_r(h)$: Relative permeability curve "
    
    cols = st.columns(2)
    with cols[0]:
        r"### Mualem model"
    with cols[1]:
        st.latex(r"k_r = \sqrt{\theta_e} \left( 1 - \left( 1 - \left( \theta_e \right)^{1/m} \right) \right)^2")

    fig, ax = plt.subplots(figsize=[5,4])

    ax.plot(h, kr, lw=2, c='darkgreen', label="$k_r$")
    # ax.scatter([h_krfc],[krfc], label=f"Field capacity\n$h=${h_krfc:.3f} m", c="purple")
    # ax.axhline(y=krfc, ls=":", lw=1, c="purple")
    ax.set_yscale('log')
    ax.set_ylabel("Rel. permeability $k_r(h)$")
    ax.set_xlabel("Capillary head \t $h$ [m]")
    ax.legend()
    st.pyplot(fig)

with tabs[3]:
    
    r"## $K(h)$: Hydraulic conductivity "
    
    cols = st.columns(2)
    with cols[0]:
        st.latex("K(h) = K_s k_r(h)")
    with cols[1]:
        st.caption(r"With $K_s$ [m/s] as the hydraulic conductivity under saturated conditions.")

    fig, ax = plt.subplots(figsize=[5,4])
    ax.plot(h, Ks*kr, lw=2, c='darkgreen', label="$K(h)$ [m/s]")
    ax.axhline(y=Ks, ls=":", lw=1, c="purple")
    ax.set_yscale('log')
    ax.set_ylabel("Hydraulic conductivity $K(h)$")
    ax.set_xlabel("Capillary head \t $h$ [m]")
    ax.legend()
    st.pyplot(fig)

with tabs[4]:
    r"## $C(h)$: Specific moisture capacity function "

    st.latex(
        r"""
        \begin{array}{rl}
        C(h) =& \dfrac{d\theta}{dh} \\ 
        \\
        =& \left( \theta_s - \theta_r \right)  \dfrac{d\theta_e}{dh} \\
        \\
        =& - \left( \theta_s - \theta_r \right) \left( \alpha m n \left( \alpha h \right)^{n-1} \left( 1 + \left( \alpha h \right)^n \right)^{-m-1} \right)
        \end{array}""")

    fig, ax = plt.subplots(figsize=[5,4])
    ax.plot(h, np.gradient(θ,h), lw=2.5, c='red', alpha=0.5, label="Numerical")
    ax.plot(h, Ch, lw=1, c='blue', label="Analytical")
    ax.set_ylabel("Capillarity \t $dθ/dh$ [1/m]")
    ax.set_xlabel("Capillary head \t $h$ [m]")
    ax.legend()
    st.pyplot(fig)

# ax = axs[3]
# ax.plot(h, qin, lw=2, c='darkgreen', label="$k_r$")
# ax.scatter([h_qtarget],[qtarget], label=f"Top bc\n$h=${h_qtarget:.3f} m", c="pink", marker="X", s=100)
# # ax.axhline(y=krfc, ls=":", lw=1, c="purple")
# ax.set_yscale('log')
# ax.set_ylabel("Infiltration rate ($q_{in}$)")
# ax.set_xlabel("Capillary head (-h)")
# ax.legend()

# for ax in axs:
#     ax.axvline(x=h_krfc, ls=":", lw=1, c="purple")

# st.pyplot(fig)


with tabs[5]:
    
    with st.expander("Given θe, solve for h"):

        cols = st.columns([1,2])
        with cols[0]:
            θe_goal = st.number_input("Goal θe", 0.001, 1.00, 0.05, format="%.2f")
        with cols[1]:
            h_seekθe = h[np.argmin(np.abs(θe - θe_goal))]
            st.latex(fr"\theta_e(h = {h_seekθe:.3f} \mathsf{{ m}}) = {vanGenuchten(h_seekθe):.3f}")

        fig, ax = plt.subplots(figsize=[5,4])
        ax.plot(h, θe, alpha=0.5)
        ax.set_ylabel("Relative water content \t $θ_e$")
        ax.set_xlabel("Capillary head \t $h$ [m]")
        ax.scatter([h_seekθe],[θe_goal], c="purple", marker="X", s=100, 
            label=f"$h$ = {h_seekθe:.3f}" + "\n" + fr"$\theta$ = {θe_goal:.3f}")
        
        ax.legend()
        st.pyplot(fig)

    with st.expander("Given kr, solve for h"):

        cols = st.columns([1,2])
        with cols[0]:
            kr_goal = st.number_input("Goal kr", 1.0E-6, 1.00, 1.0E-4, format="%.2e")
        with cols[1]:
            h_seekkr = h[np.argmin(np.abs(np.log10(kr) - np.log10(kr_goal)))]
            st.latex(fr"k_r(h = {h_seekkr:.3f} \mathsf{{ m}}) = {mualemPermeability(h_seekkr):.2E}")

        fig, ax = plt.subplots(figsize=[5,4])

        ax.plot(h, kr, alpha=0.5, label="$k_r$")
        ax.scatter([h_seekkr],[kr_goal], c="purple", marker="X", s=100)
        ax.set_yscale('log')
        ax.set_ylabel("Rel. permeability $k_r(h)$")
        ax.set_xlabel("Capillary head \t $h$ [m]")
        ax.legend()
        st.pyplot(fig)

    with st.expander("Given K, solve for h"):

        cols = st.columns([1,2])
        with cols[0]:
            K_goal = st.number_input("Goal K [m/s]", 1.0E-12, 1.00E-2, 1.0E-9, format="%.2e")
        with cols[1]:
            h_seekK = h[np.argmin(np.abs(np.log10(kr*Ks) - np.log10(K_goal)))]
            st.latex(fr"K(h = {h_seekK:.3f} \mathsf{{ m}}) = {mualemPermeability(h_seekK)*Ks:.2E}")

        fig, ax = plt.subplots(figsize=[5,4])
        ax.plot(h, kr*Ks, alpha=0.5, label="$k_r$")
        ax.scatter([h_seekK],[K_goal], c="purple", marker="X", s=100)
        ax.set_yscale('log')
        ax.set_ylabel("Hydraulic conductivity $K(h)$")
        ax.set_xlabel("Capillary head \t $h$ [m]")
        ax.legend()
        st.pyplot(fig)


"""### Quick calculator """
head_input = st.number_input("$h$ [m]", h_lowerbound, 0.0, -7.297, 0.1, format="%.3f")

rf"""
$\theta_e(h={head_input:.2f} \, \textrm{{m}}) = $ {vanGenuchten(head_input):.4f}

$\theta(h={head_input:.2f} \, \textrm{{m}}) = $ {waterSaturation(head_input):.4f}

$kr(h={head_input:.2f} \, \textrm{{m}}) = $ {mualemPermeability(head_input):.4e}

$K(h={head_input:.2f} \, \textrm{{m}}) = $ {mualemPermeability(head_input)*Ks:.3e} \, \textrm{{m/s}}

$C(h={head_input:.2f} \, \textrm{{m}}) = $ {capillarity(head_input):.3e} \, \textrm{{1/m}}
"""
    # cols = st.columns([1,2])
    # with cols[0]:
    #     kr_goal = st.number_input("Goal kr", 1.0E-6, 1.00, 1.0E-4, format="%.2e")
    # with cols[1]:
    #     h_seekkr = fsolve(lambda x: (np.log10(mualemPermeability(x)) - np.log10(kr_goal)), x0 = -10.0)
    #     st.latex(fr"h = {h_seekkr[0]:.3f} \mathsf{{ m}}")


## 2) such that K(h) = 0.01 cm/d = 1.16E-9 m/s [Twarakavi(2009)]
# Kfc = 1.16E-9 
# krfc = Kfc/Ks
# h_krfc = h[np.argmin(np.abs(kr-krfc))]
# θfc = waterSaturation(h_krfc)
# hmin_Twarakavi2009 = fsolve(lambda x: (np.log10(mualemPermeability(x)) - np.log10(krfc)), x0 = -10.0)

# print(f"{hmin_Twarakavi2009[0]=:.5f}")
# print(f"{hmin_twoPercent[0]=:.5f}")

# hmin = min([hmin_twoPercent[0], hmin_Twarakavi2009[0], -10.0])

# qin = Ks * kr
# h_qtarget = h[np.argmin(np.abs(qin - qtarget))]
