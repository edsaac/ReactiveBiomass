import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import streamlit as st

REPO_PATH = subprocess.check_output(['git', 'rev-parse', '--show-toplevel']).decode('utf-8').strip()
plt.style.use(f'{REPO_PATH}/misc/edwin.mplstyle')

"# Water content"
st.info("TDR: Time Domain Reflectometry -> dielectric sensors that measure the charge storing capacity of the soil")

## Globals
α = 2.79 #1/m
θs = 0.385
θr = 0.012
n = 7.26
m = 1 - 1/n
Ks = 2.07E-4  #m/s
qtarget = 3.32E-6 #m/s

z = np.array([0.57, 0.52, 0.47, 0.42, 0.37, 0.32])
all_probes = {k:v for k,v in zip([f"θ{i}" for i in range(1,7)], z)}

probes = st.sidebar.multiselect(
    "Probes", 
    all_probes.keys(),
    all_probes.keys())

labels = ["Control", "Inoculated"]

for tab, dataset in zip(st.tabs(labels), labels):
    with tab:
        f"## {dataset}"

        θ = pd.read_excel(".hiddendata/TDR-theta.xlsx", sheet_name=dataset)
        θend = pd.read_excel(".hiddendata/theta-endexperiment.xlsx", sheet_name=dataset)
        
        θ["Time (d)"] = θ["Time(hr)"]/24.0
        for i in range(1,7): θ[f"θ{i}_roll"] = θ[f"θ{i}"].rolling(10, center=True).mean()
        cols = st.columns(2)
        with cols[0]: 
            "### TDR measurements"
            st.dataframe(θ, height=200)
        with cols[1]:
            "### Gravimetric end-time"
            st.dataframe(θend, height=200)

        "**TDR timeseries**"
        fig, ax = plt.subplots(figsize=[10,6])
        colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33']

        for probe, color in zip(probes, colors):
            ax.plot(θ["Time (d)"], θ[probe], lw=0.5, alpha=0.6, c=color)
            ax.plot(θ["Time (d)"], θ[probe].rolling(10).mean(), lw=3, label=f"{probe} @ {all_probes[probe]:.2f}m", marker="None", alpha=1, c=color)

        ax.legend(loc='upper right')
        ax.set_ylabel("Water content (θ)")
        ax.set_xlabel("Time (d)")
        ax.set_xlim(left=0)
        ax.set_ylim(θr,θs)
        fig.tight_layout()
        st.pyplot(fig)

        "**Water content spatial distribution**"
        t = θ["Time (d)"].to_numpy()
        θroll = θ[[f"θ{i}_roll" for i in range(2,7)]].to_numpy()
        θ = θ[["θ2","θ3","θ4","θ5","θ6"]].to_numpy()

        fig,axs = plt.subplots(1,2, 
            figsize=[10,5], 
            gridspec_kw={
                "width_ratios":[1.0,0.7],
                "wspace":0.01}, 
            sharey=True)

        ax = axs[0]
        img = ax.pcolormesh(t, z[1:], θ.T, cmap="RdPu", vmin=θr, vmax=θs)
        plt.colorbar(img, ax=ax, shrink=0.5, pad=0.01)
        ax.set_ylabel("Depth [m]")
        ax.set_xlabel("Time (d)")

        ax = axs[1]
        ## TDR readings at end-time
        θ_tmean = np.average(θ[-10:], axis=0)  ## Last 50 hr
        ax.scatter(θ_tmean, z[1:], 
            c="#4d0c69", marker="o", label="TDR")
        
        θ_dmean = np.average(θ_tmean)  ## Mean over depth
        for i,j in zip(θ_tmean, z[1:]): ax.annotate(f"{i:.3f}",[i,j], [i+0.01,j+0.005], fontsize=8)
        ax.axvline(x=θ_dmean, lw=1, alpha=0.5, ls="dotted", c="#4d0c69", 
            label=rf"$\langle\bar{{\theta}}\rangle$={θ_dmean:.3f}")

        ## Gravimetric measurement at the end of the experiment
        ax.scatter("θ", "z (m)", data=θend, 
            marker="X", label="Gravimetric", c="orange")
        for i,j in zip(θend["θ"], θend["z (m)"]): 
            ax.annotate(f"{i:.3f}",[i,j], [i+0.01,j+0.005], fontsize=8, c="orange")
        θend_dmean = θend['θ'].mean()
        ax.axvline(x=θend_dmean, lw=1, alpha=0.5, ls="dotted", c="orange",
            label=rf"$\langle\bar{{\theta}}\rangle$={θend_dmean:.3f}")

        ax.set_xlim(θr,θs)
        ax.set_xlabel(r"Μean water content ($\bar{θ}$)")
        ax.legend(fontsize=8)
        fig.tight_layout()
        st.pyplot(fig)
