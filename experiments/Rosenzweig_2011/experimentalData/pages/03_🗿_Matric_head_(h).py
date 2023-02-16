import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import streamlit as st

REPO_PATH = subprocess.check_output(['git', 'rev-parse', '--show-toplevel']).decode('utf-8').strip()
plt.style.use(f'{REPO_PATH}/misc/edwin.mplstyle')

"# Matric head"

labels = ["Control", "Inoculated"]

for tab, dataset in zip(st.tabs(labels), labels):

    with tab:
        head = pd.read_excel(".hiddendata/matricHead.xlsx", sheet_name=dataset)
        head["Time (d)"] = head["Time (min)"]/(60*24)
        st.dataframe(head, height=200)

        "**Tensiometer data**"
        fig,ax = plt.subplots(figsize=[8,5])
        for hcol in head.columns:
            if "h_Avg" in hcol:
                head[hcol] = head[hcol].mask(head[hcol] > 0.0).mask(head[hcol] < -100.0)
                ax.plot("Time (d)", hcol, data=head)

        ax.legend()
        ax.set_ylabel("Matric head (cm)")
        ax.set_xlabel("Time (d)")
        ax.set_ylim(-50, 0.0)
        fig.tight_layout()
        st.pyplot(fig)

        "**Matric head spatial distribution**"
        t = head["Time (d)"].to_numpy()
        z = np.array([0.57, 0.52, 0.47, 0.42, 0.37, 0.32])
        h = head[[i for i in head.columns if "h_" in i]].to_numpy()/100.

        fig,ax = plt.subplots(figsize=[10,5])
        img = ax.pcolormesh(t, z, h.T, cmap="winter_r", vmin=-0.4, vmax=-0.1)
        plt.colorbar(img, ax=ax, shrink=0.5, pad=0.01, label=r"$h$ (m)")
        ax.set_ylabel("Depth [m]")
        ax.set_xlabel("Time (d)")
        st.pyplot(fig)