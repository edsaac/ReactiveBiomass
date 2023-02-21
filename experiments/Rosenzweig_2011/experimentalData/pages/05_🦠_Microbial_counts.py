import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots

REPO_PATH = subprocess.check_output(['git', 'rev-parse', '--show-toplevel']).decode('utf-8').strip()
plt.style.use(f'{REPO_PATH}/misc/edwin.mplstyle')
colors = ['#e41a1c','#377eb8']
markers = ["octagon-open", "star-square-open-dot"]

"# Cell counts"
"(CFU/mL: Collony Formation Units/mL)"

labels = ["Control", "Inoculated"]

fig = go.Figure()

for m, color, dataset in zip(markers, colors, labels):

    cfu = pd.read_excel(".hiddendata/Live-counts.xlsx", sheet_name=dataset)

    fig.add_trace(go.Scatter(
        mode="markers",
        x=cfu["CFU/mL"], y=cfu["z (m)"],
        marker=dict(color=color, symbol=m, size=12, line_width=1),
        name=dataset,
        error_x=dict(
            type="data",
            array=cfu["stdev"])
    ))

fig.update_yaxes(
    title_text="Depth [m]", 
    range=[0.28,0.61],
    showgrid=False,
    **st.session_state.axis_setup)

fig.update_xaxes(
    title_text="🦠 CFU/mL",
    type="log",
    showgrid=False,
    **st.session_state.axis_setup)

fig.update_layout(
    legend=dict(
        orientation="h",
        yanchor="bottom", y=0.99,
        xanchor="center", x=0.50
    ))

st.plotly_chart(fig, use_container_width=True)
    # ax.errorbar("CFU/mL", "z (m)", xerr="stdev", fmt='none', elinewidth=1, data=cfu)
    # ax.scatter("CFU/mL", "z (m)", marker='o', s=100, data=cfu)
    # ax.set_ylabel("Depth (m)")
    # ax.set_ylim(0, 0.62)
    # ax.set_xlabel("CFU/mL")
    # ax.set_xscale("log")
    # ax.set_xlim
    # fig.tight_layout()
    # st.pyplot(fig)