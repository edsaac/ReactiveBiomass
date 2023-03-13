import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots

REPO_PATH = subprocess.check_output(['git', 'rev-parse', '--show-toplevel']).decode('utf-8').strip()
plt.style.use(f'{REPO_PATH}/misc/edwin.mplstyle')
st.markdown(st.session_state.css_style, unsafe_allow_html=True)

colors = ['#e41a1c','#377eb8']

"# ⚖️ Column weight & flow rate"

## import data
weight = pd.read_excel(".hiddendata/ColumnWeight.xlsx", sheet_name="Weight")
weight["Time (d)"] = weight["Time (min)"]/(60*24)
weight["Rel. Inoculated (g)"] = weight["Inoculated (g)"] - weight["Inoculated (g)"].iloc[0]
weight["Rel. Control (g)"] = weight["Control (g)"] - weight["Control (g)"].iloc[1] #Starts at 1 in the original excel file

flowrate = pd.read_excel(".hiddendata/ColumnWeight.xlsx", sheet_name="Flow rate")

with st.expander("🗃️ **Dataframes**"):
    cols = st.columns(2)
    
    with cols[0]: 
        "**Weight**"
        st.dataframe(weight, height=200)
    
    with cols[1]: 
        "**Flowrate**"
        st.dataframe(flowrate, height=200)

fig = make_subplots(specs=[[{"secondary_y": True}]])

for color, label in zip(colors, ["Inoculated", "Control"]):

    hovertemplate = """Δw = %{y:.1f} g <br> t = %{x:.1f} d"""
    
    fig.add_trace(
        go.Scatter(
            x=weight["Time (d)"],
            y=weight[f"Rel. {label} (g)"],
            name=f"Rel. weight change (g)",
            legendgroup=label,
            legendgrouptitle_text=label,
            line=dict(color=color, width=3),
            hovertemplate=hovertemplate
        ),
        secondary_y=False
    )

    for ls, flow in zip(["dash", "dot"], ["Inflow", "Outflow"]):
        fig.add_trace(
            go.Scatter(
                x=flowrate["Time (d)"],
                y=flowrate[f"{flow} {label} (mL/min)"],
                name=f"{flow} {label} (mL/min)",
                legendgroup=label,
                legendgrouptitle_text=label,
                line=dict(color=color, width=2, dash=ls),
                hovertemplate=hovertemplate
            ),
            secondary_y=True
        )

fig.update_layout(
    legend=dict(
        orientation="h",
        yanchor="bottom", y=0.99,
        xanchor="center", x=0.50
    ))

fig.update_xaxes(title_text="Time (d)", showgrid=False, **st.session_state.axis_setup)
fig.update_yaxes(title_text="Weight difference [g]", secondary_y=False, showgrid=False, **st.session_state.axis_setup)
fig.update_yaxes(title_text="Flow rate (mL/min)", secondary_y=True, showgrid=False, **st.session_state.axis_setup)


st.plotly_chart(fig, use_container_width=True)
    # "**Change in weight**"
    # fig,ax = plt.subplots()
    # ax.plot("Time (d)", "Rel. Inoculated (g)", data=weight, c="k")
    # ax.plot("Time (d)", "Rel. Control (g)", data=weight, c="cornflowerblue")
    # ax.legend()
    # ax.set_ylabel("Weight difference (g)")
    # ax.set_xlabel("Time (d)")
    # fig.tight_layout()
    # st.pyplot(fig)

    # "**Flow rate timeseries**"
    # fig,ax = plt.subplots()
    # ax.plot("Time (d)", "Inflow Inoculated (mL/min)", data=flowrate, c="k")
    # ax.plot("Time (d)", "Outflow Inoculated (mL/min)", data=flowrate, c="k", ls="dashed")

    # ax.plot("Time (d)", "Inflow Control (mL/min)", data=flowrate, c="cornflowerblue")
    # ax.plot("Time (d)", "Outflow Control (mL/min)", data=flowrate, c="cornflowerblue", ls="dashed")

    # ax.legend()
    # ax.set_ylabel("Flow rate (mL/min)")
    # ax.set_xlabel("Time (d)")
    # fig.tight_layout()
    # st.pyplot(fig)