import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import streamlit as st

REPO_PATH = subprocess.check_output(['git', 'rev-parse', '--show-toplevel']).decode('utf-8').strip()
plt.style.use(f'{REPO_PATH}/misc/edwin.mplstyle')

cols = st.columns(2)

with cols[0]:
    "# Column weight"

    weight = pd.read_excel(".hiddendata/ColumnWeight.xlsx", sheet_name="Weight")
    weight["Time (d)"] = weight["Time (min)"]/(60*24)
    weight["Rel. Inoculated (g)"] = weight["Inoculated (g)"] - weight["Inoculated (g)"].iloc[0]
    weight["Rel. Control (g)"] = weight["Control (g)"] - weight["Control (g)"].iloc[1] #Starts at 1 in the original excel file
    st.dataframe(weight, height=200)


    "**Change in weight**"
    fig,ax = plt.subplots()
    ax.plot("Time (d)", "Rel. Inoculated (g)", data=weight, c="k")
    ax.plot("Time (d)", "Rel. Control (g)", data=weight, c="cornflowerblue")
    ax.legend()
    ax.set_ylabel("Weight difference (g)")
    ax.set_xlabel("Time (d)")
    fig.tight_layout()
    st.pyplot(fig)

with cols[1]:
    "# Flow rate"
    flowrate = pd.read_excel(".hiddendata/ColumnWeight.xlsx", sheet_name="Flow rate")
    st.dataframe(flowrate, height=200)
    
    "**Flow rate timeseries**"
    fig,ax = plt.subplots()
    ax.plot("Time (d)", "Inflow Inoculated (mL/min)", data=flowrate, c="k")
    ax.plot("Time (d)", "Outflow Inoculated (mL/min)", data=flowrate, c="k", ls="dashed")

    ax.plot("Time (d)", "Inflow Control (mL/min)", data=flowrate, c="cornflowerblue")
    ax.plot("Time (d)", "Outflow Control (mL/min)", data=flowrate, c="cornflowerblue", ls="dashed")

    ax.legend()
    ax.set_ylabel("Flow rate (mL/min)")
    ax.set_xlabel("Time (d)")
    fig.tight_layout()
    st.pyplot(fig)