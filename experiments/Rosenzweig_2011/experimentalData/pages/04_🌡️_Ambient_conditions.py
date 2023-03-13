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

"# 🌡️ Ambient conditions"

temp = pd.read_excel(".hiddendata/Temperature.xlsx", sheet_name="Temperature")
temp["Time (d)"] = temp["Time (min)"]/(60*24)
with st.expander("🗃️ **Dataframes**"):
    st.dataframe(temp, height=200)


fig = make_subplots(specs=[[{"secondary_y": True}]])

fig.add_trace(
    go.Scatter(
        x=temp["Time (d)"], y=temp["Temperature (C)"],
        name="Temperature (C)",
        line=dict(color="#da4400", width=3),
        hovertemplate= """T = %{y:.1f} C <br> t = %{x:.1f} d"""
    ),
    secondary_y=False
)

fig.add_trace(
    go.Scatter(
        x=temp["Time (d)"], y=temp["Rel. Humidity (%)"],
        name="Rel. Humidity (%)",
        line=dict(color="#0044ff", width=3),
        hovertemplate= """RH = %{y:.1f}% <br> t = %{x:.1f} d"""
    ),
    secondary_y=True
)

fig.update_yaxes(
    title_text="Temperature [C]",
    range=[0,25],
    secondary_y=False,
    showgrid=False,
    **st.session_state.axis_setup)

fig.update_yaxes(
    title_text="Rel. Humidity (%)", 
    range=[0,100],
    secondary_y=True,
    showgrid=False,
    **st.session_state.axis_setup)

fig.update_xaxes(
    title_text="Time [d]", 
    showgrid=False,
    **st.session_state.axis_setup)

fig.update_layout(
    legend=dict(
        orientation="h",
        yanchor="bottom", y=0.99,
        xanchor="center", x=0.50
    ))

st.plotly_chart(fig, use_container_width=True)