import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import streamlit as st

REPO_PATH = subprocess.check_output(['git', 'rev-parse', '--show-toplevel']).decode('utf-8').strip()
plt.style.use(f'{REPO_PATH}/misc/edwin.mplstyle')

"# Ambient conditions"

temp = pd.read_excel(".hiddendata/Temperature.xlsx", sheet_name="Temperature")
temp["Time (d)"] = temp["Time (min)"]/(60*24)
st.dataframe(temp, height=200)

"**Change in temperature**"
fig,ax = plt.subplots()
ax.plot("Time (d)", "Temperature (C)", data=temp, c="k")
ax.set_ylabel("Temperature (C)")
ax.set_xlabel("Time (d)")
ax.set_ylim(5,25)
fig.tight_layout()
st.pyplot(fig)

"**Change in relative humidity**"
fig,ax = plt.subplots()
ax.plot("Time (d)", "Rel. Humidity (%)", data=temp, c="cornflowerblue")
ax.set_ylabel("Relative humidity (%)")
ax.set_xlabel("Time (d)")
ax.set_ylim(0,100)
fig.tight_layout()
st.pyplot(fig)
