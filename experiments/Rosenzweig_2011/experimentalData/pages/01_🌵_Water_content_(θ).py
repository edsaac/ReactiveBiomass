import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots

REPO_PATH = subprocess.check_output(['git', 'rev-parse', '--show-toplevel']).decode('utf-8').strip()
plt.style.use(f'{REPO_PATH}/misc/edwin.mplstyle')

"# 🌵 Water content (θ)"
st.info("TDR: Time Domain Reflectometry -> dielectric sensors that measure the charge storing capacity of the soil")

## Globals
α = 2.79 #1/m
θs = 0.385
θr = 0.012
n = 7.26
m = 1 - 1/n
Ks = 2.07E-4  #m/s
qtarget = 3.32E-6 #m/s

colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33']
z = np.array([0.57, 0.52, 0.47, 0.42, 0.37, 0.32])
all_probes = {k:v for k,v in zip([f"θ{i}" for i in range(1,7)], z)}

probes = [f"θ{i}" for i in range(1,7)]

labels = ["Control", "Inoculated"]

for tab, dataset in zip(st.tabs([f"**{i}**" for i in labels]), labels):
    
    with tab:

        θ = pd.read_excel(".hiddendata/TDR-theta.xlsx", sheet_name=dataset)
        θend = pd.read_excel(".hiddendata/theta-endexperiment.xlsx", sheet_name=dataset)
        
        θ["Time (d)"] = θ["Time(hr)"]/24.0
        for i in range(1,7): 
            θ[f"θ{i}_roll"] = θ[f"θ{i}"].rolling(10, center=True).mean()
        
        with st.expander("🗃️ **Dataframes**"):
            cols = st.columns(2)
            with cols[0]: 
                "### TDR measurements"
                st.dataframe(θ, height=200)
            with cols[1]:
                "### Gravimetric end-time"
                st.dataframe(θend, height=200)

        with st.expander("📈 **TDR timeseries**"):

            fig = go.Figure()
            for probe, color in zip(probes, colors):

                # Raw data
                fig.add_trace(
                    go.Scatter(
                        x=θ["Time (d)"], 
                        y=θ[probe],
                        legendgroup=probe,
                        legendgrouptitle_text=probe + f" @ {all_probes[probe]:.2f}m",
                        line=dict(color=color, width=1),
                        name="Readings"
                    )
                )

                #Rolling average
                fig.add_trace(
                    go.Scatter(
                        x=θ["Time (d)"], 
                        y=θ[probe].rolling(10).mean(), 
                        legendgroup=probe,
                        legendgrouptitle_text=probe + f" @ {all_probes[probe]:.2f}m",
                        line=dict(color=color, width=3),
                        name="Rolling average"
                    )
                )
                
            fig.update_xaxes(title_text="Time (d)")
            fig.update_yaxes(title_text="Water content (θ)", range=[θr,θs])

            st.plotly_chart(fig, use_container_width=True)

        with st.expander("🌽 **Water content spatial distribution**", expanded=True):
        
            t = θ["Time (d)"].to_numpy()
            θroll = θ[[f"θ{i}_roll" for i in range(1,7)]].to_numpy()
            θ = θ[[f"θ{i}" for i in range(1,7)]].to_numpy()
            θ_tmean = np.average(θ[-10:], axis=0)  ## Last hours

            fig = make_subplots(rows=1, cols=2, column_widths=[0.7, 0.3])

            hovertemplate = """<b>θ = %{z:.2f}</b> <br> z = %{y:.2f} m <br> t = %{x:.1f} d"""

            fig.add_trace(
                go.Heatmap(
                    x=t, y=z, z=θ.T,
                    zmin=0.01, zmax=0.35,
                    hovertemplate=hovertemplate,
                    name="θ",
                    colorbar=dict(
                        len=0.7,
                        lenmode="fraction",
                        thickness=0.02,
                        thicknessmode="fraction",
                        title_text="θ",
                        titlefont_size=18,
                        outlinewidth=0,
                        y=0.30,
                        x=1.10
                    )),
                row=1, col=1,
            )
            
            fig.update_xaxes(
                title_text="Time [d]", 
                showgrid=False,
                **st.session_state.axis_setup,
                row=1, col=1)

            fig.update_yaxes(
                title_text="Depth [m]", 
                range=[0.28,0.61], 
                showgrid=False,
                **st.session_state.axis_setup,
                row=1, col=1)
            
            # Final averaged TDR data

            hovertemplate = """θ = %{x:.2f} m <br> z = %{y:.2f} m"""

            fig.add_trace(
                go.Scatter(
                    mode="markers",
                    x=θ_tmean, y=z,
                    marker=dict(color="#6495ED", size=5),
                    name="TDR data",
                    hovertemplate=hovertemplate
                ),
                row=1, col=2
            )
            
            # Final gravimetric samples
            fig.add_trace(
                go.Scatter(
                    mode="markers",
                    x=θend["θ"], y=θend["z (m)"],
                    marker=dict(color="#BDB76B", size=5, symbol="x-dot"),
                    name="Gravimetric",
                    hovertemplate=hovertemplate
                ),
                row=1, col=2
            )

            fig.update_xaxes(
                title_text="Water content (θ)", 
                showgrid=False,
                range=[θr,θs],
                **st.session_state.axis_setup,
                row=1, col=2)

            fig.update_yaxes( 
                range=[0.28,0.61], 
                showgrid=True,
                **st.session_state.axis_setup,
                row=1, col=2)

            st.plotly_chart(fig, use_container_width=True)