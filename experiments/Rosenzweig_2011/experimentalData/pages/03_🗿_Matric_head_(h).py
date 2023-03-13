import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import streamlit as st
import plotly.graph_objects as go


REPO_PATH = subprocess.check_output(['git', 'rev-parse', '--show-toplevel']).decode('utf-8').strip()
plt.style.use(f'{REPO_PATH}/misc/edwin.mplstyle')
st.markdown(st.session_state.css_style, unsafe_allow_html=True)

colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33']
z = np.array([0.57, 0.52, 0.47, 0.42, 0.37, 0.32])

"# 🗿 Matric head"

labels = ["Control", "Inoculated"]

for tab, dataset in zip(st.tabs([f"**{i}**" for i in labels]), labels):

    with tab:
        head = pd.read_excel(".hiddendata/matricHead.xlsx", sheet_name=dataset)
        head["Time (d)"] = head["Time (min)"]/(60*24)
        
        with st.expander("🗃️ **Dataframes**"):
            st.dataframe(head, height=200)


        with st.expander("📈 **TDR timeseries**"):
            
            fig = go.Figure()
            
            for color, zval, hcol in zip(colors, z, [k for k in head.columns if "h_Avg" in k]):

                head[hcol] = head[hcol].mask(head[hcol] > 0.0).mask(head[hcol] < -100.0)

                fig.add_trace(
                    go.Scatter(
                        x=head["Time (d)"],
                        y=head[hcol],
                        name=f"{hcol} @ {zval} m",
                        line=dict(color=color, width=3),
                    )
                )
            fig.update_yaxes(
                title_text="Matric head [cm]", 
                range=[-50,0.0], 
                showgrid=False,
                **st.session_state.axis_setup)

            fig.update_xaxes(
                title_text="Time [d]", 
                showgrid=False,
                **st.session_state.axis_setup)
            
            st.plotly_chart(fig, use_container_width=True)

        with st.expander("🌽 **Matric head spatial distribution**", expanded=True):
        
            t = head["Time (d)"].to_numpy()
            h = head[[i for i in head.columns if "h_" in i]].to_numpy()/100.

            fig = go.Figure()

            hovertemplate = """<b>h = %{z:.2f} m</b> <br> z = %{y:.2f} m <br> t = %{x:.1f} d"""

            fig.add_trace(
                go.Heatmap(
                    x=t, y=z, z=h.T,
                    zmin=-0.40, zmax=-0.10,
                    hovertemplate=hovertemplate,
                    name="h [m]",
                    colorbar=dict(
                        len=0.7,
                        lenmode="fraction",
                        thickness=0.02,
                        thicknessmode="fraction",
                        title_text="h [m]",
                        titlefont_size=18,
                        outlinewidth=0,
                        y=0.50,
                        x=1.04
                    ))
            )
            
            fig.update_xaxes(
                title_text="Time [d]",
                showgrid=False,
                **st.session_state.axis_setup)

            fig.update_yaxes(
                title_text="Depth [m]", 
                range=[0.28,0.61], 
                showgrid=False,
                **st.session_state.axis_setup)

            st.plotly_chart(fig, use_container_width=True)