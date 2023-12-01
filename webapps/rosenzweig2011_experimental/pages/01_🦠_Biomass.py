import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from plotly.graph_objs.layout import XAxis

st.set_page_config(initial_sidebar_state="collapsed")
st.markdown(st.session_state.css_style, unsafe_allow_html=True)

st.title("🦠")

## Bring functionality from main page
casepath = st.session_state.casepath
ModelResult = st.session_state.ModelResult
get_results_blob = st.session_state.get_results_blob

"### 💻 Numerical simulation results"
results = get_results_blob("biomass")

selection = st.selectbox(
    "Select biomass:", 
    results.keys(), 
    format_func=lambda x:results[x].long_name
    )

results[selection].plot_heatmap()
st.plotly_chart(results[selection].figure, use_container_width=True)

"### 💻 Post-processing options"

postproc = get_results_blob("postprocessed_biomass")

postproc["Xactive"].data = results["XAR"].data + results["XN"].data + results["XDN"].data
postproc["Xinactive"].data = results["XI"].data + results["EPS"].data
postproc["Xtotal"].data = postproc["Xactive"].data + postproc["Xinactive"].data

select_post = st.selectbox(
    "Select a postprocessed quantity:", 
    postproc.keys(), 
    format_func=lambda x:postproc[x].long_name
    )

postproc[select_post].plot_heatmap()
st.plotly_chart(postproc[select_post].figure, use_container_width=True)

"****"
"### 🏁 End of experiment results"
cfu = pd.read_excel("./experimentalData/.hiddendata/Live-counts.xlsx", sheet_name="Inoculated")

layout = go.Layout(
    height=700,
    xaxis=XAxis(
        title=dict(
            text="Simulation results <br> X<sub>AR</sub> [kg/m³]",
            font_color="#444"),
        type="log",
        showline=True,
        color="RGBA(0, 68, 255, 0.3)",
        tickcolor="RGBA(0, 68, 255, 0.3)",
        showgrid=True,
        griddash="dash",
        linewidth=1,
        gridcolor="RGBA(0, 68, 255, 0.3)"),

    xaxis2 = XAxis(
        overlaying= 'x',
        title=dict(
            text="Experiment results <br> 🧫 [CFU/mL]",
            font_color="#444"),
        type="log",
        side= 'top',
        showline=True,
        griddash="dashdot",
        color="RGBA(228, 26, 28, 0.3)",
        showgrid=True,
        gridcolor="RGBA(228, 26, 28, 0.3)",
        tickcolor="RGBA(228, 26, 28, 0.3)",
        linewidth=1
        ),
    
    yaxis=dict(
        title=dict(
            text="Depth [m]",
            font_color="#444"),
        range=[0, 0.65],
        showgrid=False,
        tickcolor="#bbb",
        showline=False
    ),
    legend=dict(
        orientation="v",
        bordercolor="gainsboro",
        borderwidth=1,
        yanchor="bottom", y=0.10,
        xanchor="center", x=0.50
    )
)

fig = go.Figure(layout=layout)

fig.add_trace(
    go.Scatter(
        mode="lines",
        x=results["XAR"].data.isel(t=-1), y=results["XAR"].data.z,
        line=dict(width=3.0, color="RGBA(0, 68, 255, 1)"),
        name="Simulations",
    )
)

fig.add_trace(
    go.Scatter(
        mode="lines+markers",
        x=cfu["CFU/mL"], y=cfu["z (m)"],
        marker=dict(color="RGBA(228, 26, 28, 1.0)", symbol="star-square-open-dot", size=10, line_width=2),
        line=dict(width=0.5),
        name="Experiments",
        xaxis='x2',
        error_x=dict(
            type="data",
            array=cfu["stdev"])
        )
    )

cols = st.columns([0.5,2,0.5])
with cols[1]: 
    st.plotly_chart(fig, use_container_width=True)