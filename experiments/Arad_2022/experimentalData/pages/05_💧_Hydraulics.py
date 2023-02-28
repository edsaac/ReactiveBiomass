import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots

"# 💧 Inflow and hydraulic head"

## Globals
colors = st.session_state.colors
z = np.array([25.0, 65.0, 105.0, 145.0])/100.
probes = [*"ABCD"]
add_flood_and_dry_periods = st.session_state["add_flood_and_dry_periods"]

@st.cache_data
def read_data():
    head = pd.read_excel(
        ".hiddendata/ORP.xlsx", 
        usecols=["Time (min)", "INFLOW", "HEAD (cm)"]
    )
    head["Time (d)"] = head["Time (min)"]/(60*24)
    return head

head = read_data()

## Dataframe render
with st.expander("🗃️ **Dataframes**"):
    st.dataframe(head, height=200, use_container_width=True)

## Timeseries render
with st.expander("📈 **Timeseries**", expanded=True):

    fig = make_subplots(specs=[[{"secondary_y": True}]])

    # Raw data
    fig.add_trace(
        go.Scatter(
            x=head["Time (d)"],
            y=head[f"INFLOW"],
            line=dict(color=colors[0], width=3),
            name=f"Qin"
        ),
        secondary_y=False
    )

    fig.add_trace(
        go.Scatter(
            x=head["Time (d)"],
            y=head[f"HEAD (cm)"],
            line=dict(color=colors[1], width=3),
            name=f"Head [cm]"
        ),
        secondary_y=True
    )

    add_flood_and_dry_periods(fig)
    fig.update_xaxes(title_text="Time (d)")
    fig.update_yaxes(title_text="Water inflow [?/?]", range=[600,790], secondary_y=False)
    fig.update_yaxes(title_text="Head [cm]", range=[-0.9, 25], secondary_y=True)

    fig.update_layout(
        title="Hydraulic conditions",
        legend=dict(
            orientation="h",
            yanchor="bottom", y=0.99,
            xanchor="center", x=0.50
        )
    )

    st.plotly_chart(fig, use_container_width=True)