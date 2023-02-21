import streamlit as st

st.set_page_config(layout="wide")

st.session_state["axis_setup"] = dict(
    mirror=False,
    ticks="outside",
    color="#bbb",
    showline=True)

cols = st.columns(2)

with cols[0]:
    st.image(".hiddendata/column_schematic.png")

with cols[1]:
    "# Rosenzweig (2011) - The effect of biofilms on the hydraulic properties of unsaturated soils"
