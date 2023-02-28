import streamlit as st

st.set_page_config(layout="wide")

st.markdown(
    """
    <style>
        [data-testid=stImage]{
            text-align: center;
            display: block;
            margin-left: auto;
            margin-right: auto;
            width: 90%;
        }
    </style>
    """, unsafe_allow_html=True
)

st.session_state["axis_setup"] = dict(
    mirror=False,
    ticks="outside",
    color="#bbb",
    showline=True)

st.image(".hiddendata/column_schematic.png")

"# Rosenzweig (2011) - The effect of biofilms on the hydraulic properties of unsaturated soils"
