import streamlit as st

st.set_page_config(layout="wide")

st.session_state.css_style = """
    <style>
        [data-testid=stImage]{
            text-align: center;
            display: block;
            margin-left: auto;
            margin-right: auto;
            width: 90%;
        }

        [data-testid="stSidebar"] ul li:first-child {
            font-size: 1.2rem;
            font-style: italic;
        }
    </style>
    """

st.markdown(st.session_state.css_style, unsafe_allow_html=True)

st.session_state["axis_setup"] = dict(
    mirror=False,
    ticks="outside",
    color="#bbb",
    showline=True)

st.image(".hiddendata/column_schematic.png")

"# Rosenzweig (2011) - The effect of biofilms on the hydraulic properties of unsaturated soils"
