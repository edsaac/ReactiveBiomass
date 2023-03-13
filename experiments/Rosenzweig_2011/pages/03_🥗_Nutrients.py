import pandas as pd
import numpy as np
import streamlit as st

st.set_page_config(initial_sidebar_state="collapsed")
st.markdown(st.session_state.css_style, unsafe_allow_html=True)

st.title("🥗")

casepath = st.session_state.casepath
ModelResult = st.session_state.ModelResult
get_results_blob = st.session_state.get_results_blob


####### Simulation data #########################
"## **💻 Numerical simulation**"
results = get_results_blob("chemical_species")

selection = st.selectbox(
    "Select biomass:", 
    results.keys(), 
    format_func=lambda x:results[x].long_name
    )

results[selection].plot_heatmap()
st.plotly_chart(results[selection].figure, use_container_width=True)
