import streamlit as st

st.set_page_config(initial_sidebar_state="collapsed")
st.markdown(st.session_state.css_style, unsafe_allow_html=True)

casepath = st.session_state.casepath
ModelResult = st.session_state.ModelResult
get_results_blob = st.session_state.get_results_blob

st.title("💣")

####### Simulation data #########################
"## **💻 Numerical simulation**"
particulates = get_results_blob("particulates")

selection = st.selectbox(
    "Select rate:", 
    particulates.keys(), 
    format_func=lambda x:particulates[x].long_name
    )

particulates[selection].plot_heatmap()
st.plotly_chart(particulates[selection].figure, use_container_width=True)
