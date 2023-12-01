import streamlit as st

st.set_page_config(initial_sidebar_state="collapsed")
st.markdown(st.session_state.css_style, unsafe_allow_html=True)

st.title("💱")

casepath = st.session_state.casepath
ModelResult = st.session_state.ModelResult
get_results_blob = st.session_state.get_results_blob


####### Simulation data #########################
"## **💻 Numerical simulation**"
rates = get_results_blob("rates")

selection = st.selectbox(
    "Select rate:", 
    rates.keys(), 
    format_func=lambda x:rates[x].long_name
    )

rates[selection].plot_heatmap()
st.plotly_chart(rates[selection].figure, use_container_width=True)
