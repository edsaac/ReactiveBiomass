import streamlit as st
from pandas import read_excel
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib as mpl
import numpy as np


r"""
## List of experiments 

|Experiment →|AS150|SW150|SW240|RW150|RW240|
|---------:|:---:|:---:|:---:|:---:|:---:|
|Surface head| - | Fig2a | Fig S1a | - | - |
|Water content| - | Fig 2b-e | Fig S1b-e | - | - | 
|Dissolved oxygen| - | Fig 3 | Fig 3 | Fig S4 | - |
|Diss. org. carbon| - | - | - | Fig 4 | Fig 4 |
|NH4+| - | - | - | Fig 4 | Fig 4 |
|TKN| - | - | - | Fig 4 | Fig 4 |
|NO3-| - | - | - | Fig 4 | Fig 4 |
"""

r"""
## Soil characterization
"""
cols = st.columns([2,1])

with cols[0]:
    r"""

    ### Size distribution (Table S1 HESS)

    |Layer|Porosity|Gravel (%)|Sand (&)|Silt (%)|Clay (%)|TOC (%)|
    |:---:|:-----:|:----:|:-----:|:----:|:----:|:-----:|
    | 1 | 0.48 | 0.7 | 93.9 | 2.4 | 3.0 | 0.25 |
    | 2 | 0.48 | 0.6 | 94.4 | 2.1 | 2.9 | 0.10 |
    | 3 | 0.47 | 0.0 | 99.2 | 0.1 | 0.7 | 0.01 |
    | 4 | 0.42 | 0.0 | 97.4 | 0.6 | 2.0 | 0.02 |
    | 5 |    - |   - |    - |   - |   - |    - |
    | 6 | 0.42 | 0.4 | 86.4 | 5.4 | 7.6 | 0.05 |

    &nbsp;

    ### Calibrated hydraulic parameters (Table S2 VZJ)
    |Layer|$\theta_s$|$\theta_r$|$\alpha$ [1/m]|$n$|$K_s$ [m/min]|
    |:---:|:-----:|:----:|:-----:|:----:|:----:|
    |Thin?|0.45| 0.070 | 10.4  | 1.08 | 8E-7 |
    | 1 | 0.31 | 0.010 | 0.912 | 2.3  | 1.7E-3  |
    | 2 | 0.27 | 0.045 | 3.28  | 2.25 | 3.8E-3  |
    | 3 | 0.22 | 0.043 | 3.28  | 2.7  | 2.37E-3 |
    | 4 | 0.28 | 0.045 | 3.28  | 5.2  | 3.8E-3  |
    | 5 | 0.28 | 0.045 | 3.28  | 3.9  | 2.8E-3  |
    | 6 | 0.29 | 0.043 | 4.88  | 2.0  | 2.25E-3 |

    """ 

with cols[1]:
    column_diameter = 0.15
    soil_layers = [0, -1.53, -2.58, -3.52, -4.62, -4.80, -6.00]
    cmap = mpl.cm.tab20c
    colors = cmap(np.linspace(0, 1, len(soil_layers) - 1))

    fig, ax = plt.subplots(figsize=[2, 6])

    for i, (yi, yj, c) in enumerate(zip(soil_layers[:-1], soil_layers[1:], colors), start=1):
        ax.add_artist(Rectangle([0, yi], column_diameter, yj, fc=c, label=f"Layer {i}"))

    ax.set_ylim(-6.1, 0.1)
    ax.legend()
    fig



r"""
## Collected data

"""

book1 = read_excel("./experimental-data/.hiddendata/redatarequestfrom10_5194hess244172020/Book1.xlsx")
book1

## Time
time = book1.index.to_numpy()
deltatime = time[1] - time[0]
time_plus_one = np.concatenate([time, [time[-1] + deltatime]])

## Depth
sensor_depths_centers = [-25, -75, -175, -275]
sensor_depths_edges = [0, -50, -125, -225, -325]

## Parameters
water_content = book1[["WC25", "WC75", "WC175", "WC275"]].to_numpy().T
surface_head = np.where(book1["SH"].to_numpy() > 0.5, book1["SH"].to_numpy(), 0.5)
dissolved_oxygen = book1[["DO25", "DO75", "DO175", "DO275"]].to_numpy().T


fig, axs = plt.subplots(
    2, 1, figsize=[20,5], sharex=True,
    gridspec_kw=dict(height_ratios=[0.2, 1.0]))
ax = axs[0]
ax.plot(time, surface_head)
ax = axs[1]
img = ax.pcolormesh(time_plus_one, sensor_depths_edges, water_content)
plt.colorbar(img, orientation='horizontal')
fig

if st.button("Save water_content"):
    data_to_save = np.vstack([time*60, surface_head/100]).T
    data_to_save
    np.savetxt("h_dataTable.csv", data_to_save, header="Time[s], Head[m]", fmt=["%i", "%.2E"], delimiter=',')

# z_depths = np.array([-25, -75, -175, -275])

# st.dataframe(book1)

# fig, ax = plt.subplots()
# ax.plot(book1["75SH"])
# st.pyplot(fig)


# fig = plt.figure(figsize=[6,9])
# X,Y = fig.get_figwidth(), fig.get_figheight()

# ax = fig.add_axes((0,0,1,1))
# ax.scatter([0]*4, sensor_depths)
# ax.set_xlim(-1, 10)
# ax.spines.top.set_visible(False)
# ax.spines.right.set_visible(False)
# ax.set_ylim(-300, 0)

# COLY = ax.get_ylim()[1] - ax.get_ylim()[0]
# print(COLY)
# # ax.set_ylim(0, 300)

# for (i,y), (k,v) in zip(enumerate(sensor_depths), water_content.items()) :
    
#     ax2 = fig.add_axes(
#         (0.20, 1+y/COLY - 0.07, 0.6, 0.14),  ## left, bottom, width, height
#         # sharex = ax2 if i > 0 else None
#     )

#     if i < 3:
#         ax2.xaxis.set_ticklabels([])
    
#     ax2.plot(v)
#     ax2.set_title(k)
#     # ax2.set_xlim(0,4500)

# # fig.patches.extend([Rectangle((0,0), 0.2, 0.5, transform=fig.transFigure, figure=fig)])
# st.pyplot(fig)


