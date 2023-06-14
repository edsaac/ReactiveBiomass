import streamlit as st
from pandas import read_excel
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np


r"""
|Experiment|AS150|SW150|RW150|RW240|
|:---------|:---:|:---:|:---:|:---:|
|Surface head| | Fig2a | Fig S1a | | |
|Water content| | Fig 2b-e | Fig S1b-e | | |
|Dissolved oxygen| | Fig 3 | Fig 3 | Fig S4 | |
|Diss. org. carbon| | | | Fig 4 | Fig 4 |
|NH4+| | | | Fig 4 | Fig 4 |
|TKN| | | | Fig 4 | Fig 4 |
|NO3-| | | | Fig 4 | Fig 4 |


"""

book1 = read_excel("./experimental-data/redatarequestfrom10_5194hess244172020/Book1.xlsx")
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


