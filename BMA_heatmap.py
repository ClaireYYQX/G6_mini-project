import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import openpyxl
import pandas as pd
from matplotlib.ticker import AutoMinorLocator


pairwise = pd.read_excel("BMA.xlsx", sheet_name = "pairwise")
print(pairwise.dtypes)
pairwise = pairwise.drop(columns=['Unnamed: 0'], errors='ignore')
print(pairwise.dtypes)
colors = ["white", (0, 0, 1, 0.3), (1,0,0,0.3), "purple"]

pairwise_array = pairwise.values
cmap = mcolors.ListedColormap(colors)

fig, ax = plt.subplots()
heatmap = ax.imshow(pairwise_array, cmap=cmap, interpolation = "None")
ax.set_xticks(ticks=np.arange(len(pairwise.columns)), labels = pairwise.columns, fontsize = 15)
ax.set_yticks(ticks=np.arange(len(pairwise.columns)), labels = pairwise.columns, fontsize = 15)

bounds = np.arange(len(colors)+1)
tick_positions = (bounds[:-1] + bounds[1:]) / 2.7
tick_labels = ["Exhausted Phenotype", "T-Cell Expansion Only", "Cytokine Production Only", " Both (Effector Phenotype)"]


cbar=fig.colorbar(heatmap, ticks = tick_positions)
cbar.ax.set_yticklabels(tick_labels, fontsize = 15)

##########################################################################
#Comment out region to get rid of grid lines
minor_locator = AutoMinorLocator(2)
plt.gca().xaxis.set_minor_locator(minor_locator)
plt.gca().yaxis.set_minor_locator(minor_locator)
##########################################################################

plt.grid(color='black', linestyle='dashed', linewidth=1, which = "minor")
plt.xlabel("KO1", fontsize = 30, fontweight ="bold", labelpad = 30)
plt.ylabel("KO2", fontsize = 30, fontweight ="bold")
fig.set_size_inches(18.5, 10.5)
plt.show()


fig.savefig("Heatmap_BMA.png", dpi = 100)