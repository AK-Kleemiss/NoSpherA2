import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

## This is just for bond lengths
#CH3 CH3 CH3 CH NH3 NH3 NH3
bonds_neutron = [
1.086,
1.092,
1.093,
1.099,
1.036,
1.050,
1.038
]

bonds_HAR = [
1.085,
1.089,
1.074,
1.100,
1.028, 
1.033,
1.018
]

bonds_TFVC = [
1.079,
1.085,
1.086,
1.102,
1.051,
1.053,
1.034
]

uncertainties = [
    0.005,
    0.006,
    0.006,
    0.004,
    0.004,
    0.004,
    0.003
]

uncertainties_HAR = [
    0.007,
    0.007,
    0.007,
    0.006,
    0.008,
    0.008,
    0.008
]

uncertainties_TFVC = [
    0.005,
    0.004,
    0.006,
    0.004,
    0.006,
    0.006,
    0.006
]

#make a bar plot using matplotlib, where uncertainties is the same for all bars, and the x-axis is the bond type (CH3, CH, NH3) and the y-axis is the bond length. Use different colors for each method (neutron, HAR, TFVC) and include a legend. Save the plot as "bond_lengths.png" with a dpi of 300.
#The first three bonds are CH3, the next one is CH, and the last three are NH3. So the x-axis should have three categories: CH3, CH, NH3, and each category should have three bars for the different methods. The y-axis should be labeled "Bond Length (Å)" and the title of the plot should be "Bond Lengths from Different Methods". The legend should indicate which color corresponds to which method.
bond_types = ["C-H1",  "C-H2", "C-H3" ,"C-H4", "N-H5", "N-H6", "N-H7"]
x = np.arange(len(bond_types))
width = 0.25

fig, ax = plt.subplots(figsize=(10, 6))
ax.bar(x - width, bonds_neutron, width, yerr=uncertainties,
         label="Neutron", color="blue", capsize=5)
ax.bar(x, bonds_HAR, width, yerr=uncertainties_HAR,
         label="HAR", color="orange", capsize=5)
ax.bar(x + width, bonds_TFVC, width, yerr=uncertainties_TFVC,
            label="TFVC", color="green", capsize=5)
ax.set_xlabel("Bond")
ax.set_ylabel("Bond Length (Å)")
ax.set_ylim(0.95, 1.15)
ax.set_title("Bond Lengths from Different Methods")
ax.set_xticks(x)
ax.set_xticklabels(bond_types)
ax.legend()
plt.tight_layout()
plt.savefig("bond_lengths.png", dpi=300)
#exit(0)

# Determine the path to the data file relative to this script
data_path = Path(__file__).resolve().parent / "densities.dat"

# Load the data using the explicit path, has a header row with column names
data = pd.read_table(data_path, sep="\s+", header=1)
data.columns = ["x", "HF", "DFT", "C_dens", "H_dens", "total_dens", "B_weight_C", "B_weight_H", "TFVC_weight_C", "TFVC_weight_H", "TFVC_weight_C_DFT", "TFVC_weight_H_DFT"]

data["x"] = pd.to_numeric(data["x"], errors="coerce")
data["TFVC_weight_H"] = pd.to_numeric(data["TFVC_weight_H"], errors="coerce")
data["TFVC_weight_C"] = pd.to_numeric(data["TFVC_weight_C"], errors="coerce")

# calculate weights for each
weight_C = data["C_dens"] / data["total_dens"]
weight_H = data["H_dens"] / data["total_dens"]
offset = 1.2
xmax = 1.8 + offset

diff = data["DFT"] - data["HF"]

# Plot the data
fig, axs = plt.subplots(2, 2, figsize=(16, 9))
axs[0][0].plot(data["x"], data["HF"], label="HF")
axs[0][0].plot(data["x"], data["DFT"], label="DFT")
axs[0][0].set_xlabel("d(C) /a.u.")
axs[0][0].set_ylabel("electron density")
axs[0][0].legend()
axs[0][0].set_xlim(-offset, xmax)

axs[0][1].axhline(0, color="gray", linestyle="--")
axs[0][1].plot(data["x"], diff, label="DFT - HF")
axs[0][1].set_xlabel("d(C) /a.u.")
axs[0][1].set_ylabel(r"$\Delta \rho$ /a.u.")
axs[0][1].legend()
axs[0][1].set_xlim(-offset, xmax)

axs[1][0].axhline(0, color="gray", linestyle="--")
axs[1][0].plot(data["x"], weight_C * diff, label="Hirshfeld")
# axs[1][0].plot(data["x"], data["B_weight_C"] * diff, label="Becke")
axs[1][0].plot(data["x"], data["TFVC_weight_C"] * diff, label="TFVC")
# axs[1][0].plot(data["x"], weight_C * data["HF"], label="C (HF)")
axs[1][0].set_xlabel("d(C) /a.u.")
axs[1][0].set_ylabel(r"$\omega_{C} \cdot \Delta \rho$ /a.u.")
axs[1][0].legend()
axs[1][0].set_xlim(-offset, xmax)

axs[1][1].axhline(0, color="gray", linestyle="--")
axs[1][1].plot(data["x"], weight_H * diff, label="Hirshfeld")
# axs[1][1].plot(data["x"], data["B_weight_H"] * diff, label="Becke")
axs[1][1].plot(data["x"], data["TFVC_weight_H"] * diff, label="TFVC")
# axs[1][1].plot(data["x"], weight_H * data["HF"], label="H (HF)")
axs[1][1].set_xlabel("d(C) /a.u.")
axs[1][1].set_ylabel(r"$\omega_{H} \cdot \Delta \rho$ /a.u.")
axs[1][1].legend()
axs[1][1].set_xlim(-offset, xmax)

plt.tight_layout()
plt.savefig("differences.png", dpi=300)
