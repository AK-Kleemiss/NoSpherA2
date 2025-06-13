import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

# Determine the path to the data file relative to this script
data_path = Path(__file__).resolve().parent / "densities.dat"

# Load the data using the explicit path
data = pd.read_table(data_path, sep="\s+", header=None)
data.columns = ["x", "HF", "DFT", "C_dens", "H_dens", "total_dens", "B_weight_C", "B_weight_H"]

# calculate weights for each
weight_C = data["C_dens"] / data["total_dens"]
weight_H = data["H_dens"] / data["total_dens"]
offset = 0.5
xmax = 2.0598014859741647 + offset

diff = data["DFT"] - data["HF"]

# Plot the data
fig, axs = plt.subplots(2, 2, figsize=(16, 9))
axs[0][0].plot(data["x"], data["HF"], label="HF")
axs[0][0].plot(data["x"], data["DFT"], label="DFT")
axs[0][0].set_xlabel("d(C) /a.u.")
axs[0][0].set_ylabel("electron density")
axs[0][0].legend()

axs[0][1].axhline(0, color="gray", linestyle="--")
axs[0][1].plot(data["x"], diff, label="DFT - HF")
axs[0][1].set_xlabel("d(C) /a.u.")
axs[0][1].set_ylabel(r"$\Delta \rho$ /a.u.")
axs[0][1].legend()
axs[0][1].set_xlim(-offset, xmax)

axs[1][0].axhline(0, color="gray", linestyle="--")
axs[1][0].plot(data["x"], weight_C * diff, label="Hirshfeld")
axs[1][0].plot(data["x"], data["B_weight_C"] * diff, label="Becke")
# axs[1][0].plot(data["x"], weight_C * data["HF"], label="C (HF)")
axs[1][0].set_xlabel("d(C) /a.u.")
axs[1][0].set_ylabel(r"$\omega_{C} \cdot \Delta \rho$ /a.u.")
axs[1][0].legend()
axs[1][0].set_xlim(-offset, xmax)

axs[1][1].axhline(0, color="gray", linestyle="--")
axs[1][1].plot(data["x"], weight_H * diff, label="Hirshfeld")
axs[1][1].plot(data["x"], data["B_weight_H"] * diff, label="Becke")
# axs[1][1].plot(data["x"], weight_H * data["HF"], label="H (HF)")
axs[1][1].set_xlabel("d(C) /a.u.")
axs[1][1].set_ylabel(r"$\omega_{H} \cdot \Delta \rho$ /a.u.")
axs[1][1].legend()
axs[1][1].set_xlim(-offset, xmax)

plt.tight_layout()
plt.savefig("differences.png", dpi=300)
