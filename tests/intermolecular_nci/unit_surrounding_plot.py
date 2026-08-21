from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

def read_float_arg(name, default):
    if name not in sys.argv[1:]:
        return default
    idx = sys.argv.index(name)
    if idx + 1 >= len(sys.argv):
        raise SystemExit(f"{name} requires a numeric value")
    return float(sys.argv[idx + 1])

dat_path = Path(__file__).with_name("unit_surrounding_values.dat")
data = np.loadtxt(dat_path, comments="#")
if data.ndim == 1:
    data = data.reshape(1, -1)

signed_rho = data[:, 0]
rdg = data[:, 1]
xmax = 0.5
ymax = 1
xmax = read_float_arg("-xmax", None if xmax < 0.0 else xmax)
ymax = read_float_arg("-ymax", None if ymax < 0.0 else ymax)
mask = np.ones_like(signed_rho, dtype=bool)
if xmax is not None:
    mask &= np.abs(signed_rho) <= xmax
if ymax is not None:
    mask &= rdg <= ymax
signed_rho = signed_rho[mask]
rdg = rdg[mask]
if signed_rho.size == 0:
    raise SystemExit("No points remain after applying plot limits")
rho_min = float(np.min(signed_rho))
rho_max = float(np.max(signed_rho))

cmap = LinearSegmentedColormap.from_list("nci_bgr", ["blue", "green", "red"])
if rho_min < 0.0 < rho_max:
    norm = TwoSlopeNorm(vmin=rho_min, vcenter=0.0, vmax=rho_max)
else:
    norm = None

fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
scatter = ax.scatter(signed_rho, rdg, c=signed_rho, s=4, cmap=cmap, norm=norm, linewidths=0)
ax.axvline(0.0, color="0.65", linewidth=0.8)
ax.set_xlabel("signed rho")
ax.set_ylabel("RDG")
if xmax is not None:
    ax.set_xlim(-xmax, xmax)
if ymax is not None:
    ax.set_ylim(0.0, ymax)
ax.set_title(dat_path.stem.replace("_values", ""))
cbar = fig.colorbar(scatter, ax=ax)
cbar.set_label("signed rho")
png_path = dat_path.with_name(dat_path.stem.replace("_values", "_plot") + ".png")
fig.savefig(png_path, dpi=300)
print(f"Wrote {png_path}")
if "-show" in sys.argv[1:]:
    plt.show()
else:
    plt.close(fig)
