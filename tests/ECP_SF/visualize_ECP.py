import os
import pandas as pd
import matplotlib.pyplot as plt
dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir)
lf = os.path.join(dir, 'sfacs.dat')
# Read the file
# 0 = k
# 1 = Thakakr FF
# 2 = Thakkar Core FF
# 3 = Orca Valence FF
# 4 = ORCA ECP FF
# 5 = ORCA_sphere Valence FF
# 6 = ORCA_sphere ECP FF
# 7 = ORCA_all FF
# 8 = ORCA_all Valence FF
# 9 = ORCA_ZORA all FF
# 10 = ORCA_ZORA Valence FF
# 11 = ORCA_X2C all FF
# 12 = ORCA_X2C Valence
df = pd.read_csv(lf, delim_whitespace=True, header=None)

plt.rcParams['lines.linewidth'] = 1.0
# Create a figure with 2x2 subplots
fig, axs = plt.subplots(2,2, figsize=(12, 8))

axs[0][0].plot(df[0], df[1], label='Thakkar')
axs[0][0].plot(df[0], df[4], label='HAR + Thakkar Core')
axs[0][0].plot(df[0], df[6], label='def2 + Thakkar Core')
axs[0][0].plot(df[0], df[7], label='DKH-Jorge ORCA')
axs[0][0].plot(df[0], df[9], label='ZORA ORCA')
axs[0][0].plot(df[0], df[11], label='X2C ORCA')
axs[0][0].set_title('All SF')
axs[0][0].legend()

axs[1][1].plot(df[0], df[1]-df[4], label = 'Thakkar - HAR', lw =0, marker='o', markersize=0.5)
axs[1][1].plot(df[0], df[7]-df[4], label = 'ORCA DKH-Jorge - HAR', lw = 0, marker = 'o', markersize=0.5)
axs[1][1].legend()

axs[0][1].plot(df[0], df[1]-df[2], label = 'Thakkar - Thakkar Core', lw = 1.0)
axs[0][1].plot(df[0], df[3], label='HAR Valence')
axs[0][1].plot(df[0], df[5], label='def2 Valence', lw = 1.0)
axs[0][1].plot(df[0], df[8], label='DKH-Jorge Valence', lw = 1.0)
axs[0][1].plot(df[0], df[10], label='ORCA_ZORA Valence', lw = 1.0)
axs[0][1].plot(df[0], df[12], label='ORCA_X2C Valence', lw = 1.0)
axs[0][1].set_title('Valence SF')
axs[0][1].legend()


axs[1][0].plot(df[0], df[7]-df[1], label = 'ORCA DKH-Jorge - Thakkar', lw = 1.0)
axs[1][0].plot(df[0], df[6]-df[1], label = 'ORCA def2 - Thakkar', lw = 1.0)
axs[1][0].plot(df[0], df[9]-df[1], label = 'ORCA ZORA - Thakkar', lw = 1.0)
axs[1][0].plot(df[0], df[11]-df[1], label = 'ORCA X2C - Thakkar', lw = 1.0)
axs[1][0].plot(df[0], df[7]-df[6], label = 'ORCA DKH-Jorge - ORCA def2', lw = 1.0)
axs[1][0].plot(df[0], df[7]-df[9], label = 'ORCA DKH-Jorge - ORCA ZORA', lw = 1.0)
axs[1][0].plot(df[0], df[7]-df[11], label = 'ORCA DKH-Jorge - ORCA X2C', lw = 1.0)
axs[1][0].legend(fontsize='x-small')

fig.savefig("ECP_SF_Au.png", dpi=600, bbox_inches='tight') 