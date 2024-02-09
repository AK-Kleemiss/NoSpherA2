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
df = pd.read_csv(lf, delim_whitespace=True, header=None)

plt.rcParams['lines.linewidth'] = 0.3
# Create a figure with 2x2 subplots
fig, axs = plt.subplots(2,2, figsize=(12, 8))

axs[0][0].plot(df[0], df[1], label='Thakkar')
axs[0][0].plot(df[0], df[4], label='HAR + Thakkar Core')
axs[0][0].plot(df[0], df[6], label='def2 + Thakkar Core')
axs[0][0].plot(df[0], df[7], label='DKH-Jorge ORCA')
axs[0][0].set_title('All SF')
axs[0][0].legend()

axs[1][1].plot(df[0], df[1]-df[2] - df[3], label='Thakkar - ORCA Valence')
axs[1][1].legend()

axs[0][1].plot(df[0], df[1]-df[2], label = 'Thakkar - Thakkar Core')
axs[0][1].plot(df[0], df[3], label='HAR Valence')
axs[0][1].plot(df[0], df[5], label='def2 Valence')
axs[0][1].plot(df[0], df[8], label='DKH-Jorge Valence')
axs[0][1].set_title('Valence SF')
axs[0][1].legend()

axs[1][0].plot(df[0], df[1]-df[4], label = 'Thakkar - HAR')
axs[1][0].plot(df[0], df[7]-df[4], label = 'ORCA Jorge - HAR')
axs[1][0].plot(df[0], df[7]-df[1], label = 'ORCA Jorge - Thakkar')
axs[1][0].plot(df[0], df[6]-df[1], label = 'ORCA def2 - Thakkar')
axs[1][0].plot(df[0], df[7]-df[6], label = 'ORCA Jorge - ORCA def2')
axs[1][0].legend()

fig.savefig("ECP_SF_Au.png", dpi=300, bbox_inches='tight') 