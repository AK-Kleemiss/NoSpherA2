import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import struct
import re
import scipy.special as sp
import json
dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir)
a0 = 0.52917721067

lf = os.path.join(dir, 'sfacs.dat')
# Read the file
# 0 = k
# 1 = Thakakr FF
# 2 = Thakkar Core FF
# 3 = Orca Valence FF
# 4 = ORCA ECP FF
df = pd.read_csv(lf, delim_whitespace=True, header=None)

plt.rcParams['lines.linestyle'] = ''
plt.rcParams['lines.marker'] = '.'
plt.rcParams['lines.markersize'] = 0.5
# Create a figure with 2x2 subplots
fig, axs = plt.subplots(2,2, figsize=(12, 8))

axs[0][0].plot(df[0], df[1], label='Thakkar')
axs[0][0].plot(df[0], df[4], label='ORCA Full')
axs[0][0].legend()

axs[1][1].plot(df[0], df[1]-df[2] - df[3], label='Thakkar - ORCA Valence')
axs[1][1].legend()

axs[0][1].plot(df[0], df[1]-df[2], label = 'Thakkar - Thakkar Core')
axs[0][1].plot(df[0], df[3], label='ORCA Valence')
axs[0][1].legend()

axs[1][0].plot(df[0], df[1]-df[4], label = 'Thakkar - ORCA')
axs[1][0].legend()

fig.savefig("ECP_SF_Au.png", dpi=300, bbox_inches='tight')
exit(0)

