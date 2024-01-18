import os

import struct
dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir)
lf = os.path.join(dir, 'Rb_adf.adf')

# Open the binary file in read mode
with open(lf, 'rb') as f:
   data1 = f.read(4)
   number = struct.unpack('i', data1)[0]
   print(number)
   counter = 0
   while True:
      # Read 4 bytes (size of a float)
      data = f.read(8)

      # If data is empty, we've reached the end of the file
      if not data:
         break

      # Unpack the bytes to a float
      number = struct.unpack('d', data)[0]

      # Print the number
      print(number)
      counter += 1
   print("number: ", counter)

exit(0)

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir)
lf = os.path.join(dir, 'NoSpherA2.log')

# Read the file
df = pd.read_csv(lf, delim_whitespace=True, header=None)

dx = 0.1
nr = 20001
upper = (nr-1) * dx / 1000
a0 = 0.52917721067
x = np.linspace(0, upper, nr)
def Rb_func(r):
   s = 89.5001980 * r**0 * np.exp(-5.0365510*4*np.pi * r**2) \
      + 0.4937610 * r**0 * np.exp(-1.9708490*4*np.pi * r**2) \
      +12.3169000 * r**0 * np.exp(-3.8431140*4*np.pi * r**2)
   p = 58.5689740 * r**1 * np.exp(-4.2583410*4*np.pi * r**2) \
      + 0.4317910 * r**1 * np.exp(-1.4707090*4*np.pi * r**2) \
      +12.3169000 * r**1 * np.exp(-3.8431140*4*np.pi * r**2)  
   d = 26.2248980 * r**2 * np.exp(-3.0231270*4*np.pi * r**2) \
      + 0.9628390 * r**2 * np.exp(-0.6503830*4*np.pi * r**2) \
      +12.3169000 * r**2 * np.exp(-3.8431140*4*np.pi * r**2)  
   f =-12.3169000 * r**3 * np.exp(-3.8431140*4*np.pi * r**2)  
   return s**2 + p**2 + d**2 + f**2
y = 2 * Rb_func(x) / 4 / np.pi
integral = sum(x**2 * y * x[1] * 4 * np.pi)
print(integral)

integral2 = sum(df[0][1]*dx * df[5] * (df[0]*dx)**2 * 4 * np.pi)
print(integral2)


# Create a figure with 2 subplots
fig, axs = plt.subplots(2,2, figsize=(12, 8))

axs[0][0].set_title('ECP SF')
axs[0][0].plot(df[0], df[1], label='Thakkar')
axs[0][0].plot(df[0], df[2], label='Gaussian')
axs[0][0].legend()

axs[0][1].set_title('ECP Density')
axs[0][1].plot(x, y, label='Gaussian python')
axs[0][1].plot(df[0]*dx, df[3], label='Thakkar')
axs[0][1].plot(df[0]*dx, df[4], label='Gaussian')
axs[0][1].plot(df[0]*dx, df[5], label='Thakkar full')
axs[0][1].plot(df[0]*dx, df[6], label='ORCA')
axs[0][1].legend()

axs[1][0].set_title('Valence Density')
axs[1][0].plot(df[0]*dx, df[5] - df[3], label='Thakkar Valence')
axs[1][0].plot(df[0]*dx, df[5] - df[4], label='Gaussian Valence')
axs[1][0].plot(df[0]*dx, df[6] - df[7], label='ORCA Valence')
axs[1][0].plot(df[0]*dx, df[5] - df[6], label='Thakkar vs ORCA')
axs[1][0].set_ylim(-10,63)
axs[1][0].legend()

axs[1][1].set_title('ECP radial Density')
axs[1][1].plot(x, x**2 * y * x[1], label='Gaussian python')
axs[1][1].plot(df[0]*dx, (df[0]*dx)**2 * df[3] * 4 *np.pi * df[0][1]*dx, label='Thakkar')
axs[1][1].plot(df[0]*dx, (df[0]*dx)**2 * df[4] * 4 *np.pi * df[0][1]*dx, label='Gaussian')
axs[1][1].plot(df[0]*dx, (df[0]*dx)**2 * df[5] * 4 *np.pi * df[0][1]*dx, label='Thakkar full')
axs[1][1].plot(df[0]*dx, (df[0]*dx)**2 * df[6] * 4 *np.pi * df[0][1]*dx, label='ORCA')
axs[1][1].legend()

# Show the figure
plt.show()
exit(0)
fig.savefig("ECP_SF.png", dpi=300, bbox_inches='tight')
exit(0)

