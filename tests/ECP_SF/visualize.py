import pandas as pd
import matplotlib.pyplot as plt

# Read the file
df = pd.read_csv('NoSpherA2.log', delim_whitespace=True, header=None)

# Create a figure with 2 subplots
fig, axs = plt.subplots(2)

# Plot the second and third columns in the first subplot
axs[0].set_title('ECP SF')
axs[0].plot(df[0], df[1], label='Thakkar')
axs[0].plot(df[0], df[2], label='Gaussian')
axs[0].legend()

axs[1].set_title('ECP Density')
# Plot the fourth and fifth columns in the second subplot
axs[1].plot(df[0], df[3], label='Thakkar')
axs[1].plot(df[0], df[4], label='Gaussian')
axs[1].plot(df[0], df[5], label='Thakkar full')
axs[1].legend()

# Show the figure
plt.show()
exit(0)

