import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the file
df = pd.read_csv('NoSpherA2.log', delim_whitespace=True, header=None)

x = np.linspace(0, 10, 1001)
def Rb_func(r):
    #r /= 0.52917721092
    s = 89.5001980 * r**0 * np.exp(-5.0365510 * r**2) \
       + 0.4937610 * r**0 * np.exp(-1.9708490 * r**2) \
       +12.3169000 * r**0 * np.exp(-3.8431140 * r**2)
    p = 58.5689740 * r**1 * np.exp(-4.2583410 * r**2) \
       + 0.4317910 * r**1 * np.exp(-1.4707090 * r**2) \
       +12.3169000 * r**1 * np.exp(-3.8431140 * r**2)  
    d = 26.2248980 * r**2 * np.exp(-3.0231270 * r**2) \
       + 0.9628390 * r**2 * np.exp(-0.6503830 * r**2) \
       +12.3169000 * r**2 * np.exp(-3.8431140 * r**2)  
    f =-12.3169000 * r**3 * np.exp(-3.8431140 * r**2)  
    return s**2 + p**2 + d**2 + f**2
y = Rb_func(x)
integral = sum(x**2 * y * x[1])
print(integral)



# Create a figure with 2 subplots
fig, axs = plt.subplots(3)

# Plot the second and third columns in the first subplot
axs[0].set_title('ECP SF')
axs[0].plot(df[0], df[1], label='Thakkar')
axs[0].plot(df[0], df[2], label='Gaussian')
axs[0].legend()

axs[1].set_title('ECP Density')
# Plot the fourth and fifth columns in the second subplot
axs[1].plot(df[0], df[3], label='Thakkar')
axs[1].plot(df[0], df[4]*4, label='Gaussian')
axs[1].plot(df[0], df[5], label='Thakkar full')
axs[1].plot(x, y, label='Gaussian python')
axs[1].legend()

axs[2].plot(df[0], df[5] - df[3], label='Thakkar Difference')
axs[2].plot(df[0], df[5] - df[4], label='Gaussian Difference')
axs[2].legend()

# Show the figure
plt.show()
exit(0)

