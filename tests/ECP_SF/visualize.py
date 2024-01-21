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
exps = [
1310.720000,
 655.360000,
 327.680000,
 163.840000,
  81.820000,
  40.960000,
  20.480000,
  10.240000,
   5.120000,
   2.560000,
   1.280000,
   0.640000,
   0.320000,
   0.160000,
   0.080000,
   0.040000,
   0.020000,
   0.010000
]

def Rb_func(r, coefs, exponent = None):
    res = 0
    for _i, c in enumerate(coefs):
        if exponent is None:
            res += abs(c) * np.exp(-exps[_i] * r**2)
        else :
            res += abs(c) * np.exp(-exponent[_i] * r**2)
    return res

def get_numbers_from_filename(filename):
    # Use a regular expression to find all groups of digits in the filename
    numbers = re.findall(r'\d+', filename)
    str = ''
    str = str.join(numbers)

    # Convert the numbers to integers and return them
    return str

atom_dict = {}

def read_adf_file(lf_n):
    # Open the binary file in read mode
    with open(lf_n, 'rb') as f:
        data1 = f.read(4)
        number_vals = struct.unpack('i', data1)[0]
        at = get_numbers_from_filename(lf_n)
        atom_dict[at] = []
        for i in range(number_vals):
            # Read 4 bytes (size of a float)
            data = f.read(8)
            data2 = f.read(8)
    
            # If data is empty, we've reached the end of the file
            if not data or not data2:
                break
    
            # Unpack the bytes to a float
            number = struct.unpack('d', data)[0]
            number2 = struct.unpack('d', data2)[0]
    
            # Save the number
            atom_dict[at].append(number)
        
dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir)

def process_adf_files():
    # Directory containing the .sdf files
    directory = os.path.join(dir, 'Atomic_densities')

    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        # Check if the file is a .sdf file
        if filename.endswith('.adf'):
            # Full path to the file
            file_path = os.path.join(directory, filename)

            # Call the read_adf_file function on the file
            read_adf_file(file_path)


def write_dict_to_file(dictionary, filename):
    with open(filename, 'w') as f:
        json.dump(dictionary, f, indent=2)

def read_dict_from_file(filename):
    with open(filename, 'r') as f:
        dictionary = json.load(f)
    return dictionary

#!!This updates the json file!!
#process_adf_files()
#write_dict_to_file(atom_dict, 'atom_dict.json')
atom_dict = read_dict_from_file('atom_dict.json')


lf = os.path.join(dir, 'NoSpherA2.log')
# Read the file
df = pd.read_csv(lf, delim_whitespace=True, header=None)

#dx = 0.1
upper = 2.0
nr = 200001

x = np.linspace(0, upper, nr)
dx = x[1]
y = np.zeros_like(x)
z = np.zeros_like(x)

coeffies = atom_dict['240']
def int_check():
    inty = 0
    for i in range(len(exps)):
        inty_loc = 0
        for v in x:
            norm = pow(pow(2, 7) * pow(exps[i], 3 ) / np.pi, 0.25)
            coef = coeffies[i]
            inty_loc += v**2 * (Rb_func(v, [coef*norm], [exps[i]])**2) * x[1]
                                   #[coef*norm], [exps[i]]) * x[1]
        inty += inty_loc
        print(i, inty_loc)

    print("Final: ", inty)
    exit(0)

for i in range(len(coeffies)):
    coeffies[i] *= pow(pow(2, 7) * pow(exps[i], 3 ) / np.pi, 0.25)
for i, v in enumerate(x):
    y[i] = Rb_func(v, coeffies) #/ 4 / np.pi
    z[i] = y[i] * (v)**2 * x[1]

integral = sum(x**2 * y * x[1])
print(integral)

integral2 = sum(df[0][1]*dx * df[5] * (df[0]*dx)**2 * 4 * np.pi)
print(integral2)

# Create a figure with 2 subplots
fig, axs = plt.subplots(2,2, figsize=(12, 8))

#axs[0][0].set_title('ECP SF')
#axs[0][0].plot(df[0], df[1], label='Thakkar')
#axs[0][0].plot(df[0], df[2], label='Gaussian')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[0])], [exps[0]]), label='Gaussian1 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[1])], [exps[1]]), label='Gaussian2 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[2])], [exps[2]]), label='Gaussian3 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[3])], [exps[3]]), label='Gaussian4 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[4])], [exps[4]]), label='Gaussian5 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[5])], [exps[5]]), label='Gaussian6 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[6])], [exps[6]]), label='Gaussian7 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[7])], [exps[7]]), label='Gaussian8 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[8])], [exps[8]]), label='Gaussian9 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[9])], [exps[9]]), label='Gaussian10 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[10])], [exps[10]]), label='Gaussian11 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[11])], [exps[11]]), label='Gaussian12 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[12])], [exps[12]]), label='Gaussian13 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[13])], [exps[13]]), label='Gaussian14 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[14])], [exps[14]]), label='Gaussian15 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[15])], [exps[15]]), label='Gaussian16 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[16])], [exps[16]]), label='Gaussian17 python')
axs[0][0].plot(x,Rb_func(x, [abs(coeffies[17])], [exps[17]]), label='Gaussian18 python')
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
#axs[1][0].plot(df[0]*dx, df[5] - df[4], label='Gaussian Valence')
axs[1][0].plot(df[0]*dx, df[6] - df[7], label='ORCA Valence')
axs[1][0].plot(df[0]*dx, df[5] - df[6], label='Thakkar vs ORCA')
axs[1][0].set_ylim(-10,63)
axs[1][0].legend()

axs[1][1].set_title('ECP radial Density')
axs[1][1].plot(x, z, label='Gaussian python')
axs[1][1].plot(df[0]*dx, (df[0]*dx)**2 * df[3] * 4 *np.pi, label='Thakkar')
#axs[1][1].plot(df[0]*dx, (df[0]*dx)**2 * df[4] * 4 *np.pi, label='Gaussian')
axs[1][1].plot(df[0]*dx, (df[0]*dx)**2 * df[5] * 4 *np.pi, label='Thakkar full')
axs[1][1].plot(df[0]*dx, (df[0]*dx)**2 * df[6] * 4 *np.pi, label='ORCA')
axs[1][1].legend()

# Show the figure
plt.show()
exit(0)
fig.savefig("ECP_SF.png", dpi=300, bbox_inches='tight')
exit(0)

