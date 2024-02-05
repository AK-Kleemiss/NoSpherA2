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

exps  = [5.0365510000,
         1.9708490000,
         3.8431140000,
         4.2583410000,
         1.4707090000,
         3.8431140000,
         3.0231270000,
         0.6503830000,
         3.8431140000,
         3.8431140000
]

coefs = [ 89.5001980000,
           0.4937610000,
          12.3169000000,
          58.5689740000,
           0.4317910000,
          12.3169000000,
          26.2248980000,
           0.9628390000,
          12.3169000000,
         -12.3169000000
]

rad_exps = [0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0
]

def Rb_func(r, coefs = None, exponent = None, rad_exp = None):
    res = 0
    if exponent == None:
        raise ValueError("exponent is None")
    if coefs == None:
        raise ValueError("coefs is None")
    if rad_exp == None:
        raise ValueError("rad_exp is None")
    for _i, c in enumerate(coefs):
        p1 = c * r**rad_exp[_i]
        p2 = np.exp(-exponent[_i] * r**2)
        res += p1 * p2
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
    # Check if the file exists and raise an error otherwise
    if not os.path.isfile(lf_n):
        raise ValueError('File {} does not exist'.format(lf_n))
    # Open the binary file in read mode
    with open(lf_n, 'rb') as f:
        data1 = f.read(4)
        number_vals = struct.unpack('i', data1)[0]
        print(data1)
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
            print(number,number2)
    
            # Save the number
            atom_dict[at].append(number)

dir = os.path.dirname(os.path.realpath(__file__))            
os.chdir(dir)

orca_read = json.load(open(os.path.join(dir,'Atomic_densities','atom_atom37.json')))

def process_adf_files():
    # Directory containing the .sdf files
    directory = os.path.join(dir, 'Atomic_densities')

    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        # Check if the file is a .sdf file
        if filename.endswith('.adf'):
            # Full path to the file
            file_path = os.path.join(directory, filename)
            
            # Check if the file exists and raise an error otherwise
            if not os.path.isfile(file_path):
                raise ValueError('File {} does not exist'.format(file_path))

            # Call the read_adf_file function on the file
            read_adf_file(file_path)


def write_dict_to_file(dictionary, filename):
    with open(filename, 'w') as f:
        json.dump(dictionary, f, indent=2)

def read_dict_from_file(filename):
    # Check if the file exists and raise an error otherwise
    if not os.path.isfile(filename):
        raise ValueError('File {} does not exist'.format(filename))
    with open(filename, 'r') as f:
        dictionary = json.load(f)
    return dictionary

#!!This updates the json file!!
#process_adf_files()
#write_dict_to_file(atom_dict, 'atom_dict.json')
######read_adf_file(os.path.join(dir,'Atomic_densities','atom_atom37.densities'))
atom_dict = read_dict_from_file('atom_dict.json')

lf = os.path.join(dir, 'NoSpherA2.log')
# Read the file
# 0 = k
# 1 = Thakakr FF
# 2 = Gaussian FF
# 3 = Thakkar Core Density
# 4 = Gaussian Core Density
# 5 = Thakkar Density
# 6 = ORCA Density
# 7 = ORCA Core Density
# 8 = r
df = pd.read_csv(lf, delim_whitespace=True, header=None)

upper = 13.0
nr = 2000001

x = np.linspace(0, upper, nr)
dx = x[1]
y = np.zeros_like(x)
z = np.zeros_like(x)

#coeffies = atom_dict['237']
coeffies = coefs
#for i in range(len(coeffies)):
#    coeffies[i] *= pow(pow(2, 7) * pow(exps[i], 3 ) / np.pi, 0.25)
#for i, v in enumerate(x):
#    y[i] = Rb_func(v, coeffies, exps, rad_exps)
#    z[i] = y[i] * (v)**2 * 4 * np.pi

#integral = sum(x**2 * y * x[1] * 4 * np.pi)
#print(integral)

integral2 = sum(df[8][1] * df[5] * (df[8])**2 * 4 * np.pi)
print("Thakkar Integral:",integral2)
integral2 = sum(df[8][1] * df[3] * (df[8])**2 * 4 * np.pi)
print("Thakkar Core Integral:",integral2)
integral2 = sum(df[8][1] * df[6] * (df[8])**2 * 4 * np.pi)
print("ORCA Integral:",integral2)
integral2 = sum(df[8][1] * df[7] * (df[8])**2 * 4 * np.pi)
print("ORCA Core Integral:",integral2)


integral2 = sum(df[8][1] * df[9] * (df[8])**2 * 4 * np.pi)
print("ORCA ECP Integral:",integral2)


# Create a figure with 2 subplots
fig, axs = plt.subplots(2,2, figsize=(12, 8))

#axs[0][0].set_title('ECP SF')
#axs[0][0].plot(df[0], df[1], label='Thakkar')
#axs[0][0].plot(df[0], df[2], label='Gaussian')


#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[0])], [exps[0]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[1])], [exps[1]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[2])], [exps[2]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[3])], [exps[3]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[4])], [exps[4]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[5])], [exps[5]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[6])], [exps[6]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[7])], [exps[7]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[8])], [exps[8]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[9])], [exps[9]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[10])], [exps[10]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[11])], [exps[11]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[12])], [exps[12]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[13])], [exps[13]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[14])], [exps[14]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[15])], [exps[15]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[16])], [exps[16]]))
#axs[0][0].plot(x,Rb_func(x, [abs(coeffies[17])], [exps[17]]))

axs[0][0].set_title('Radial dens diff')
axs[0][0].plot(df[8], (df[8])**2 * (df[9] + df[3] - df[5]) * 4 *np.pi, label='ECP +  Thakkar Core - Thakkar full')
axs[0][0].plot(df[8], (df[8])**2 * (df[9] + df[3] - df[6]) * 4 *np.pi, label='ECP +  Thakkar Core - ORCA full')
axs[0][0].plot(df[8], (df[8])**2 * (df[5] - df[6]) * 4 *np.pi, label='Thakkar - ORCA full')
axs[0][0].set_xlim(-0.1,3)
axs[0][0].legend()

axs[0][1].set_title('ECP Density')
#axs[0][1].plot(x, y, label='Gaussian python')
axs[0][1].plot(df[8], df[3], label='Thakkar')
#axs[0][1].plot(df[0]*dx_log, df[4], label='Gaussian')
axs[0][1].plot(df[8], df[5], label='Thakkar full')
axs[0][1].plot(df[8], df[6], label='ORCA')
axs[0][1].plot(df[8], df[9] + df[3], label='ECP + Thakkar Core')
axs[0][1].set_ylim(0,200)
axs[0][1].set_xlim(-0.1,2)
axs[0][1].legend()

axs[1][0].set_title('Valence Density')
#axs[1][0].plot(x, y, label='Gaussian python')
axs[1][0].plot(df[8], df[5] - df[3], label='Thakkar Valence')
#axs[1][0].plot(df[0]*dx_log, df[5] - df[4], label='Gaussian Valence')
axs[1][0].plot(df[8], df[6] - df[7], label='ORCA Valence')
axs[1][0].plot(df[8], df[9], label='ORCA ECP Valence')
axs[1][0].plot(df[8], df[5] - df[6], label='Thakkar vs ORCA')
axs[1][0].set_ylim(-63,63)
axs[1][0].set_xlim(-0.1,2)
axs[1][0].legend()

axs[1][1].set_title('ECP radial Density')
#axs[1][1].plot(x, z, label='Gaussian python')
axs[1][1].plot(df[8], (df[8])**2 * df[3] * 4 *np.pi, label='Thakkar')
#axs[1][1].plot(df[0]*dx_log, (df[0]*dx_log)**2 * df[4] * 4 *np.pi, label='Gaussian')
axs[1][1].plot(df[8], (df[8])**2 * df[5] * 4 *np.pi, label='Thakkar full')
axs[1][1].plot(df[8], (df[8])**2 * df[6] * 4 *np.pi, label='ORCA')
axs[1][1].plot(df[8], (df[8])**2 * df[9] * 4 *np.pi, label='ECP')
axs[1][1].plot(df[8], (df[8])**2 * (df[9] + df[3]) * 4 *np.pi, label='ECP +  Thakkar Core')
axs[1][1].set_xlim(-0.1,3)
axs[1][1].legend()

# Show the figure
#plt.show()
#exit(0)
fig.savefig("ECP_SF.png", dpi=300, bbox_inches='tight')
exit(0)

