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
           1,
           1,
           1,
           2,
           2,
           2,
           3
]

def Rb_func(r, Z, coefs = None, exponent = None, rad_exp = None):
    res = 0
    if exponent == None:
        raise ValueError("exponent is None")
    if coefs == None:
        raise ValueError("coefs is None")
    if rad_exp == None:
        raise ValueError("rad_exp is None")
    if Z < 0:
        raise ValueError("Z is negative")
    
    UL = -4*coefs[-1]*np.exp(-exponent[-1]*r**2) * \
        (1-3*exponent[-1]*r**2 + exponent[-1]**2 * r**4)
    res = UL
    
    for _i in range(3):
       Ul = 0#-UL
       fac = (2*_i+1)
       for k in range(3):
          Ul += 4*coefs[_i*3+k]*np.exp(-exponent[_i*3+k]*r**2) \
              * (1-3*exponent[_i*3+k]*r**2 + exponent[_i*3+k]**2 * r**4) 
       res += Ul * fac
    if (Z>0):
        res -= 4*Z/r**4#/(4*np.pi)
    return -res/(4*np.pi)**2
    for _i, c in enumerate(coefs):
        p1 = c
        p2 = np.exp(-exponent[_i] * r**2)
        if rad_exp[_i] == 0:
            p3 = -4 * (1-3*exponent[_i]*r**2 + exponent[_i]**2 * r**4)
        elif rad_exp[_i] == 1:
            p3 = -r * (9-16*exponent[_i]*r**2 + 4 * exponent[_i]**2 * r**4)
        elif rad_exp[_i] == 2:
            p3 = -4 * r**2 * (4-5*exponent[_i]*r**2 + exponent[_i]**2 * r**4)
        elif rad_exp[_i] == 3:
            p3 = -r**3 * (25-24*exponent[_i]*r**2+4*exponent[_i]**2*r**4)
        res += p1 * p2 * p3
    #for _i in range(4):
    #    c1 = coeffies[3*_i]
    #    ex1 = exponent[3*_i]
    #    e1 = np.exp(-ex1 * r**2)
    #    ep1 = np.exp(ex1 * r**2)
    #    if not _i == 3:
    #        c2 = coeffies[3*_i+1]
    #        ex2 = exponent[3*_i+1]
    #        e2 = np.exp(-ex2 * r**2)
    #        ep2 = np.exp(ex2 * r**2)
    #        c3 = coeffies[3*_i+2]
    #        ex3 = exponent[3*_i+2]
    #        e3 = np.exp(-ex3 * r**2)
    #        ep3 = np.exp(ex3 * r**2)
    #    else:
    #        c2 = 1
    #        ex2 = 1
    #        e2 = 1
    #        ep2 = 1
    #        c3 = 1
    #        ex3 = 1
    #        e3 = 1
    #        ep3 = 1
    #    #if _i==0:
    #    res += 4*e1*e2*e3*(c3*ep1*ep2*(1-3*ex3*r**2+ex3**2*r**4) +
    #                          ep3*(c1*ep2*(1-3*ex1*r**2+ex1**2*r**4) +
    #                               c2*ep1*(1-3*ex2*r**2+ex2**2*r**4)))
    #    #elif _i==1:
    #    #    res += r*e1*e2*e3*(c3*ep1*ep2*(9-16*ex3*r**2+4*ex3**2*r**4) +
    #    #                    ep3*(c1*ep2*(9-16*ex1*r**2+4*ex1**2*r**4) +
    #    #                         c2*ep1*(9-16*ex2*r**2+4*ex2**2*r**4)))
    #    #elif _i==2:
    #    #    res += 4*r**2*e1*e2*e3*(c3*ep1*ep2*(4-5*ex3*r**2+ex3**2*r**4) +
    #    #                    ep3*(c1*ep2*(4-5*ex1*r**2+ex1**2*r**4) +
    #    #                         c2*ep1*(4-5*ex2*r**2+ex2**2*r**4)))
    #    #elif _i==3:
    #    #    res += e1*e2*e3*r**3*(c3*ep1*ep2*(25-24*ex3*r**2+4*ex3**2*r**4) +
    #    #                    ep3*(c1*ep2*(25-24*ex1*r**2+4*ex1**2*r**4) +
    #    #                         c2*ep1*(25-24*ex2*r**2+4*ex2**2*r**4)))
    if (Z>0):
        res -= 4*Z/r**4#/(4*np.pi)
    return -res/(4*np.pi)**2

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

upper = 4.0
nr = 200000

x = np.linspace(0.00001, upper, nr)
dx = x[1]
y = np.zeros_like(x)
z = np.zeros_like(x)

for i, v in enumerate(x):
    y[i] = Rb_func(v, 28, coefs, exps, rad_exps)
    z[i] = y[i] * (v)**2 * 4 * np.pi

integral = sum(x**2 * y * x[1] * 4 * np.pi)
print("ECP Integral: ", integral)

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


# Create a figure with 3x2 subplots
fig, axs = plt.subplots(2,2, figsize=(12, 8))

#axs[0][0].set_title('ECP SF')
#axs[0][0].plot(df[0], df[1], label='Thakkar')
#axs[0][0].plot(df[0], df[2], label='Gaussian')


#axs[2][0].plot(x,Rb_func(x, 0, [coefs[0]],  [exps[0]] , [0]))
#axs[2][0].plot(x,Rb_func(x, 0, [coefs[1]],  [exps[1]] , [0]))
#axs[2][0].plot(x,Rb_func(x, 0, [coefs[2]],  [exps[2]] , [0]))
#axs[2][0].plot(x,Rb_func(x, 0, [coefs[3]],  [exps[3]] , [1]))
#axs[2][0].plot(x,Rb_func(x, 0, [coefs[4]],  [exps[4]] , [1]))
#axs[2][0].plot(x,Rb_func(x, 0, [coefs[5]],  [exps[5]] , [1]))
#axs[2][0].plot(x,Rb_func(x, 0, [coefs[6]],  [exps[6]] , [2]))
#axs[2][0].plot(x,Rb_func(x, 0, [coefs[7]],  [exps[7]] , [2]))
#axs[2][0].plot(x,Rb_func(x, 0, [coefs[8]],  [exps[8]] , [2]))
#axs[2][0].plot(x,Rb_func(x, 0, [coefs[9]],  [exps[9]] , [3]))
#axs[2][0].plot(x,Rb_func(x, 9, [0],  [exps[9]] , [3]))
#axs[2][0].set_ylim(-0.5,5)

axs[0][0].set_title('Radial dens diff')
axs[0][0].plot(df[8], (df[8])**2 * (df[9] + df[3] - df[5]) * 4 *np.pi, label='ECP +  Thakkar Core - Thakkar full')
axs[0][0].plot(df[8], (df[8])**2 * (df[9] + df[3] - df[6]) * 4 *np.pi, label='ECP +  Thakkar Core - ORCA full')
axs[0][0].plot(df[8], (df[8])**2 * (df[5] - df[6]) * 4 *np.pi, label='Thakkar - ORCA full')
axs[0][0].set_xlim(-0.1,3)
axs[0][0].legend()

axs[0][1].set_title('ECP Density')
axs[0][1].plot(x, y, label='Gaussian python')
axs[0][1].plot(df[8], df[3], label='Thakkar')
#axs[0][1].plot(df[0]*dx_log, df[4], label='Gaussian')
axs[0][1].plot(df[8], df[5], label='Thakkar full')
axs[0][1].plot(df[8], df[6], label='ORCA')
axs[0][1].plot(df[8], df[9], label='ECP')
axs[0][1].plot(df[8], df[9] + df[3], label='ECP + Thakkar Core')
axs[0][1].set_ylim(-5,200)
axs[0][1].set_xlim(-0.1,3)
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
axs[1][1].plot(x, z, label='Gaussian python')
axs[1][1].plot(df[8], (df[8])**2 * df[3] * 4 *np.pi, label='Thakkar')
#axs[1][1].plot(df[0]*dx_log, (df[0]*dx_log)**2 * df[4] * 4 *np.pi, label='Gaussian')
axs[1][1].plot(df[8], (df[8])**2 * df[5] * 4 *np.pi, label='Thakkar full')
axs[1][1].plot(df[8], (df[8])**2 * df[6] * 4 *np.pi, label='ORCA')
axs[1][1].plot(df[8], (df[8])**2 * df[9] * 4 *np.pi, label='ECP')
axs[1][1].plot(df[8], (df[8])**2 * (df[9] + df[3]) * 4 *np.pi, label='ECP +  Thakkar Core')
axs[1][1].set_xlim(-0.1,3)
axs[1][1].set_ylim(-3,160)
axs[1][1].legend()

# Show the figure
#plt.show()
#exit(0)
fig.savefig("ECP_SF.png", dpi=300, bbox_inches='tight')
exit(0)

