import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import scipy.special as sp
dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir)

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

rad_exps = [2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2
]

l_facs = [
    0.28209479**2,
    0.48860251**2,
    0.63078313**2,
    0.74635267**2
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
    
    UL = l_facs[3] * 4*coefs[-1]*exponent[-1]*np.exp(-exponent[-1]*r**2) * \
        (exponent[-1]*r**2 - 1) + Z/r**3
    res = UL
    
    for l in range(3):
        Ul = 0#-UL
        for k in range(3):
            Ul += l_facs[l] * 4*coefs[l*3+k]*exponent[l*3+k]*np.exp(-exponent[l*3+k]*r**2) \
               * (exponent[l*3+k]*r**2 - 1)
        #Ul += Z/r**3
        res += Ul
    
    return res/(4*np.pi)**2

def Ul_func(r, Z, coefs = None, exponent = None, rad_exp = None):
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
    
    for l in range(3):
        Ul = -UL
        for k in range(3):
            Ul += 4*coefs[l*3+k]*exponent[l*3+k]*np.exp(-exponent[l*3+k]*r**2) \
               * (exponent[l*3+k]*r**2 - 1)
        Ul += Z/r
        for m in range(-l,l+1,1):
            fac = np.sqrt(2*l+1)*sp.factorial(l-abs(m))/sp.factorial(l+abs(m))
            res += (Ul) * fac
    return res/(4*np.pi)

def get_numbers_from_filename(filename):
    # Use a regular expression to find all groups of digits in the filename
    numbers = re.findall(r'\d+', filename)
    str = ''
    str = str.join(numbers)

    # Convert the numbers to integers and return them
    return str

def Slm(l,m,theta,phi):
    if m == 0:
        return np.real(sp.sph_harm(m,l,phi,theta))
    elif m > 0:
        if l == 1:
            return np.real(sp.sph_harm(-m,l,phi,theta) - sp.sph_harm(m,l,phi,theta)) / np.sqrt(2)
        else:
            if m%2 == 1:
                return np.real(sp.sph_harm(-m,l,phi,theta) - sp.sph_harm(m,l,phi,theta)) / np.sqrt(2)
            else:
                return np.real(sp.sph_harm(-m,l,phi,theta) + sp.sph_harm(m,l,phi,theta)) / np.sqrt(2)
    elif m < 0:
        if l == 1:
            return np.imag(sp.sph_harm(-m,l,phi,theta) + sp.sph_harm(m,l,phi,theta)) / np.sqrt(2)
        else:
            if m%2 == 1:
                return np.imag(sp.sph_harm(-m,l,phi,theta) + sp.sph_harm(m,l,phi,theta)) / np.sqrt(2)
            else:
                return np.imag(sp.sph_harm(-m,l,phi,theta) - sp.sph_harm(m,l,phi,theta)) / np.sqrt(2)
    
def Slm_self(l,m,theta,phi):
    if l==0:
        if m != 0:
            return 0.0
        else:
            return np.sqrt(1/4/np.pi)*sp.lpmv(m,l,np.cos(theta))
    elif l==1:
        if m == 0:
            return np.sqrt(3/4/np.pi)*sp.lpmv(m,l,np.cos(theta))
        elif m == 1:
            return -np.sqrt(3/4/np.pi)*np.cos(phi)*sp.lpmv(m,l,np.cos(theta))
        elif m == -1:
            return -np.sqrt(3/4/np.pi)*np.sin(phi)*sp.lpmv(m,l,np.cos(theta)) *2
        else:
            return 0.0
    elif l==2:
        if m == 0:
            return np.sqrt(5/4/np.pi)*sp.lpmv(m,l,np.cos(theta))
        elif m == 1:
            return -np.sqrt(5/12/np.pi)*np.cos(phi)*sp.lpmv(m,l,np.cos(theta))
        elif m == -1:
            return -np.sqrt(5/12/np.pi)*np.sin(phi)*sp.lpmv(m,l,np.cos(theta))
        elif m == 2:
            return np.sqrt(5/48/np.pi)*np.cos(2*phi)*sp.lpmv(m,l,np.cos(theta))
        elif m == -2:
            return np.sqrt(5/48/np.pi)*np.sin(2*phi)*sp.lpmv(m,l,np.cos(theta))
        else:
            return 0.0

def check_SLM():
    for l in range(4):
        for m in range(-l,l+1,1):
            print("l: ", l, "m: ", m)
            #Check if Slm and Slm_self are the same
            for theta in np.linspace(0, np.pi, 1):
                for phi in np.linspace(0, 2*np.pi, 1):
                        print("theta: ", theta, "phi: ", phi)
                        print(f"Scipy:     {Slm(l,m,theta,phi):14.8f}")
                        #print(f"selfbuild: {Slm_self(l,m,theta,phi):14.8f}")
                        
#check_SLM()

lf = os.path.join(dir, 'core_dens.dat')
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

upper = 3.0
nr = 500000

x = np.linspace(0.000001, upper, nr)
dx = x[1]
y = np.zeros_like(x)
z = np.zeros_like(x)

def br_test():
    y = np.zeros_like(x)
    z = np.zeros_like(x)
    a = np.zeros_like(x)
    y_lanl = np.zeros_like(x)
    z_lanl = np.zeros_like(x)
    a_lanl = np.zeros_like(x)
    c_lanl = np.zeros_like(x)
    lanl_density = np.zeros_like(x)
    KG_dens = np.zeros_like(x)
    e_br = [1.42  ,      3.56    ,    6.35    ,    10.4    ,    20.7, 32.4, 76.5, 142.1, 341.871, 1000.205, 220.0000]
    d_br0 = [4.202568   ,    13.194570  ,    10.795897  ,    26.648345  ,    33.130248  ,    45.359127  ,    35.346799  ,    -31.766999 ,    -633.120024,    -47.285358 , 3.0]
    d_br1 = [3.611982   ,    10.720320  ,    5.828031   ,    34.390992  ,    61.018209  ,    115.349011 ,    257.006175 ,    161.631632 ,    -144.676033,    -9.413723  , 5.0]
    d_br2 = [3.484330   ,    -0.211214  ,    9.199904   ,    -3.542208  ,    43.713385  ,    77.076843  ,    245.662800 ,    509.898585 ,    -330.399360,    -21.536241 , 7.0]
    n_br = [2,2,2,2,2,2,2,2,2,1,0]
    lanl_e1 = [54.1980682,32.9053558,13.6744890,3.0341152]
    lanl_d1 = [3.0000000,27.3430642,118.8028847,43.4354876]
    lanl_n1 = [0,1,2,2]
    lanl_e2 = [54.2563340,26.0095593,28.2012995,9.4341061,2.5321764]
    lanl_d2 = [5.0000000,25.0504252,92.6157463,95.8249016,26.2684983]
    lanl_n2 = [0,1,2,2,2]
    lanl_e3 = [87.6328721,61.7373377,32.4385104,8.7537199,1.6633189]
    lanl_d3 = [3.0000000,22.5533557,178.1241988,76.9924162,9.4818270]
    lanl_n3 = [0,1,2,2,2]
    lanl_e4 = [213.6143969,41.0585380,8.7086530,2.6074661]
    lanl_d4 = [-28.0000000,-134.9268852,-41.9271913,-5.9336420]
    lanl_n4 = [1,2,2,2]

    for i in range(11):
        y += 0.28209479**2 *d_br0[i] * np.exp(-e_br[i]*x**2) * x**(n_br[i]-2)
        z += 0.48860251**2 *d_br1[i] * np.exp(-e_br[i]*x**2) * x**(n_br[i]-2)
        a += 0.63078313**2 *d_br2[i] * np.exp(-e_br[i]*x**2) * x**(n_br[i]-2)
        if n_br[i] == 2:
            KG_dens += 0.28209479**2 * 4 * d_br0[i] *e_br[i] * np.exp(-e_br[i]*x**2) * (-1 + e_br[i]*x**2) + \
                       0.48860251**2 * 4 * d_br1[i] *e_br[i] * np.exp(-e_br[i]*x**2) * (-1 + e_br[i]*x**2) + \
                       0.63078313**2 * 4 * d_br2[i] *e_br[i] * np.exp(-e_br[i]*x**2) * (-1 + e_br[i]*x**2)
        elif n_br[i] == 1:
            KG_dens += 0.28209479**2 *d_br0[i] * np.exp(-e_br[i]*x**2) * (1 + 4*e_br[i]**2 * x**4) / x**3 + \
                       0.48860251**2 *d_br1[i] * np.exp(-e_br[i]*x**2) * (1 + 4*e_br[i]**2 * x**4) / x**3 + \
                       0.63078313**2 *d_br2[i] * np.exp(-e_br[i]*x**2) * (1 + 4*e_br[i]**2 * x**4) / x**3
        elif n_br[i] == 0:#  used to be 2
            KG_dens += 0.28209479**2 *4*d_br2[i] * np.exp(-e_br[i]*x**2) * (1 + e_br[i]*x**2 + e_br[i]**2 * x**4) / x**4 + \
                       0.48860251**2 *4*d_br1[i] * np.exp(-e_br[i]*x**2) * (1 + e_br[i]*x**2 + e_br[i]**2 * x**4) / x**4 + \
                       0.63078313**2 *4*d_br0[i] * np.exp(-e_br[i]*x**2) * (1 + e_br[i]*x**2 + e_br[i]**2 * x**4) / x**4
    
    for i in range(4):
        y_lanl += 0.28209479**2 *lanl_d1[i] * np.exp(-lanl_e1[i]*x**2) * x**(lanl_n1[i])
        c_lanl += 0.74635267**2 *lanl_d4[i] * np.exp(-lanl_e4[i]*x**2) * x**(lanl_n4[i])
        
        
        if lanl_n1[i] == 0:
            lanl_density -= 0.28209479**2 *4*lanl_d1[i] *lanl_e1[i] * np.exp(-lanl_e1[i]*x**2) * (-1 + lanl_e1[i]*x**2)
        elif lanl_n1[i] == 1:
            lanl_density -= 0.28209479**2 *lanl_d1[i] * np.exp(-lanl_e1[i]*x**2) * (1 - 8*lanl_e1[i]*x**2 + 4*lanl_e1[i]**2 * x**4) / x
        elif lanl_n1[i] == 2:
            lanl_density -= 0.28209479**2 *4*lanl_d1[i] * np.exp(-lanl_e1[i]*x**2) * (1 - 3*lanl_e1[i]*x**2 + lanl_e1[i]**2 * x**4)
        
        if lanl_n4[i] == 0:
            lanl_density -= 0.74635267**2 *4*lanl_d4[i] *lanl_e4[i] * np.exp(-lanl_e4[i]*x**2) * (-1 + lanl_e4[i]*x**2)
        elif lanl_n4[i] == 1:
            lanl_density -= 0.74635267**2 *lanl_d4[i] * np.exp(-lanl_e4[i]*x**2) * (1 - 8*lanl_e4[i]*x**2 + 4*lanl_e4[i]**2 * x**4) / x
        elif lanl_n4[i] == 2:
            lanl_density -= 0.74635267**2 *4*lanl_d4[i] * np.exp(-lanl_e4[i]*x**2) * (1 - 3*lanl_e4[i]*x**2 + lanl_e4[i]**2 * x**4)
    for i in range(5):
        z_lanl += 0.48860251**2 *lanl_d2[i] * np.exp(-lanl_e2[i]*x**2) * x**(lanl_n2[i])
        a_lanl += 0.63078313**2 *lanl_d3[i] * np.exp(-lanl_e3[i]*x**2) * x**(lanl_n3[i])
        
        
        if lanl_n2[i] == 0:
            lanl_density -= 0.48860251**2 *4 * lanl_d2[i] *lanl_e2[i] * np.exp(-lanl_e2[i]*x**2) * (-1 + lanl_e2[i]*x**2)
        elif lanl_n2[i] == 1:
            lanl_density -= 0.48860251**2 *lanl_d2[i] * np.exp(-lanl_e2[i]*x**2) * (1 - 8*lanl_e2[i]*x**2 + 4*lanl_e2[i]**2 * x**4) / x
        elif lanl_n2[i] == 2:
            lanl_density -= 0.48860251**2 *4*lanl_d2[i] * np.exp(-lanl_e2[i]*x**2) * (1 - 3*lanl_e2[i]*x**2 + lanl_e2[i]**2 * x**4)
        
        if lanl_n3[i] == 0:
            lanl_density -= 0.63078313**2 *4*lanl_d3[i] *lanl_e3[i] * np.exp(-lanl_e3[i]*x**2) * (-1 + lanl_e3[i]*x**2)
        elif lanl_n3[i] == 1:
            lanl_density -= 0.63078313**2 *lanl_d3[i] * np.exp(-lanl_e3[i]*x**2) * (1 - 8*lanl_e3[i]*x**2 + 4*lanl_e3[i]**2 * x**4) / x
        elif lanl_n3[i] == 2:
            lanl_density -= 0.63078313**2 *4*lanl_d3[i] * np.exp(-lanl_e3[i]*x**2) * (1 - 3*lanl_e3[i]*x**2 + lanl_e3[i]**2 * x**4)

    b = a+z+y + 28/x
    y += 28/x
    z += 28/x
    a += 28/x
    KG_dens += 28/x**3
    
    KG_radial_dens = KG_dens * x**2 * 4 * np.pi
    lanl_radial_dens = lanl_density * x**2 * 4 * np.pi
    
    fig, axs = plt.subplots(4,1, figsize=(14, 8))
    
    axs[0].plot(x, y_lanl, linestyle='-', label='l=0')
    axs[0].plot(x, z_lanl, linestyle='-', label='l=1')
    axs[0].plot(x, a_lanl, linestyle='-', label='l=2')
    axs[0].plot(x, y_lanl + z_lanl + a_lanl - c_lanl, linestyle='-', label='sum')
    axs[0].plot(x, -c_lanl, linestyle='-', label='(-) l=3')
    
    axs[0].set_xlabel('r')
    axs[0].set_ylabel(r'$U_l(r)$ LANL')
    axs[0].legend()
    #axs[0].set_ylim(0,15.0)
    axs[0].set_xlim(0,2.5)

    axs[1].plot(x, y, linestyle='-', label='l=0')
    axs[1].plot(x, z, linestyle='-', label='l=1')
    axs[1].plot(x, a, linestyle='-', label='l=2')
    axs[1].plot(x, b, linestyle='-', label='sum')
    
    axs[1].set_xlabel('r')
    axs[1].set_ylabel(r'$(U_l(r) - \frac{N_C}{r})$ KG')
    axs[1].legend()
    axs[1].set_ylim(0.0,150.0)
    axs[1].set_xlim(0,2.5)

    axs[2].plot(x, KG_dens, linestyle='-', label='KG')
    axs[2].plot(x, lanl_density, linestyle='-', label='LANL')
    
    axs[2].set_ylim(-50,300.0)
    axs[2].set_xlim(0,2.5)
    axs[2].legend()
    
    axs[3].plot(x, KG_radial_dens, linestyle='-', label='KG')
    axs[3].plot(x, lanl_radial_dens, linestyle='-', label='LANL')
    
    axs[3].set_ylim(-20,100.0)
    axs[3].set_xlim(0,2.5)
    axs[3].legend()
    
    
    
    fig.savefig("Br_test.png", dpi=300, bbox_inches='tight')
    exit(0)
    
#br_test()

for i, v in enumerate(x):
    y[i] = Rb_func(v, 28, coefs, exps, rad_exps)
    z[i] = y[i] * (v)**2 * 4 * np.pi

integral = sum(x**2 * y * (x[1]-x[0]) * 4 * np.pi)
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
axs[1][1].set_ylim(-3,170)
axs[1][1].legend()

# Show the figure
#plt.show()
#exit(0)
fig.savefig("ECP_SF.png", dpi=300, bbox_inches='tight')
exit(0)

