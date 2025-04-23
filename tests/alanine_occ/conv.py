import numpy as np
from scipy.linalg import lstsq
import math

d_sph = [[-2.23631590e-07, -3.53832159e-07, -2.58237987e-08,  2.31980290e-08,  2.78416882e-08],
         [-1.01831319e-04, -1.60205144e-04, -1.07348696e-05,  1.12237407e-05,  1.29138533e-05],
         [-1.50466912e-11, -3.06924162e-12,  3.40847592e-12,  1.89421761e-11, -6.22555026e-13],
         [-4.72251792e-12, -2.63507740e-13, -8.20445649e-12, -1.53928068e-12, -2.16501938e-11],
         [-2.55296136e-12, -2.28554956e-11,  3.70994203e-12, -1.06395400e-11,  1.99902497e-12],
         [ 3.37295847e-02, -2.79691189e-01, -3.61839495e-01,  5.21558638e-01,  7.19497413e-01],
         [-1.04104745e-01,  4.26810177e-01, -3.10575972e-01, -6.75468303e-01,  5.04256052e-01],
         [ 1.91122376e-10, -1.09088261e-10,  6.88387031e-11, -3.59681349e-10, -4.45798841e-11],
         [ 3.30307194e-03,  5.29720339e-03,  2.91701031e-04, -6.36271854e-04, -4.81584577e-04]
         ]
d_cart = [
[
  [0.31436209E-05, 0.21860350E-05,-0.53296559E-05, 0.11492704E-05,-0.14605753E-04,-0.10659744E-05],
  [0.11740160E-05, 0.81639618E-06,-0.19904122E-05, 0.42920627E-06,-0.54546614E-05,-0.39809858E-06],
  [0.96939490E-07, 0.67410520E-07,-0.16435001E-06, 0.35439924E-07,-0.45039599E-06,-0.32871335E-07]
],
[
  [0.14450885E-02, 0.98178638E-03,-0.24268749E-02,0.53306785E-03,-0.66130697E-02,-0.44312211E-03],
  [0.53968245E-03, 0.36665774E-03,-0.90634019E-03,0.19907941E-03,-0.24697156E-02,-0.16548829E-03],
  [0.44562034E-04, 0.30275239E-04,-0.74837273E-04,0.16438155E-04,-0.20392650E-03,-0.13664508E-04]
],
[
  [0.57025342E-09,-0.21165612E-09,-0.35859731E-09,-0.25698299E-10,-0.12669449E-09,0.14069766E-09],
  [0.21296673E-09,-0.79045050E-10,-0.13392168E-09,-0.95972816E-11,-0.47315297E-10,0.52544919E-10],
  [0.17584842E-10,-0.65268164E-11,-0.11058026E-10,-0.79245563E-12,-0.39068640E-11,0.43386782E-11]
],
[
  [0.24504441E-10,0.88044038E-10,-0.11254848E-09,-0.89369315E-09,-0.10877273E-10,-0.33866979E-09],
  [0.91514233E-11,0.32880908E-10,-0.42032332E-10,-0.33375846E-09,-0.40622240E-11,-0.12647955E-09],  
  [0.75564074E-12,0.27150043E-11,-0.34706451E-11,-0.27558717E-10,-0.33542126E-12,-0.10443523E-10]
],
[
 [-0.18917203E-09,0.25001498E-09,-0.60842949E-10,0.82517272E-10,-0.94344652E-09,0.15314181E-09],
 [-0.70648149E-10,0.93370544E-10,-0.22722395E-10,0.30816883E-10,-0.35233934E-09,0.57192308E-10],
 [-0.58334772E-11,0.77096845E-11,-0.18762073E-11,0.25445760E-11,-0.29092956E-10,0.47224170E-11]
],
[
  [0.10362720E+02, -0.11166574E+02, 0.80385369E+00, 0.29699961E+02, -0.11545305E+02, -0.14936286E+02],
  [0.38700593E+01, -0.41702663E+01, 0.30020703E+00, 0.11091741E+02, -0.43117074E+01, -0.55781022E+01],
  [0.31955406E+00, -0.34434241E+00, 0.24788348E-01, 0.91585443E+00, -0.35602131E+00, -0.46058860E+00]
],
[
 [-0.12700719E+02, 0.15181775E+02, -0.24810558E+01, 0.20815065E+02, 0.17618195E+02, -0.12820191E+02],
 [-0.47432079E+01, 0.56697825E+01, -0.92657460E+00, 0.77735897E+01, 0.65796874E+01, -0.47878259E+01],
 [-0.39165068E+00, 0.46815872E+00, -0.76508046E-01, 0.64187185E+00, 0.54329034E+00, -0.39533482E+00]
],
[
 [-0.97010432E-08, 0.51461568E-08, 0.45548864E-08, -0.18402023E-08, -0.45030281E-08, 0.28415763E-08],
 [-0.36229496E-08, 0.19218827E-08, 0.17010669E-08, -0.68724158E-09, -0.16817000E-08, 0.10612145E-08],
 [-0.29915000E-09, 0.15869147E-09, 0.14045853E-09, -0.56746117E-10, -0.13885938E-09, 0.87625376E-10]
],
[
 [-0.52492163E-01, -0.26227649E-01, 0.78719812E-01, -0.19879214E-01, 0.21866199E+00, 0.12041057E-01],
 [-0.19603712E-01, -0.97949725E-02, 0.29398685E-01, -0.74240872E-02, 0.81661461E-01, 0.44968506E-02],
 [-0.16186950E-02, -0.80877914E-03, 0.24274742E-02, -0.61301314E-03, 0.67428556E-02, 0.37130874E-03]
]
]

d_exp = [1.09473708e+01,  3.33929702e+00,  1.28840460e+00]
d_con = [2.19767951e-01,  6.55547363e-01,  2.86573259e-01]

all_d_sph_norm = []
all_d_cart = []

for set_index in range(len(d_sph)):
    d_s = d_sph[set_index]
    d_c = d_cart[set_index]

    for i,e in enumerate(d_exp):
        n_fac = pow(2048*pow(e,7)/pow(np.pi,3),0.25)
        row = [n_fac * d_con[i] * sph for sph in d_s]
        all_d_sph_norm.append(row)
        all_d_cart.append(d_c[i])
    
# Convert to numpy arrays
all_d_sph_norm = np.array(all_d_sph_norm)
all_d_cart = np.array(all_d_cart)

print(f"Combined data shape: {all_d_sph_norm.shape} -> {all_d_cart.shape}")

# Solve for the transformation matrix
# For each Cartesian component i, find coefficients that map spherical components to it
transformation_matrix = np.zeros((5, 6))

for i in range(6):  # For each Cartesian component (xx, yy, zz, xy, xz, yz)
    b = all_d_cart[:, i]
    transformation_matrix[:, i], _, _, _ = lstsq(all_d_sph_norm, b)

np.set_printoptions(formatter={'float': lambda x: f"{x:.9f}" if abs(x) < 1e10 else f"{x:.4e}"})
print("derived transformation matrix:")
print(transformation_matrix.T)

# Test the transformation
# d_cart_calc = np.zeros((3, 6))

# for i in range(3):
#     # Apply the transformation matrix to each normalized spherical component
#     d_cart_calc[i] = np.dot(all_d_sph_norm[i], transformation_matrix)
    
#     # Compare with the given Cartesian components
#     print(f"Exponent {d_exp[i]}:")
#     print(f"  Original: {d_cart[i]}")
#     print(f"  Calculated: {d_cart_calc[i]}")
#     print(f"  difference: {np.array(d_cart[i]) - d_cart_calc[i]}")
#     print()

def find_closest_sqrt_fractions(target, count=20, max_denominator=2048):
    """
    Find the closest square root fractions to a target value.
    
    Args:
        target: The target floating point value
        count: Number of closest approximations to return
        max_denominator: Maximum allowed denominator in the fractions
    
    Returns:
        list of tuples: [(a, b, sqrt_val, error), ...] representing √(a/b)
    """
    target_squared = target ** 2
    approximations = []
    
    # Try different denominators and numerators
    for b in range(1, max_denominator + 1):
        for a in range(1, max_denominator * 2):  # Try a range of numerators
            sqrt_val = math.sqrt(a/b)
            error = abs(sqrt_val - target)
            
            # Add to our list of approximations
            approximations.append((a, b, sqrt_val, error))
            
            # If we've collected enough, sort and keep only the best ones
            if len(approximations) > count * 10:
                approximations.sort(key=lambda x: x[3])  # Sort by error
                approximations = approximations[:count * 2]  # Keep top results
    
    # Final sort and trim
    approximations.sort(key=lambda x: x[3])
    
    # Process the results to simplify fractions and extract square factors
    simplified_results = []
    for a, b, sqrt_val, error in approximations[:count*5]:
        # Find reduced form
        gcd = math.gcd(a, b)
        a //= gcd
        b //= gcd
        
        # Find perfect squares
        a_sqrt_factor = 1
        b_sqrt_factor = 1
        
        # Check for perfect square factors in numerator
        for i in range(2, int(math.sqrt(a)) + 1):
            while a % (i**2) == 0:
                a //= i**2
                a_sqrt_factor *= i
        
        # Check for perfect square factors in denominator
        for i in range(2, int(math.sqrt(b)) + 1):
            while b % (i**2) == 0:
                b //= i**2
                b_sqrt_factor *= i
        
        # Calculate final values
        final_a = a * (a_sqrt_factor**2)
        final_b = b * (b_sqrt_factor**2)
        sqrt_val = math.sqrt(final_a/final_b)
        error = abs(sqrt_val - target)
        
        # Add to results if not already there
        result = (final_a, final_b, a_sqrt_factor, b_sqrt_factor, a, b, sqrt_val, error)
        if not any(r[0] == final_a and r[1] == final_b for r in simplified_results):
            simplified_results.append(result)
    
    # Final sort and trim
    simplified_results.sort(key=lambda x: x[7])
    return simplified_results[:count]

known = [
    [-0.5/np.sqrt(3.0),	0,	    0,	    0.5,	0],
    [-0.5/np.sqrt(3.0),	0,	    0,	    -0.5,	0],
    [np.sqrt(1.0/3.0),  0,	    0,	    0,	    0],
    [0,	                0,	    0,	    0,	    1.0],
    [0,	                1.0,	0,	    0,	    0],
    [0,	                0,	    1.0,	0,	    0]
]

known = np.array(known)
ratio = np.zeros((6, 5))
for i in range(6):
    for j in range(5):
        if known[i][j] == 0:
            ratio[i][j] = 0.0
        else:
            ratio[i][j] = transformation_matrix[j][i] / known[i][j]
        if abs(ratio[i][j]) < 1E-4:
            ratio[i][j] = 0.0
non_zero_values = ratio[ratio != 0]
target = np.average(non_zero_values)
results = find_closest_sqrt_fractions(target)

print(f"Target value: {target}")
print("Top 10 approximations:")
for i, (final_a, final_b, a_sqrt, b_sqrt, a, b, approx, error) in enumerate(results):
    if a_sqrt > 1 or b_sqrt > 1:
        # Complex representation with extracted square roots
        print(f"{i+1}. √({final_a}/{final_b}) = {a_sqrt}√({a})/{b_sqrt}√({b}) = {approx}, error: {error:.10e}")
    else:
        # Simple representation
        print(f"{i+1}. √({final_a}/{final_b}) = {approx}, error: {error:.10e}")
exit(0)