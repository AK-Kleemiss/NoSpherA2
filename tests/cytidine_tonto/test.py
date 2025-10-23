import math
a = -1.9395659680943534e-05

c1 = -0.00538935
coef_1 = -0.53431809926E-02
c_1 = 2266.1767785000002
norm_c1 = 45614.363993120329
target_1 = 0.24469292E-04

c2 = -0.04023472
coef_2 = -0.39890039230E-01
c_2 = 340.87010191000002
norm_c2 = 4272.8783706497243
target_2 = 0.44122217E-04

c3 = -0.18008184
coef_3 = -0.17853911985
c_3 = 77.363135170000007
norm_c3 = 669.34799877167518
target_3 = 0.64935938E-04

c4 = -0.46828858
coef_4 = -0.46427684959
c_4 = 21.479644940000000
norm_c4 = 134.90196175243631
target_4 = 0.64587509E-04

c5 = -0.44692617
coef_5 = -0.44309745172
c_5  = 6.6589433099999997
norm_c5 = 31.206210921346205
target_5 = 0.25609696E-04


b = 3.9566960145875465e-05

d1 = 1.0000000000000000
d_1 = 0.80975976000000005
norm_d1 = 2.2409369190886355
target_n2 = 0.24071875E-04

def double_factorial(n:int) -> int:
    if n <= 0:
        return 1
    else:
        return n * double_factorial(n - 2)
l = 0
#pow(pow(2, 7 + 4 * l) * pow(exp, 3 + 2 * l) / math.pi / pow(double_factorial(2 * l + 1), 2), 0.25);

#pow(pow(2 ,3 + 4 * l) * pow(exp, 2 * l + 3) / math.pi**3 / pow(double_factorial(l),2),0.25)

difference = pow(pow(2,4) * math.pi * math.pi / (double_factorial(l)**2 / double_factorial(2*l +1)**2),0.25)
print("faktor             ", f"{difference:<.8e}", f"{1/difference:<.8e}")
#fac = pow(4.0*exp, 0.5*.l+0.75) * (1.0/(.norm*sqrt(.l.double_factorial)))



def compare(n: int) -> None:
    # Select the appropriate variables based on n
    c_dict = {1: c1, 2: c2, 3: c3, 4: c4, 5: c5}
    coef_dict = {1: coef_1, 2: coef_2, 3: coef_3, 4: coef_4, 5: coef_5}
    norm_dict = {1: norm_c1, 2: norm_c2, 3: norm_c3, 4: norm_c4, 5: norm_c5}
    target_dict = {1: target_1, 2: target_2, 3: target_3, 4: target_4, 5: target_5}
    exp_dict = {1: c_1, 2: c_2, 3: c_3, 4: c_4, 5: c_5}
    
    c = c_dict[n]
    coef = coef_dict[n]
    norm_c = norm_dict[n]
    target = target_dict[n]
    exp = exp_dict[n]
    
    calc = pow(pow(2 , 3 + 4 * l) * pow(exp, 2 * l + 3) / math.pi**3 / pow(double_factorial(l),2),0.25)

    print(f"---- Target: {n} ----", f"{target:<.8e} ------------")
    print(f"c{n} from NSA2       {c:<.8e}  c{n} from basis {coef:<.8e} ratio {c/coef:<.8e}")

    print(f"norm from Tonto     {calc:<.8e}")
    print(f"norm from tonto fac {calc/difference:<.8e}")
    print(f"norm from tonto fac {calc*difference:<.8e}")
    print(f"norm NSA2           {norm_c:<.8e}")
 
    print(f"no norm             {a*c:<.8e} ratio: {(a*c)/target:<.8f}")
    print(f"no_norm / faktor    {a*c/difference:<.8e} ratio: {(a*c/difference)/target:<.8f}")
    print(f"no_norm * faktor    {a*c*difference:<.8e} ratio: {(a*c*difference)/target:<.8f}")
    print(f"no norm_2           {a*coef:<.8e} ratio: {(a*coef)/target:<.8f}")
    print(f"no norm_2 / faktor  {a*coef/difference:<.8e} ratio: {(a*coef/difference)/target:<.8f}")
    print(f"no norm_2 * faktor  {a*coef*difference:<.8e} ratio: {(a*coef*difference)/target:<.8f}")
    print(f"incl norm           {a*c*norm_c:<.8e} ratio: {(a*c*norm_c)/target:<.8f}")
    print(f"incl norm / factor  {a*c*norm_c/difference:<.8e} ratio: {(a*c*norm_c/difference)/target:<.8f}")
    print(f"incl norm * factor  {a*c*norm_c*difference:<.8e} ratio: {(a*c*norm_c*difference)/target:<.8f}")
    print(f"incl norm           {a*c*calc:<.8e} ratio: {(a*c*calc)/target:<.8f}")
    print(f"incl norm / factor  {a*c*calc/difference:<.8e} ratio: {(a*c*calc/difference)/target:<.8f}")
    print(f"incl norm * factor  {a*c*calc*difference:<.8e} ratio: {(a*c*calc*difference)/target:<.8f}")



for i in range(1,5):
    compare(i)