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



#for i in range(1,5):
#    compare(i)
    
    

d4 = -1.15134806E-07  
d5 = 1.10898632E-06 
d6 = -2.44742291E-06

tonto_d4 = -0.19941933E-06
tonto_d5 = 0.19208206E-05
tonto_d6 = -0.42390608E-05

print(" D functions ----   "
      f"{d4:<.8e} {d5:<.8e} {d6:<.8e}")
print("Tonto D functions   "
        f"{tonto_d4:<.8e} {tonto_d5:<.8e} {tonto_d6:<.8e}")
print("Ratio               "
      f"{d4/tonto_d4:<.8e} {d5/tonto_d5:<.8e} {d6/tonto_d6:<.8e}")

print((d4/tonto_d4)**2)
print(math.sqrt(0.3333))

print(math.sqrt(1.5))


f1 = 3.72056418E-04 
f2 = -9.80555572E-14 
f3 = -2.35484325E-13
f4 = -7.39529468E-14 
f5 = -1.89159647E-13 
f6 = -2.20064316E-03 
f7 = -2.41091929E-13 
f8 = -2.20064316E-03
f9 = -9.54164514E-14 
f10=  3.61091530E-15

tonto_f1 = 0.16638869E-03 
tonto_f2 =-0.43851778E-13 
tonto_f3 =-0.10531179E-12
tonto_f4 = -0.33072763E-13 
tonto_f5 = -0.84594766E-13 
tonto_f6 = -0.98415754E-03 
tonto_f7 = -0.10781959E-12 
tonto_f8 = -0.98415754E-03
tonto_f9 = -0.42671534E-13  
tonto_f10 = 0.16148504E-14

print(" F functions ----   "
      f"{f1:<.8e} {f2:<.8e} {f3:<.8e} {f4:<.8e} {f5:<.8e} "
      f"{f6:<.8e} {f7:<.8e} {f8:<.8e} {f9:<.8e} {f10:<.8e}")
print("Tonto F functions   "
      f"{tonto_f1:<.8e} {tonto_f2:<.8e} {tonto_f3:<.8e} {tonto_f4:<.8e} {tonto_f5:<.8e} "
      f"{tonto_f6:<.8e} {tonto_f7:<.8e} {tonto_f8:<.8e} {tonto_f9:<.8e} {tonto_f10:<.8e}")
print("Ratio               "
      f"{f1/tonto_f1:<.8e} {f2/tonto_f2:<.8e} {f3/tonto_f3:<.8e} {f4/tonto_f4:<.8e} {f5/tonto_f5:<.8e} "
      f"{f6/tonto_f6:<.8e} {f7/tonto_f7:<.8e} {f8/tonto_f8:<.8e} {f9/tonto_f9:<.8e} {f10/tonto_f10:<.8e}")

print((f1/tonto_f1)**2)
print((f2/tonto_f2)**2)
print((f3/tonto_f3)**2)
print((f4/tonto_f4)**2)
print((f5/tonto_f5)**2)
print((f6/tonto_f6)**2)
print((f7/tonto_f7)**2)
print((f8/tonto_f8)**2)
print((f9/tonto_f9)**2)
print((f10/tonto_f10)**2)

g11 = -2.21752490E-02
g12 = -5.14765074E-02 
g13 = -5.14765074E-02 
g14 = -2.72421014E-13  
g15 = 4.98557839E-14 
g16 = -2.27440097E-13
g17 =  5.09601293E-14 
g18 =  -1.58563739E-14  
g19 =  7.79416965E-14 
g20 =  -7.14942370E-02 
g21 =  -7.14942370E-02
g22 = -1.02953019E-01  
g23 = 2.26403626E-14 
g24 = -1.33348554E-14 
g25 = -2.49425372E-13

tonto_g11=  -0.61209502E-02
tonto_g12= -0.14208866E-01 
tonto_g13= -0.14208866E-01 
tonto_g14= -0.75195344E-13  
tonto_g15= 0.13761504E-13 
tonto_g16= -0.62779431E-13
tonto_g17=  0.14066332E-13 
tonto_g18=  -0.43767750E-14  
tonto_g19=  0.21513952E-13 
tonto_g20=  -0.19734284E-01 
tonto_g21=  -0.19734284E-01
tonto_g22= -0.28417733E-01 
tonto_g23= 0.62493338E-14 
tonto_g24= -0.36807698E-14 
tonto_g25= -0.68847944E-13

print(" G functions ----   "
      f"{g11:<.8e} {g12:<.8e} {g13:<.8e} {g14:<.8e} {g15:<.8e} "
      f"{g16:<.8e} {g17:<.8e} {g18:<.8e} {g19:<.8e} {g20:<.8e} "
      f"{g21:<.8e} {g22:<.8e} {g23:<.8e} {g24:<.8e} {g25:<.8e} ")
print("Tonto G functions   "
        f"{tonto_g11:<.8e} {tonto_g12:<.8e} {tonto_g13:<.8e} {tonto_g14:<.8e} {tonto_g15:<.8e} "
        f"{tonto_g16:<.8e} {tonto_g17:<.8e} {tonto_g18:<.8e} {tonto_g19:<.8e} {tonto_g20:<.8e} "
        f"{tonto_g21:<.8e} {tonto_g22:<.8e} {tonto_g23:<.8e} {tonto_g24:<.8e} {tonto_g25:<.8e} ")
print("Ratio               "
      f" {g11/tonto_g11:<.8e}  {g12/tonto_g12:<.8e}  {g13/tonto_g13:<.8e}  {g14/tonto_g14:<.8e}  {g15/tonto_g15:<.8e}"
      f" {g16/tonto_g16:<.8e}  {g17/tonto_g17:<.8e}  {g18/tonto_g18:<.8e}  {g19/tonto_g19:<.8e}  {g20/tonto_g20:<.8e}"
      f" {g21/tonto_g21:<.8e}  {g22/tonto_g22:<.8e}  {g23/tonto_g23:<.8e}  {g24/tonto_g24:<.8e}  {g25/tonto_g25:<.8e}")


print((g11/tonto_g11)**2)
print((g12/tonto_g12)**2)
print((g13/tonto_g13)**2)
print((g14/tonto_g14)**2)
print((g15/tonto_g15)**2)
print((g16/tonto_g16)**2)
print((g17/tonto_g17)**2)
print((g18/tonto_g18)**2)
print((g19/tonto_g19)**2)
print((g20/tonto_g20)**2)
print((g21/tonto_g21)**2)
print((g22/tonto_g22)**2)
print((g23/tonto_g23)**2)
print((g24/tonto_g24)**2)
print((g25/tonto_g25)**2)

print(math.sqrt(13.125))