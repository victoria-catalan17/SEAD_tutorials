import numpy as np
from math import *

"""geomotric paramters"""
eta = 0.95      #given in slides, required for CLa

#fuselage
bf = 3.56

#main wing
S = 77.3
b = 26.21
A = 8.89
taper = 0.356   #Taper ratio main wing
sweep_quarter_c = np.radians(15)  #[rad]
MAC = 3.17

#horizontal tail plane
Sh = 15.61
bh = 11.09
Ah = 7.88
taperh = 0.41 #Taper ratio horizontal wing
sweep_quarter_ch = np.radians(20)  #[rad]
lh = 12.35
v_hv = 1        #T tail configuration

v0 = 404        #cruise speed,  Max cruise speed 430 [kts]
x_ac = 1        #location of aerodynamic center




"""def"""
"""stability"""
# C_La_h and C_La_w computation
M = v0*0.514444/295         #a=295 speed of sound 11000 feet
beta = sqrt(1-M**2)
sweep_half_ch = atan(tan(sweep_quarter_ch) - 1/Ah*(1-taperh)/(1+taperh))
sweep_half_c = atan(tan(sweep_quarter_c) - 1/Ah*(1-taper)/(1+taper))

CLa_h = (2*pi*Ah)/(2 + sqrt(4+((Ah*beta)/eta)**2 * (1+tan(sweep_half_ch)**2/beta**2)))
CLa_w = (2*pi*A)/(2 + sqrt(4+((A*beta)/eta)**2 * (1+tan(sweep_half_c)**2/beta**2)))

# C_La_A_h computation
cr =2*S/((1+taper)*b)
S_net = S - cr*bf
CLa_A_h = CLa_w*(1+2.15*bf/b)*(S_net/S + pi/2*bf**2/S)

# derivative of e over alpha
m_tv = 4.1/b*2      #fuselage height is 3.56m and tail height is 4.85m, 4.1 seems reasonable
r = 2*lh/b
K_esweep = 0.1124 + 0.1265*sweep_quarter_c + 0.1766*sweep_quarter_c**2/r**2 + 0.1024/r + 2
K_sweep0 = 0.1124/r**2 + 0.1024/r + 2
deda = K_esweep/K_sweep0*(r/(r**2 + m_tv**2) * 0.4876/sqrt(r**2 + 0.6319 + m_tv**2)//
        + (1 + (r**2/(r**2 + 0.7915 + 5.0734*m_tv**2))**0.3113) * (1 - sqrt(m_tv**2/(1 + m_tv**2)))) * CLa_w/(pi*A)

print(m_tv)
print(deda)

#x_cg = x_ac + CLa_h/CLa_A_h*(1-deda)*ShS*lh/MAC *v_hv**2 - SM

"""controllability"""









# if __name__ == '__main__':
#     print_hi('PyCharm')

