import numpy as np
from math import *

"""geomotric paramters"""
eta = 0.95      #given in slides, required for CLa
S = 77.3
b = 26.21
A = 8.89
MAC = 3.17
taper = 0.356   #Taper ratio main wing
sweep_quarter_c = np.radians(15)  #[rad]

#horizontal tail plane
Sh = 15.61
bh = 11.09
Ah = 7.88
taperh = 0.41 #Taper ratio horizontal wing
sweep_quarter_ch = np.radians(20)  #[rad]
lh = 12.35

v_hv = 1        #T tail configuration
x_ac = 1        #location of aerodynamic center

#def

# C_La_h computation
beta = sqrt(1-M**2)
CLa_h = (2*pi*A_h)/(2 + sqrt(4+((A_h*beta)/eta)**2 * (1+tan(sweep_half_ch)**2/beta**2)))
CLa_w = (2*pi*A)/(2 + sqrt(4+((A*beta)/eta)**2 * (1+tan(sweep_half_c)**2/beta**2)))

# C_La_A-h computation
S_net = S-cr*bf
CLa_A-h = CLa_w*(1+2.15*bf/b)*(S_net/S + pi/2*bf**2/S)

# derivative of e over alpha
m_tv = 1
r = 2*lh/b
K_esweep = 0.1124 + 0.1265*sweep_quarter_c + 0.1766*sweep_quarter_c**2/r**2 + 0.1024/r + 2
K_sweep0 = 0.1124/r**2 + 0.1024/r + 2
deda = K_esweep/K_sweep0*(r/(r**2 + m_tv**2) * 0.4876/sqrt(r**2 + 0.6319 + m_tv**2)//
        + (1 + (r**2/(r**2 + 0.7915 + 5.0734*m_tv**2)**0.3113)) * (1 - sqrt(m_tv**2/(1 + m_tv**2)))) * CLa_w/pi*A













# if __name__ == '__main__':
#     print_hi('PyCharm')

