from utils.general import *
from scissor_plot import scissor_plot
import numpy as np
from math import *
""" Geometric parameters of the BAe Avro RJ85 """

#fuselage
b_f = 3.56
h_f = 3.56
l_f = 26.5

#main wing
S_w = 77.3
b_w = 26.21
A_w = 8.89
tr_w = 0.356   #Taper ratio main wing
sweep_025_c_w = np.radians(15)  #[rad]
MAC = 3.17

#horizontal tail plane
S_h = 15.61
b_h = 11.09
A_h = 7.88
tr_h = 0.41 #Taper ratio horizontal wing
sweep_025_c_h = np.radians(20)  #[rad]
l_h = 12.35
vh_v = 1        #T tail configuration

tail_height = 4.85

# Flaps
b_f_b_w = 0.78

v_c = 404        #cruise speed,  Max cruise speed 430 [kts]
v_l = 125

M_c = 0.68
M_l = mach_number_landing(v_l*0.5144444)

#print(sqrt(1-M_c**2)*A_w,sqrt(1-M_l**2)*A_w)
x_ac = 1        #location of aerodynamic center

# Coefficients
CL_max_l = 3.62

SM = 0.05   # Stability margin


scissor_plot(M_c, M_l, vh_v, b_f, h_f, l_f, S_w, A_w, b_w, tr_w, MAC, sweep_025_c_w, S_h, A_h, tr_h, l_h, sweep_025_c_h,
             tail_height, b_f_b_w, SM)

#x_cg = x_ac + CLa_h/CLa_A_h*(1-deda)*ShS*lh/MAC *v_hv**2 - SM
