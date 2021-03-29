from utils.general import *
from scissor_plot import scissor_plot
import numpy as np
from math import *
""" Geometric parameters of the BAe Avro RJ85 """

original_design = False

if original_design == True:

    # Fuselage
    b_f = 3.56
    h_f = 3.56
    l_f = 26.5

    # Main wing
    S_w = 77.3
    b_w = 26.21
    A_w = 8.89
    tr_w = 0.356                            # Taper ratio main wing
    sweep_025_c_w = np.radians(15)          # [rad]
    MAC = 3.17

    # horizontal tail plane
    S_h = 15.61
    b_h = 11.09
    A_h = 7.88
    tr_h = 0.41                             # Taper ratio horizontal wing
    sweep_025_c_h = np.radians(20)          # [rad]
    l_h = 12.35
    vh_v = 1                                # For a T-tail configuration

    tail_height = 4.85

    # Flaps
    b_f_b_w = 0.78

    v_c = 404                               # cruise speed,  Max cruise speed 430 [kts]
    v_l = 122                               # Approach speed

    M_c = 0.68
    M_l = mach_number_landing(v_l * 0.5144444)

    # Coefficients
    CL_max_l = 3.62

else:

    # Fuselage
    b_f = 3.56
    h_f = 3.56
    l_f = 26.5

    # Main wing
    S_w = 77.3
    b_w = 28.1071                           # Changes when increasing the aspect ratio
    A_w = 8.89*1.15                         # Design modification
    tr_w = 0.2645                           # changes as the root chord stays the same
    sweep_025_c_w = np.radians(15)          # [rad]
    MAC = 3.06036                           # Changes as the span increases and the surface area is kept constant

    #horizontal tail plane
    S_h = 15.61
    b_h = 11.09
    A_h = 7.88
    tr_h = 0.41                             # Taper ratio horizontal wing
    sweep_025_c_h = np.radians(20)          # [rad]
    l_h = 12.35
    vh_v = 1                                # For a T-tail configuration

    tail_height = 4.85

    # Flaps
    b_f_b_w = 0.78

    v_c = 404                               # cruise speed,  Max cruise speed 430 [kts]
    v_l = 122                               # Approach speed

    M_c = 0.68
    M_l = mach_number_landing(v_l*0.514444)

    # Coefficients
    CL_max_l = 3.62


# SCISSOR PLOT

SM = 0.05                                   # Stability margin

scissor_plot(M_c, M_l, vh_v, b_f, h_f, l_f, S_w, A_w, b_w, tr_w, MAC, sweep_025_c_w, S_h, A_h, tr_h, l_h, sweep_025_c_h,
             tail_height, b_f_b_w, SM, original_design)
