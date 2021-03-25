# Imports
import numpy as np
import matplotlib.pyplot as plt

from general import *

SM = 0.05   # Safety margin of a 5%


def stability(M_c, b_f, S_w, A_w, b_w, tr_w, sw_025_c_w, A_h, tr_h, l_h, sw_025_c_h, ht_h, x_cg, SM):

    sw_05_c_h = sweep_half_chord(sw_025_c_h, A_h, tr_h)
    sw_05_c_w = sweep_half_chord(sw_025_c_w, A_w, tr_w)

    CLa_h = alpha_lift_coefficient(M_c, A_h, sw_05_c_h)
    # C_La_A_h computation
    CLa_A_h = alpha_lift_coefficient_aircraft_minus_tail(M_c, A_w, sw_05_c_w, S_w, b_w, tr_w, b_f)

    # # derivative of e over alpha
    deda = wing_downwash_gradient(ht_h, l_h, b_w, A_w, tr_w, sw_025_c_w, M_c)

    Sh_S = np.ones(x_cg.shape) - SM

    return Sh_S


def controllability(M_l, Vh_V, b_f, h_f, S_w, A_w, b_w, tr_w, MAC, sw_025_c_w, S_h, A_h, tr_h, l_h, sw_025_c_h, ht_h, x_cg):

    # Geometry
    sw_05_c_w = sweep_half_chord(sw_025_c_w, A_w, tr_w)
    V_h = tail_volume_coefficient(S_h, S_w, l_h, MAC)

    CL_h = -0.35 * A_h**(1/3)
    CL_A_h = 1.35               # Will change depending on the values of the plot.

    CLa_A_h = alpha_lift_coefficient_aircraft_minus_tail(M_l, A_w, sw_05_c_w, S_w, b_w, tr_w, b_f)

    Sh_S = np.ones(x_cg.shape)*2

    return Sh_S


def scissor_plot(M_c, M_l, Vh_V, b_f, h_f, S_w, A_w, b_w, tr_w, MAC, sw_025_c_w, S_h, A_h, tr_h, l_h, sw_025_c_h, ht_h, SM):

    x_cg = np.linspace(0, 20, 1000)

    Sh_S_stab = stability(M_c, b_f, S_w, A_w, b_w, tr_w, sw_025_c_w, A_h, tr_h, l_h, sw_025_c_h, ht_h, x_cg, SM)
    Sh_S_stab_SM0 = stability(M_c, b_f, S_w, A_w, b_w, tr_w, sw_025_c_w, A_h, tr_h, l_h, sw_025_c_h, ht_h, x_cg, 0)

    Sh_S_control = controllability(M_l, Vh_V, b_f, h_f, S_w, A_w, b_w, tr_w, MAC, sw_025_c_w, S_h, A_h, tr_h, l_h, sw_025_c_h, ht_h, x_cg)

    plt.figure(figsize=(5, 5))
    plt.plot(x_cg/MAC, Sh_S_stab)
    plt.plot(x_cg/MAC, Sh_S_stab_SM0)
    plt.plot(x_cg/MAC, Sh_S_control)

    plt.grid(color='lightgray', linestyle='dashdot')

    # x-axis
    plt.xlabel(r'$\frac{x_{cg}}{MAC}$ [-]', fontsize=15)
    plt.xticks(fontsize=12)

    # y-axis
    plt.ylabel(r'$\frac{S_h}{S}$ [-]', fontsize=15)
    plt.yticks(fontsize=12)

    plt.tight_layout()
    plt.show()
    plt.close()