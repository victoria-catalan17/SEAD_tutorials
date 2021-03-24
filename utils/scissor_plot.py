# Imports
import numpy as np
import matplotlib.pyplot as plt

from general import *

SM = 0.05   # Safety margin of a 5%


def stability(M, Vh_V, bf, S_w, A_w, b_w, tr_w, sw_025_c_w, A_h, tr_h, l_h, sw_025_c_h, ht_h, x_cg, SM):

    sweep_05_ch = sweep_half_chord(sw_025_c_h, A_h, tr_h)
    sweep_05_c = sweep_half_chord(sw_025_c_w, A_w, tr_w)

    CLa_h = alpha_lift_coefficient(M, A_h, sweep_05_ch)
    # C_La_A_h computation
    CLa_A_h = alpha_lift_coefficient_aircraft_minus_tail(M, A_w, sweep_05_c, S_w, b_w, tr_w, bf)

    # # derivative of e over alpha
    deda = wing_downwash_gradient(ht_h, l_h, b_w, A_w, tr_w, sw_025_c_w, M)

    Sh_S = np.ones(x_cg.shape) - SM

    return Sh_S


def controllability(M, Vh_V, bf, S_w, A_w, b_w, tr_w, sw_025_c_w, S_h, A_h, tr_h, l_h, sw_025_c_h, ht_h, MAC, x_cg):

    CL_h = -0.35 * A_h**(1/3)
    CL_A_h = 1.35               # Will change depending on the values of the plot.

    Vh = tail_volume_coefficient(S_h, S_w, l_h, MAC)
    print(Vh)

    Sh_S = np.ones(x_cg.shape)*2

    return Sh_S


def scissor_plot(M, Vh_V, bf, S_w, A_w, b_w, tr_w, sw_025_c_w, S_h, A_h, tr_h, l_h, sw_025_c_h, ht_h, SM, MAC):

    x_cg = np.linspace(0, 20, 1000)

    Sh_S_stab = stability(M, Vh_V, bf, S_w, A_w, b_w, tr_w, sw_025_c_w, A_h, tr_h, l_h, sw_025_c_h, ht_h, x_cg, SM)
    Sh_S_stab_SM0 = stability(M, Vh_V, bf, S_w, A_w, b_w, tr_w, sw_025_c_w, A_h, tr_h, l_h, sw_025_c_h, ht_h, x_cg, 0)

    Sh_S_control = controllability(M, Vh_V, bf, S_w, A_w, b_w, tr_w, sw_025_c_w, S_h, A_h, tr_h, l_h, sw_025_c_h, ht_h, MAC, x_cg)

    plt.figure(figsize=(5, 5))
    plt.plot(x_cg/MAC, Sh_S_stab)
    plt.plot(x_cg/MAC, Sh_S_stab_SM0)
    plt.plot(x_cg/MAC, Sh_S_control)

    plt.grid(color='lightgray', linestyle='dashdot')

    plt.tight_layout()
    plt.show()
