# Imports
import numpy as np
import matplotlib.pyplot as plt

from general import *

SM = 0.05  # Safety margin of a 5%


def stability(M_c, b_f, h_f, S_w, A_w, b_w, tr_w, sw_025_c_w, A_h, tr_h, l_h, sw_025_c_h, ht_h, x_cg, SM, MAC, vh_v):

    sw_05_c_h = sweep_half_chord(sw_025_c_h, A_h, tr_h)
    sw_05_c_w = sweep_half_chord(sw_025_c_w, A_w, tr_w)

    CLa_h = alpha_lift_coefficient(M_c, A_h, sw_05_c_h)
    # C_La_A_h computation
    CLa_A_h = alpha_lift_coefficient_aircraft_minus_tail(M_c, A_w, sw_05_c_w, S_w, b_w, tr_w, b_f)

    x_ac = x_aerodynamic_center(b_f, h_f, S_w, MAC, S_w, b_w, tr_w, sw_025_c_w, CLa_A_h)
    # # derivative of e over alpha
    deda = wing_downwash_gradient(ht_h, l_h, b_w, A_w, tr_w, sw_025_c_w, M_c)

    #x_mac_le = x_MAC_leading_edge_loc(sw_025_c_w, tr_w, A_w, b_w)

    Sh_S = (x_cg - x_ac + SM) / (CLa_h/CLa_A_h*(1-deda)*l_h/MAC*vh_v**2)

    return Sh_S


def controllability(M_l, Vh_V, b_f, h_f, l_f, S_w, A_w, b_w, tr_w, MAC, sw_025_c_w, S_h, A_h, tr_h, l_h, sw_025_c_h, ht_h,
                    b_f_b_w, x_cg):
    # Geometry
    sw_05_c_w = sweep_half_chord(sw_025_c_w, A_w, tr_w)
    # V_h = tail_volume_coefficient(S_h, S_w, l_h, MAC)

    CL_h = -0.35 * A_h ** (1 / 3)
    CL_A_h = 3.43  # Maximum lift coefficient of the aircraft for landing

    CLa_A_h = alpha_lift_coefficient_aircraft_minus_tail(M_l, A_w, sw_05_c_w, S_w, b_w, tr_w, b_f)

    x_ac = x_aerodynamic_center(b_f, h_f, S_w, MAC, S_w, b_w, tr_w, sw_025_c_w, CLa_A_h)

    Cm_ac = aircraft_aerodynamic_pitching_moment(M_l, A_w, sw_025_c_w, S_w, b_w, tr_w, MAC, b_f, h_f, l_f, b_f_b_w, x_ac)

    Sh_S = (x_cg + (Cm_ac / CL_A_h) - x_ac) / ((CL_h * l_h * Vh_V*Vh_V) / (CL_A_h * MAC))

    return Sh_S


def scissor_plot(M_c, M_l, Vh_V, b_f, h_f, l_f, S_w, A_w, b_w, tr_w, MAC, sw_025_c_w, S_h, A_h, tr_h, l_h, sw_025_c_h, ht_h,
                 b_f_b_w, SM, original_design):

    x_cg = np.linspace(0, 20, 1000)
    x_cg = (x_cg - x_MAC_leading_edge_loc(sw_025_c_w, tr_w, A_w, b_w)) / MAC

    Sh_S_stab = stability(M_c, b_f, h_f, S_w, A_w, b_w, tr_w, sw_025_c_w, A_h, tr_h, l_h, sw_025_c_h, ht_h, x_cg, SM, MAC, Vh_V)
    Sh_S_stab_SM0 = stability(M_c, b_f, h_f, S_w, A_w, b_w, tr_w, sw_025_c_w, A_h, tr_h, l_h, sw_025_c_h, ht_h, x_cg, 0, MAC, Vh_V)

    Sh_S_control = controllability(M_l, Vh_V, b_f, h_f, l_f, S_w, A_w, b_w, tr_w, MAC, sw_025_c_w, S_h, A_h, tr_h, l_h,
                                   sw_025_c_h, ht_h, b_f_b_w, x_cg)

    plt.figure(figsize=(7, 5))
    plt.plot(x_cg, Sh_S_stab, label= "Stability curve" )
    plt.plot(x_cg, Sh_S_stab_SM0, label= "Stability curve (SM = 0)")
    plt.plot(x_cg, Sh_S_control, label= "Controllability curve")
    plt.axvline(0.556)
    plt.axvline(0.3616)

    if original_design:
        plt.plot(np.linspace(0,1,10), np.ones(10)*0.202, label= r"$\frac{S_h}{S}$ for the BAe Avro RJ85")

    # x-axis
    plt.xlabel(r'$\frac{x_{cg}}{MAC}$ [-]', fontsize=15)
    plt.xticks(fontsize=12)
    plt.xlim([0, 1])

    # y-axis
    plt.ylabel(r'$\frac{S_h}{S}$ [-]', fontsize=15)
    plt.yticks(fontsize=12)
    plt.ylim([0, 1])

    plt.grid(color='lightgray', linestyle='dashdot')
    plt.legend()

    plt.tight_layout()
    plt.show()
    plt.close()

