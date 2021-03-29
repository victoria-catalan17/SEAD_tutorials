# Imports
import numpy as np


# FLIGHT CONDITIONS

def mach_number_cruise(airspeed):

    return airspeed/295         #a=295 speed of sound 11000 feet


def mach_number_landing(airspeed):

    return airspeed/344         #a=340 speed of sound at sea level


# GEOMETRY

def sweep_half_chord(sweep_quarter_chord, aspect_ratio, taper_ratio):

    sw05ch = np.arctan(np.tan(sweep_quarter_chord) - 1/aspect_ratio * (1-taper_ratio)/(1+taper_ratio))

    return sw05ch


def root_chord(surface_area, wingspan, taper_ratio):

    cr = 2 * surface_area / ((1 + taper_ratio) * wingspan)

    return cr


def any_chord(surface_area, wingspan, taper_ratio, position_wrt_center):

    cr = root_chord(surface_area, wingspan, taper_ratio)
    ct = cr * taper_ratio

    return cr - (2 * (cr - ct)) / (wingspan) * position_wrt_center


def net_wing_surface_area(surface_area_wing, wingspan_wing, taper_ratio_wing, fuselage_width):

    cr = root_chord(surface_area_wing, wingspan_wing, taper_ratio_wing)
    net_surface = surface_area_wing - cr * fuselage_width

    return net_surface


def surface_area_part_wing(surface_area, wingspan, taper_ratio, position_from_root_chord, fuselage_width):

    cr = root_chord(surface_area, wingspan, taper_ratio)
    c_other = any_chord(surface_area, wingspan, taper_ratio, position_from_root_chord+fuselage_width/2)

    return (cr + c_other) / 2 * position_from_root_chord


def m_tv_ratio(tail_height, wingspan):          # Computed assuming that the root chord of the wing
                                                # is parallel to the ground (a0 = 0)
    return (tail_height * 2) / wingspan


def r_ratio(tail_arm, wingspan):

    return (2 * tail_arm) / wingspan


def tail_volume_coefficient(surface_area_tail, surface_area_wing, tail_arm, mean_aerodynamic_chord):

    return (surface_area_tail * tail_arm) / (surface_area_wing * mean_aerodynamic_chord)


def x_MAC_leading_edge_loc(surface_area, sweep_quarter_chord, taper_ratio, aspect_ratio, wingspan_wing, MAC):

    sweep_LE = np.arctan(np.tan(sweep_quarter_chord) + 4/aspect_ratio*0.25*(1-taper_ratio)/(1+taper_ratio))
    cr = root_chord(surface_area, wingspan_wing, taper_ratio)
    x_t = (cr - MAC) * wingspan_wing / (cr * (1 - taper_ratio) * 2)

    return np.tan(sweep_LE)*x_t + 10 - 0.7491052114    # Got this with x_t = b_f/2


#COEFFICIENTS

def lift_coefficient(W, V_e, S_w):

    rho_0 = 1.225

    L = (2* W) / (rho_0 * V_e*V_e * S_w)

    return L


def alpha_lift_coefficient(mach_number, aspect_ratio, sweep_half_chord):

    beta = np.sqrt(1-mach_number*mach_number)
    eta = 0.95         # As given in the slides

    cla = (2 * np.pi * aspect_ratio) / \
          (2 + np.sqrt(4 + (aspect_ratio * beta / eta)**2 * (1 + np.tan(sweep_half_chord)**2 / (beta*beta))))

    return cla


def alpha_lift_coefficient_aircraft_minus_tail(mach_number, aspect_ratio_wing, sweep_half_chord_wing,
                                               surface_area_wing, wingspan_wing, taper_ratio_wing, fuselage_width):

    cla_w = alpha_lift_coefficient(mach_number, aspect_ratio_wing, sweep_half_chord_wing)
    surface_area_wing_net = net_wing_surface_area(surface_area_wing, wingspan_wing, taper_ratio_wing, fuselage_width)

    cla_a_h = cla_w * (1 + 2.15 * fuselage_width / wingspan_wing) * surface_area_wing_net/surface_area_wing + \
              np.pi/2 * fuselage_width * fuselage_width / surface_area_wing

    return cla_a_h


def wing_downwash_gradient(tail_height, tail_arm, wingspan, aspect_ratio_wing,
                           taper_ratio_wing, sweep_quarter_chord, mach_number):

    sweep_half_chord_wing = sweep_half_chord(sweep_quarter_chord, aspect_ratio_wing, taper_ratio_wing)

    m_tv = m_tv_ratio(tail_height, wingspan)
    r = r_ratio(tail_arm, wingspan)
    cla_w = alpha_lift_coefficient(mach_number, aspect_ratio_wing, sweep_half_chord_wing)

    k_e_sweep = (0.1124 + 0.1265 * sweep_quarter_chord + 0.1766 * sweep_quarter_chord*sweep_quarter_chord) / (r*r) + \
                 0.1024 / r + 2
    k_e_sweep0 = 0.1124 / (r*r) + 0.1024 / r + 2

    deda = (k_e_sweep / k_e_sweep0) * (r / (r*r + m_tv*m_tv) * 0.4876 / np.sqrt(r*r + 0.6319 + m_tv*m_tv) +
                                       (1 + ((r*r) / (r*r + 0.7915 + 5.0734 * m_tv*m_tv))**0.3113) *
                                       (1 - np.sqrt((m_tv*m_tv) / (1 + m_tv*m_tv)))) * cla_w / (np.pi*aspect_ratio_wing)
    #print(k_e_sweep,k_e_sweep0,r,m_tv)
    return deda


# AERODYNAMIC CENTER

def x_ac_fuselage_contribution_1(fuselage_width, fuselage_height, S_w, mean_aerodynamic_chord, cla_a_h):

    l_fn = 10   #nose till leading edge root chord position
    x_ac_fuselage_1 = - 1.8/cla_a_h * fuselage_width * fuselage_height * l_fn/(S_w*mean_aerodynamic_chord)
    return x_ac_fuselage_1


def x_ac_fuselage_contribution_2(surface_area_wing, wingspan_wing, taper_ratio_wing, fuselage_width,
                                 mean_aerodynamic_chord, sweep_quarter_chord):
    c_g = surface_area_wing/wingspan_wing
    x_ac_fuselage_2 = 0.273/(1 + taper_ratio_wing)* fuselage_width*c_g*(wingspan_wing-fuselage_width)\
                      /(mean_aerodynamic_chord**2*(wingspan_wing + 2.15*fuselage_width))*np.tan(sweep_quarter_chord)
    return x_ac_fuselage_2


def x_ac_nacelles_contribution(wingspan_wing, taper_ratio_wing, surface_area_wing, mean_aerodynamic_chord,
                               sweep_quarter_chord, cla_a_h, original_design):

    if original_design:

        k_n = -4
        b_n = 1.40      #Nacelle width
        MAC_y = wingspan_wing/2*(1 + 2*taper_ratio_wing)/(3 + 3*taper_ratio_wing)       #MAC y location
        engine_1_percentage_Y = 0.315       # engine y location
        engine_2_percentage_Y = 0.5         # engine y location
        l_n1 = 0.25*mean_aerodynamic_chord + (MAC_y-wingspan_wing/2*engine_1_percentage_Y)*np.tan(sweep_quarter_chord) + 2.6
        l_n2 = 0.25*mean_aerodynamic_chord + (MAC_y-wingspan_wing/2*engine_2_percentage_Y)*np.tan(sweep_quarter_chord) + 2.6

        x_ac_nacelles = 2*k_n*b_n**2*l_n1/(surface_area_wing*mean_aerodynamic_chord*cla_a_h)\
                       + 2*k_n*b_n**2*l_n2/(surface_area_wing*mean_aerodynamic_chord*cla_a_h)

    else:

        k_n = -4
        b_n = 1.40*1.15      #Nacelle width
        MAC_y = wingspan_wing/2*(1 + 2*taper_ratio_wing)/(3 + 3*taper_ratio_wing)       #MAC y location
        engine_1_percentage_Y = 0.315       # engine y location
        l_n1 = 0.25*mean_aerodynamic_chord + (MAC_y-wingspan_wing/2*engine_1_percentage_Y)*np.tan(sweep_quarter_chord) + 2.6

        x_ac_nacelles = 2*k_n*b_n**2*l_n1/(surface_area_wing*mean_aerodynamic_chord*cla_a_h)
    return x_ac_nacelles


def x_aerodynamic_center(fuselage_width, fuselage_height, S_w, mean_aerodynamic_chord, surface_area_wing,
                         wingspan_wing, taper_ratio_wing, sweep_quarter_chord, cla_a_h, original_design):
    if original_design:
        x_ac_wing = 0.26        #/MAC
    else:
        x_ac_wing = 0.29        #/MAC
    x_ac = x_ac_wing + x_ac_fuselage_contribution_1(fuselage_width, fuselage_height, S_w, mean_aerodynamic_chord, cla_a_h)\
           + x_ac_fuselage_contribution_2(surface_area_wing, wingspan_wing, taper_ratio_wing, fuselage_width,
                                          mean_aerodynamic_chord, sweep_quarter_chord)\
           + x_ac_nacelles_contribution(wingspan_wing, taper_ratio_wing, surface_area_wing, mean_aerodynamic_chord,
                                        sweep_quarter_chord, cla_a_h, original_design)

    return x_ac


# Moment coefficients
def aircraft_aerodynamic_pitching_moment_wing(aspect_ratio_wing, sweep_quarter_chord_wing):

    Cm_0 = -0.3365992129

    Cm_ac_w = Cm_0 * (aspect_ratio_wing * np.cos(sweep_quarter_chord_wing) * np.cos(sweep_quarter_chord_wing)) / \
                     (aspect_ratio_wing + 2 * np.cos(sweep_quarter_chord_wing))

    print("-------------------Cm_ac WING-------------------")
    print("Cm_0:         ", Cm_0)
    print("Cm_ac_wing:   ", Cm_ac_w)

    return Cm_ac_w


def aircraft_aerodynamic_pitching_moment_flaps(A_w, sw_025_c_w, S_w, b_w, tr_w, MAC, b_f, b_f_b_w, x_ac, original_design):

    if original_design:
        # Triple-slotted Fowler flaps, landing setting: d_f = 33 [deg]
        # dc_cf = 0.7
        # cf = 0.9543235602       # Computed from the area and span of the flaps
        dc = 0.6680265
        c_prime = MAC + dc
        # cf_c_prime = 0.24864955
        c_prime_c = c_prime / MAC

        # Coefficients (found by using the data above)
        mu_1 = 0.215
        mu_2 = 0.99
        mu_3 = 0.033

    else:
        # Triple-slotted Fowler flaps, landing setting: d_f = 33 [deg]
        # dc_cf = 0.7
        # cf = Sf / (b_f_b_w*b_w) = 0.8899111      # Computed from the area and span of the flaps
        dc = 0.6229378
        c_prime = MAC + dc
        # cf_c_prime = 0.241607
        c_prime_c = c_prime / MAC

        # Coefficients (found by using the data above)
        mu_1 = 0.215
        mu_2 = 1.058
        mu_3 = 0.028

    dcl_max = 1.9 * c_prime_c
    CL_l = 3.43

    S_wf = 2*surface_area_part_wing(S_w, b_w, tr_w, b_f_b_w*b_w/2, b_f)

    Cm_025_c = mu_2 * (-mu_1 * dcl_max * c_prime_c - (CL_l + dcl_max * (1-S_wf/S_w)) * c_prime_c/8 * (c_prime_c-1)) + \
               0.7 * (A_w / (1 + 2/A_w)) * mu_3 * dcl_max * np.tan(sw_025_c_w)

    Cm_ac_f = Cm_025_c - CL_l * (0.25 - x_ac)

    print("-------------------Cm_ac FLAPS-------------------")
    print("x_ac:         ", x_ac)
    print("Swf:          ", S_wf)
    print("Cm_1/4:       ", Cm_025_c)
    print("Cm_ac_flaps:  ", Cm_ac_f)

    return Cm_ac_f


def aircraft_aerodynamic_pitching_moment_fuselage(M, A_w, sw_025_c_w, S_w, b_w, tr_w, MAC, b_f, h_f, l_f):

    sw_05_c_w = sweep_half_chord(sw_025_c_w, A_w, tr_w)

    CL_0 = 0.6          # Obtained from comparison with other aircraft
    CLa_A_h = alpha_lift_coefficient_aircraft_minus_tail(M, A_w, sw_05_c_w, S_w, b_w, tr_w, b_f)

    Cm_ac_fus = -1.8 * (1 - 2.5*b_f / l_f) * ((np.pi * b_f * h_f * l_f * CL_0) / (4 * S_w * MAC * CLa_A_h))

    print("-------------------Cm_ac FLAPS-------------------")
    print("CL_A-h:        ", CLa_A_h)
    print("Cm_ac_fuselage:", Cm_ac_fus)

    return Cm_ac_fus


def aircraft_aerodynamic_pitching_moment(M, A_w, sw_025_c_w, S_w, b_w, tr_w, MAC, b_f, h_f, l_f, b_f_b_w, x_ac, original_design):

    print("----------------PITCHING MOMENT COEFFICIENT----------------")

    Cm_ac_w = aircraft_aerodynamic_pitching_moment_wing(A_w, sw_025_c_w)
    Df_Cm_ac = aircraft_aerodynamic_pitching_moment_flaps(A_w, sw_025_c_w, S_w, b_w, tr_w, MAC, b_f, b_f_b_w, x_ac, original_design)
    Dfus_Cm_ac = aircraft_aerodynamic_pitching_moment_fuselage(M, A_w, sw_025_c_w, S_w, b_w, tr_w, MAC, b_f, h_f, l_f)

    Cm_ac = Cm_ac_w + Df_Cm_ac + Dfus_Cm_ac

    print("---------Cm_ac: ", Cm_ac)

    return Cm_ac
