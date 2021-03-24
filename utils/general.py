# Imports
import numpy as np


# FLIGHT CONDITIONS

def mach_number_cruise(airspeed):

    return airspeed/295         #a=295 speed of sound 11000 feet


# GEOMETRY

def sweep_half_chord(sweep_quarter_chord, aspect_ratio, taper_ratio):

    sw05ch = np.arctan(np.tan(sweep_quarter_chord) - 1/aspect_ratio * (1-taper_ratio)/(1+taper_ratio))

    return sw05ch


def root_chord(surface_area, wingspan, taper_ratio):

    cr = 2 * surface_area / ((1 + taper_ratio) * wingspan)

    return cr


def net_wing_surface_area(surface_area_wing, wingspan_wing, taper_ratio_wing, fuselage_width):

    cr = root_chord(surface_area_wing, wingspan_wing, taper_ratio_wing)
    net_surface = surface_area_wing - cr * fuselage_width

    return net_surface


def m_tv_ratio(tail_height, wingspan):          # Computed assuming that the root chord of the wing
                                                # is parallel to the ground (a0 = 0)
    return (tail_height * 2) / wingspan


def r_ratio(tail_arm, wingspan):

    return (2 * tail_arm) / wingspan


def tail_volume_coefficient(surface_area_tail, surface_area_wing, tail_arm, mean_aerodynamic_chord):

    return (surface_area_tail * tail_arm) / (surface_area_wing * mean_aerodynamic_chord)


#COEFFICIENTS

def lift_coefficient():
    pass


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

    return deda
