
def x_ac_fuselage_contribution_1(fuselage_width, fuselage_height, S_w, mean_aerodynamic_chord):

    cla_a_h = alpha_lift_coefficient_aircraft_minus_tail(mach_number, aspect_ratio_wing, sweep_half_chord_wing,
                                                         surface_area_wing, wingspan_wing, taper_ratio_wing,
                                                         fuselage_width)
    l_fn = 10   #nose till leading edge root chord postion

    x_ac_fuselage_1 = - 1.8/cla_a_h * fuselage_width * fuselage_height * l_fn/(S_w*mean_aerodynamic_chord)

    return x_ac_fuselage_1


def x_ac_fuselage_contribution_2(surface_area_wing, wingspan_wing, taper_ratio_wing, fuselage_width, wingspan,
                                 mean_aerodynamic_chord, sweep_quarter_chord):

    c_g = surface_area_wing/wingspan_wing
    x_ac_fuselage_2 = 0.273/(1 + taper_ratio_wing)* fuselage_width*c_g*(wingspan-fuselage_width) / \
                      (mean_aerodynamic_chord**2*(wingspan + 2.15*fuselage_width))*tan(sweep_quarter_chord)

    return x_ac_fuselage_2


def x_ac_nacelles_contribution(wingspan_wing, taper_ratio_wing, surface_area_wing, mean_aerodynamic_chord,
                               sweep_quarter_chord):

    k_n = -4
    b_n = 1.40      #Nacelle width
    MAC_y = wingspan_wing/2*(1 + 2*taper_ratio_wing)/(3 + 3*taper_ratio_wing)       #MAC y location
    engine_1_percentage_Y = 0.315       # engine y location
    engine_2_percentage_Y = 0.5         # engine y location
    l_n1 = 0.25*mean_aerodynamic_chord + (MAC_y-wingspan_wing/2*engine_1_percentage_Y)*np.tan(sweep_quarter_chord) + 2.6
    l_n2 = 0.25*mean_aerodynamic_chord + (MAC_y-wingspan_wing/2*engine_2_percentage_Y)*np.tan(sweep_quarter_chord) + 2.6
    cla_a_h = alpha_lift_coefficient_aircraft_minus_tail(mach_number, aspect_ratio_wing, sweep_half_chord_wing,
                                                         surface_area_wing, wingspan_wing, taper_ratio_wing,
                                                         fuselage_width)
    x_ac_nacelles = 2*k_n*b_n**2*l_n1/(surface_area_wing*mean_aerodynamic_chord*cla_a_h) \
                    + 2*k_n*b_n**2*l_n2/(surface_area_wing*mean_aerodynamic_chord*cla_a_h)

    return x_ac_nacelles


def x_aerodynamic_center():
    x_ac_wing = 0.26        #/MAC
    x_ac = x_ac_wing + x_ac_fuselage_contribution_1(fuselage_width, fuselage_height, S_w, mean_aerodynamic_chord) \
           + x_ac_fuselage_contribution_2(surface_area_wing, wingspan_wing, taper_ratio_wing, fuselage_width, wingspan,
                                          mean_aerodynamic_chord, sweep_quarter_chord) \
           + x_ac_nacelles_contribution(wingspan_wing, taper_ratio_wing, surface_area_wing, mean_aerodynamic_chord,
                                        sweep_quarter_chord)
    return x_ac

