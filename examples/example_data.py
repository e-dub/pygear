# coding: utf-8

"""Example data for pyGear"""

# EXAMPLE GEARS
# external gears
extgear_1 = {'m_n': 5.0, 'z': 20, 'beta': 30.0, 'alpha_n': 20.0, 'x': 0.5, 'b': 50.0, 'z_2': 63, 'a': 245.0, 'h_k': 0.5,
             'd_Ff': 112.0, 'd_s': 60.0, 'A_s': -0.05}
extgear_2 = {'m_n': 5.0, 'z': 63, 'beta': -30.0, 'alpha_n': 20.0, 'b': 50.0, 'a': 245.0, 'h_k': 0.0, 'd_s': 130.0}
extgear_3 = {'m_n': 5.0, 'z': 20, 'beta': 0.0, 'alpha_n': 20.0, 'x': 0.5, 'b': 50.0, 'h_k': 0.5, 'rho_fP': 1.7,
             'A_s': -0.4}
extgear_4 = {'m_n': 5.0, 'z': 20, 'beta': 30.0, 'alpha_n': 20.0, 'x': 0.5, 'b': 50.0, 'z_2': 63, 'a': 245.0, 'h_k': 0.5,
             'd_Ff': 113.0, 'd_s': 60.0, 'A_s': -1.0}
extgear_5 = {'m_n': 5.0, 'z': 20, 'beta': 0.0, 'alpha_n': 20.0, 'x': -0.3, 'b': 50.0, 'z_2': 63, 'h_k': 0.5,
             'rho_fP': 1.0, 'A_s': -0.8}
extgear_6 = {'m_n': 8.0, 'z': 15, 'b': 50.0, 'A_s': -1.0}
extgear_7 = {'m_n': 5.0, 'z': 20, 'beta': 30.0, 'alpha_n': 20.0, 'x': -0.9, 'b': 50.0, 'z_2': 63, 'a': 245.0,
             'h_k': 0.5, 'd_s': 30.0, 'A_s': -1.0}
# internal gears
intgear_1 = {'m_n': 5.0, 'z': -100, 'b': 50.0, 'd_s': -550, 'h_k': 0.5, 'A_s': -1.5}
intgear_2 = {'m_n': 5.0, 'z': -80, 'beta': 20.0, 'b': 50.0, 'd_s': -550, 'h_k': 0.5, 'rho_f': 1.0, 'A_s': -0.4}
# gears for planetary stages
ringwheel_A = {'m_n': 16.0, 'z': -98, 'beta': 4.0, 'alpha_n': 20.0, 'b': 250.0, 'd_s': -1700.0}
sunwheel_A = {'m_n': 16.0, 'z': 30, 'beta': 4.0, 'alpha_n': 20.0, 'b': 250.0}
planetwheel_A = {'m_n': 16.0, 'z': 34, 'beta': -4.0, 'alpha_n': 20.0, 'b': 250.0}
ringwheel_B = {'m_n': 10.0, 'z': -140, 'beta': -8.0, 'alpha_n': 20.0, 'a': -460.0, 'b': 150.0, 'd_s': -1600.0}
sunwheel_B = {'m_n': 10.0, 'z': 25, 'beta': 8.0, 'alpha_n': 20.0, 'b': 200.0}
planetwheel_B1 = {'m_n': 10.0, 'z': 25, 'beta': 8.0, 'alpha_n': 20.0, 'b': 150.0}
planetwheel_B2 = {'m_n': 10.0, 'z': 90, 'beta': -8.0, 'x': 0.0, 'alpha_n': 20.0, 'b': 200.0}

ringwheel = {'m_n': 20.0, 'z': -119, 'beta': 0.0, 'alpha_n': 20.0, 'b': 300.0, 'd_s': -2500.0}
planetwheel_a1 = {'m_n': 20.0, 'z': 51, 'beta': 0.0, 'alpha_n': 20.0, 'b': 300.0}
planetwheel_a2 = {'m_n': 20.0, 'z': 17, 'beta': 0.0, 'alpha_n': 20.0, 'b': 300.0}

# EXAMPLE TOOLS
prottool = {'m': 5.0, 'h_aP0': 7.25, 'h_fP0': 5.5, 'alpha_P0': 20.0, 'alpha_prP0': 8.0, 'rho_aP0': 1.25, 'rho_fP0': 0.5,
            'pr_P0': 0.21, 's_P0': 7.3858}
stantool = {'m': 5.0, 'h_aP0': 7.25, 'h_fP0': 5.5, 'alpha_P0': 20.0, 'rho_aP0': 1.25, 'rho_fP0': 0.5}
chamtool = {'m': 5.0, 'h_aP0': 6.25, 'h_fP0': 6.0, 'alpha_P0': 20.0, 'rho_aP0': 1.00, 'h_FfP0': 4.8, 'alpha_KP0': 60.0}

stantool_small = {'m': 5.0E-1, 'h_aP0': 7.25E-1, 'h_fP0': 5.5E-1, 'alpha_P0': 20.0, 'rho_aP0': 1.25E-1,
                  'rho_fP0': 0.5E-1}
stantool_large = {'m': 5.0E+1, 'h_aP0': 7.25E+1, 'h_fP0': 5.5E+1, 'alpha_P0': 20.0, 'rho_aP0': 1.25E+1,
                  'rho_fP0': 0.5E+1}

# EXAMPLE GEAR MANUFACTURING MACHINES
hobber = {'z': 20, 'A_s': 0.0, 'x': 0.0, 'beta': 30.0}
