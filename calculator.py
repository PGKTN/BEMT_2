import numpy as np
import configparser as cfp
import sys

import read_data as rd
import airfoil_storage as af

rpm2omega = 2 * np.pi / 60
rad2deg = 180 / np.pi
deg2rad = np.pi / 180

# def get_alpha(rotor_path, phi_arr):
#     rpm, v_inf, Nb, diameter, radius_hub, section, radius, chord, pitch, rho, mu = rd.get_rotor_data(rotor_path)
#     alpha_arr = np.zeros(len(section))
#     if len(alpha_arr) != len(phi_arr):
#         raise ValueError('alpha_arr과 phi_arr의 길이가 일치하지 않습니다.')
#     else:
#         for i in range(len(section)):
#             alpha_arr[i] = pitch[i] - phi_arr[i]
#     return alpha_arr

# def get_phi(rotor_path):
#     rpm, v_inf, Nb, diameter, radius_hub, section, radius, chord, pitch, rho, mu = rd.get_rotor_data(rotor_path)

#     phi_new = 0
#     phi_arr = np.zeros(len(section))
    
#     Cl_alpha = 2 * np.pi
#     alpha_zero = 0
#     beta = 0
    
#     for i in range(len(section)):
#         phi1 = 0
#         phi2 = np.pi/2
        
#         f1 = f_solver(phi1, rpm, Nb, section[i], radius[i], chord[i], pitch[i], Cl_alpha=Cl_alpha, alpha_zero=alpha_zero, beta=beta, rho=rho, mu=mu)
#         f2 = f_solver(phi2, rpm, Nb, section[i], radius[i], chord[i], pitch[i], Cl_alpha=Cl_alpha, alpha_zero=alpha_zero, beta=beta, rho=rho, mu=mu)
        
#         print('phi1 = ', phi1)
#         print('phi2 = ', phi2)
#         print('f1 = ', f1)
#         print('f2 = ', f2)
#         # breakpoint()
        
#         if f1 * f2 <= 0:
#             while f1 * f2 <= 0:
#                 phi_new = 0.5 * (phi1 + phi2)
#                 f_new = f_solver(phi_new, rpm, Nb, section[i], radius[i], chord[i], pitch[i], Cl_alpha=Cl_alpha, alpha_zero=alpha_zero, beta=beta, rho=rho, mu=mu)
#                 if f_new * f1 > 0:
#                     phi1 = phi_new
#                     f1 = f_new
#                 else:
#                     phi2 = phi_new
#                     f2 = f_new
#                 if abs(f1) < 1e-6 and abs(f2) < 1e-6:
#                     break
#             phi_arr[i] = phi_new * rad2deg
#             print('phi_arr[', i, '] = ', phi_arr[i])
#             # breakpoint()
#         else:
#             print('f1: ', f1)
#             print('f2: ', f2)
#             print('f1 * f2 > 0')
#             sys.exit()
#     return phi_arr

def g_Solver(phi, sigma, slope_Cl, Theta_star, Cd_phi, V_rat):
    
    H = np.sin(phi) + 0.25 * sigma * Cd_phi
    E = 0.25 * sigma * slope_Cl*(Theta_star - phi)
    g = (H * np.sin(phi) - E * np.cos(phi)) - V_rat * (H * np.cos(phi) + E * np.sin(phi))
    
    return g