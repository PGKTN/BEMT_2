import numpy as np
import configparser as cfp
import pandas as pd
import xlsxwriter as xw

import calculator as ct
import read_data as rd

rpm2omega = 2 * np.pi / 60
rad2deg = 180 / np.pi
deg2rad = np.pi / 180

rotor_path = r'C:\BEMT_2\rotor.ini'
airfoil_path = r'C:\BEMT_2\airfoil.xlsx'

config = cfp.ConfigParser()
config.read(rotor_path)

ini_section_case = 'case'
ini_option_rpm = 'rpm'
ini_option_v_inf = 'v_inf'

ini_section_rotor = 'rotor'
ini_option_nblades = 'nblades'
ini_option_diameter = 'diameter'
ini_option_radius_hub = 'radius_hub'
ini_option_section = 'section'
ini_option_radius = 'radius'
ini_option_chord = 'chord'
ini_option_theta = 'pitch'

ini_section_fluid = 'fluid'
ini_option_rho = 'rho'
ini_option_mu = 'mu'

# ini파일에서 값을 읽어와 string -> float로 변환.
RPM = float(config.get(ini_section_case, ini_option_rpm))
Vinf = float(config.get(ini_section_case, ini_option_v_inf))
Nb = float(config.get(ini_section_rotor, ini_option_nblades))
D = float(config.get(ini_section_rotor, ini_option_diameter))
Rhub = float(config.get(ini_section_rotor, ini_option_radius_hub))
rho = float(config.get(ini_section_fluid, ini_option_rho))
mu = float(config.get(ini_section_fluid, ini_option_mu))

# ini파일에서 값을 읽어와 string -> list로 변환.
sections = config.get(ini_section_rotor, ini_option_section).split(' ')

# section의 개수.
Nsec = len(sections)

# radius, chord, theta 값을 split으로 나눔과 동시에 float형으로 변환
radius_ = [float(value) for value in config.get(ini_section_rotor, ini_option_radius).split(' ')]
chord_ = [float(value) for value in config.get(ini_section_rotor, ini_option_chord).split(' ')]
theta_ = [float(value) for value in config.get(ini_section_rotor, ini_option_theta).split(' ')]

# 순서가 없는 list형식에서 순서가 있는 array형식으로 변환.
section = np.array(sections)
R = np.array(radius_)
Chord= np.array(chord_)
Theta = np.array(theta_)

# Section, R, Chord, Theta의 개수가 일치하는지 확인.
if Nsec != len(R) or Nsec != len(Chord) or Nsec != len(Theta):
    raise ValueError('section, radius, chord, theta의 개수가 일치하지 않습니다.')
else :
    print('Sections 정상.')

omega = RPM * rpm2omega
epsilon = 1e-6
Rtip = D / 2
Slope_Cl = 2 * np.pi
beta = 0
alpha0 = 0

y = np.zeros(Nsec)
for i in range(len(y)):
    if i == 0:
        y[i] = Rhub
    else:
        y[i] = R[i] - R[i-1]

# Rsec = np.zeros(Nsec-1)
# for i in range(len(Rsec)):
#     Rsec[i] = Rhub + y[i]
# print('Rsec = ', Rsec)

dA = np.zeros(Nsec)
for i in range(len(dA)):
    if i == 0:
        dA[i] = 0
    else:
        dA[i] = Chord[i] * y[i]

Asec = np.zeros(Nsec)
for i in range(len(Asec)):
    for j in range(i):
        Asec[i] += dA[j+1]
A = Asec[Nsec-1]

sigma = np.zeros(Nsec)
for i in range(len(sigma)):
    if i == 0:
        sigma[i] = 0
    else:
        sigma[i] = (Nb * Asec[i]) / (np.pi * R[i]**2)

Theta_star = np.zeros(Nsec)
for i in range(len(Theta_star)):
    Theta_star[i] = (beta + Theta[i] - alpha0) * deg2rad

gamma_left = np.zeros(Nsec)
gamma_right = np.zeros(Nsec)
alpha_left = np.zeros(Nsec)
alpha_right = np.zeros(Nsec)
Reynolds = np.zeros(Nsec)
alpha = np.zeros(Nsec)
phi = np.zeros(Nsec)
Cl = np.zeros(Nsec)
Cd = np.zeros(Nsec)


slope_Cl = np.pi * 2
f_left = 0
f_right = 1

for i in range(Nsec):
    phi_left = 0
    phi_right = np.pi / 2
    phi_new = 0
    alpha_new = 0
    Cl_new = 0
    Cd_new = 0
    
    Reynolds[i] = rho * omega * R[i] * Chord[i] / mu
    
    f_not_equal_zero = True
    num = 0
    while f_not_equal_zero:
        alpha_left = Theta_star[i] - phi_left
        Cl_left = rd.Get_Coefficient(airfoil_path, section[i], 'Cl', Reynolds[i], alpha_left*rad2deg)
        Cd_left = rd.Get_Coefficient(airfoil_path, section[i], 'Cd', Reynolds[i], alpha_left*rad2deg)
        f_left = ct.f_Solver(phi_left, sigma[i], Cl_left, Cd_left)
        
        alpha_right = Theta_star[i] - phi_right
        Cl_right = rd.Get_Coefficient(airfoil_path, section[i], 'Cl', Reynolds[i], alpha_right*rad2deg)
        Cd_right = rd.Get_Coefficient(airfoil_path, section[i], 'Cd', Reynolds[i], alpha_right*rad2deg)
        f_right = ct.f_Solver(phi_right, sigma[i], Cl_right, Cd_right)
        print(num, '번째')
        print('f_left = ', f_left)
        print('f_right = ', f_right)
        num += 1
        if np.abs(f_left - f_right) < epsilon:
            phi[i] = phi_new
            alpha[i] = alpha_new
            Cl[i] = Cl_new
            Cd[i] = Cd_new
            f_not_equal_zero = False
        elif f_left * f_right > 0:
            ValueError('f_left * f_right > 0')
            exit()
        elif f_left == 0 or f_right == 0:
            phi[i] = 0
            alpha[i] = Theta_star[i] - phi[i]
            Cl[i] = 0
            Cd[i] = 0
            f_not_equal_zero = False
        else:
            phi_new = phi_left - (phi_right - phi_left)*(f_left)/(f_right - f_left)
            alpha_new = Theta_star[i] - phi_new
            Cl_new = rd.Get_Coefficient(airfoil_path, section[i], 'Cl', Reynolds[i], alpha_new*rad2deg)
            Cd_new = rd.Get_Coefficient(airfoil_path, section[i], 'Cd', Reynolds[i], alpha_new*rad2deg)
            f_new = ct.f_Solver(phi_new, sigma[i], Cl_new, Cd_new)
            if f_new * f_left > 0:
                phi_left = phi_new
                f_left = f_new
            else:
                phi_right = phi_new
                f_right = f_new

print('phi = ', phi*rad2deg)
print('alpha = ', alpha*rad2deg)
print('Cl = ', Cl)
print('Cd = ', Cd)

# for i in range(Nsec):
#     xf = XFoil()
#     xf.airfoil = af.naca0012
#     xf.max_iter = 100
#     xf.Re = rho * omega * R[i] * Chord[i] / mu
    
#     Cl, Cd, Cm, Cp = xf.a(Theta_star[i])
#     while(np.isnan(Cl) or np.isnan(Cd)):
#         xf.Re += xf.Re * 1000
#         Cl, Cd, Cm, Cp = xf.a(Theta_star[i])

#     SigCl = 0.25 * sigma[i] * Cl
#     SigCd = 0.25 * sigma[i] * Cd
    
#     G_left = -(SigCl * (Theta_star[i]*deg2rad) + Vrat[i] * SigCd)
#     G_right = 1 + SigCd + Vrat[i]*SigCl*(np.pi/2 - (Theta_star[i]*deg2rad))
    
#     Phi_left = 0
#     Phi_right = np.pi / 2
    
#     G_new = 1
#     while abs(G_new) > epsilon:
#         DelPhi = Phi_right - Phi_left
#         DelG = G_right - G_left
#         Phi_new = Phi_left - (DelPhi / DelG) * G_left
        
#         SinPhi = np.sin(Phi_new)
#         CosPhi = np.cos(Phi_new)
#         H_phi = SinPhi + SigCd
#         E_phi = SigCl * ((Theta_star[i]*deg2rad) - Phi_new)
#         F_SinPhi = H_phi * SinPhi - E_phi * CosPhi
#         G_new = F_SinPhi - Vrat[i]*(H_phi * CosPhi + E_phi * SinPhi)
        
#         if G_left * G_new < 0:
#             Phi_right = Phi_new
#             G_right = G_new
#         else:
#             Phi_left = Phi_new
#             G_left = G_new
            
#         Phi_result[i] = Phi_new * rad2deg
#     VR[i] = omega * R[i] / G_new
#     G[i] = G_new
    
#     alpha[i] = Theta_star[i] - Phi_result[i]
    
#     Cl, Cd, Cm, Cp = xf.a(alpha[i])
#     Reynolds[i] = rho * omega * R[i] * Chord[i] / mu
#     while(np.isnan(Cl) or np.isnan(Cd)):
#         Reynolds[i] += 1000
#         xf.Re = Reynolds[i]
#         Cl, Cd, Cm, Cp = xf.a(alpha[i])
        
#     SinPhi = np.sin(Phi_result[i]*deg2rad)
#     CosPhi = np.cos(Phi_result[i]*deg2rad)
    
#     Cn = Cl*CosPhi - Cd*SinPhi
#     Ct = Cl*SinPhi + Cd*CosPhi
    
#     Bcro = 0.5 * Nb * Chord[i] * rho * VR[i]**2
#     Cl = Slope_Cl * (Theta_star - Phi_result[i]) * deg2rad
    
#     dT[i] = Bcro * Cn * y[i]
#     dQ[i] = Bcro * Ct * R[i] * y[i]
    
#     if np.isnan(dT[i]) or np.isinf(dT[i]):
#         dT[i] = 0
#     if np.isnan(dQ[i]) or np.isinf(dQ[i]):
#         dQ[i] = 0

# for i in range(Nsec):
#     Thrust += 0.5 * dT[i]
#     Torque += 0.5 * dQ[i]

# print('Thrust = ', Thrust)
# print('Torque = ', Torque)