import numpy as np
import configparser as cfp
import pandas as pd
import xlsxwriter as xw
import sys
import time
import threading

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

RPM = 3000
RPM_Thrust = []
RPM_Torque = []

while RPM < 8000:
    omega = RPM * rpm2omega
    epsilon = 1e-6
    Rtip = D / 2
    beta = 0
    alpha0 = 0
    V_for = Vinf
    Vc = 0

    y = np.zeros(Nsec-1)
    for i in range(len(y)):
        y[i] = R[i+1] - R[i]

    # Rsec = np.zeros(Nsec-1)
    # for i in range(len(Rsec)):
    #     Rsec[i] = Rhub + y[i]
    # print('Rsec = ', Rsec)

    dA = np.zeros(Nsec)
    for i in range(len(dA)):
        if i == 0:
            dA[i] = 0
        else:
            dA[i] = Chord[i-1] * y[i-1]

    Asec = np.zeros(Nsec)
    for i in range(len(Asec)):
        for j in range(i):
            Asec[i] += dA[j+1]
    A = Asec[-1]

    sigma = np.zeros(Nsec)
    for i in range(len(sigma)):
        sigma[i] = (Nb * Asec[i]) / (np.pi * R[i]**2)

    Theta_star = np.zeros(Nsec)
    for i in range(len(Theta_star)):
        Theta_star[i] = (beta + Theta[i] - alpha0) * deg2rad

    alpha_left = np.zeros(Nsec)
    alpha_right = np.zeros(Nsec)
    Reynolds = np.zeros(Nsec)
    alpha = np.zeros(Nsec)
    phi = np.zeros(Nsec)
    gamma = np.zeros(Nsec)
    Cl = np.zeros(Nsec)
    Cd = np.zeros(Nsec)
    dT = np.zeros(Nsec)
    dQ = np.zeros(Nsec)

    T = 0
    Q = 0

    stop_event = threading.Event()
    dot_thread = threading.Thread(target=rd.print_dots, args=(stop_event,))
    dot_thread.start()

    slope_Cl = np.pi * 2
    g_left = 0
    g_right = 1

    for i in range(Nsec):
        phi_left = 0
        phi_right = np.pi / 2
        
        Reynolds[i] = rho * omega * R[i] * Chord[i] / mu
        
        Cl_reynolds = rd.Get_Airfoil_Xlsx(airfoil_path, section[i], 'Cl')
        Cd_reynolds = rd.Get_Airfoil_Xlsx(airfoil_path, section[i], 'Cd')
        
        V_rat = V_for / (omega * R[i])
        
        while True:
            # alpha_left = Theta_star[i] - phi_left
            # Cl_left = rd.Get_Coefficient(Cl_reynolds, Reynolds[i], alpha_left*rad2deg)
            Cd_phi_left = rd.Get_Coefficient(Cd_reynolds, Reynolds[i], phi_left*rad2deg)
            g_left = ct.g_Solver(phi_left, sigma[i], slope_Cl, Theta_star[i], Cd_phi_left, V_rat)
            
            # alpha_right = Theta_star[i] - phi_right
            # Cl_right = rd.Get_Coefficient(Cl_reynolds, Reynolds[i], alpha_right*rad2deg)
            Cd_phi_right = rd.Get_Coefficient(Cd_reynolds, Reynolds[i], phi_right*rad2deg)
            g_right = ct.g_Solver(phi_right, sigma[i], slope_Cl, Theta_star[i], Cd_phi_right, V_rat)
            # print('g_left = ', g_left)
            # print('g_right = ', g_right)
            if g_left == 0 or g_right == 0:
                phi[i] = 0
                alpha[i] = 0
                Cl[i] = 0
                Cd[i] = 0
                break
                
            elif g_left * g_right > 0:
                print('g_left = ', g_left)
                print('g_right = ', g_right)
                print('section ', i, '번째 오류, ', 'g_left * g_right > 0')
                exit()
                
            elif g_left * g_right < 0 and np.abs(g_left + g_right) < epsilon:
                phi[i] = phi_new
                alpha[i] = Theta_star[i] - phi[i]
                Cl[i] = rd.Get_Coefficient(Cl_reynolds, Reynolds[i], alpha[i]*rad2deg)
                Cd[i] = rd.Get_Coefficient(Cd_reynolds, Reynolds[i], alpha[i]*rad2deg)
                break
                
            elif g_left * g_right < 0:
                phi_new = 0.5 * (phi_left + phi_right)
                # phi_new = phi_left - (phi_right - phi_left)*(g_left)/(g_right - g_left)
                Cd_phi_new = rd.Get_Coefficient(Cd_reynolds, Reynolds[i], phi_new*rad2deg)
                g_new = ct.g_Solver(phi_new, sigma[i], slope_Cl, Theta_star[i], Cd_phi_new, V_rat)
                if g_new * g_left > 0:
                    # g_left = g_new
                    # phi_new = phi_left - (phi_right - phi_left)*(g_left)/(g_right - g_left)
                    phi_left = phi_new
                else:
                    # phi_new = phi_left - (phi_right - phi_left)*(g_left)/(g_right - g_left)
                    phi_right = phi_new
                    
            else:
                pass
        if Cl[i] != 0:
            gamma[i] = np.arctan(Cd[i] / Cl[i])
        else:
            gamma[i] = 0
            
        G_phi = np.cos(phi[i]) + 0.25 * sigma[i] * Cl[i] * np.sin(phi[i]+gamma[i]) / (np.cos(gamma[i] * np.sin(phi[i])))
        Vr = omega*R[i]/G_phi
        
        if i == 0:
            dT[i] = 0
            dQ[i] = 0
        else:
            dT[i] = 0.5 * Nb * Chord[i-1] * rho * Vr**2 * (Cl[i-1]*np.cos(phi[i-1]) - Cd[i-1]*np.sin(phi[i-1])) * y[i-1]
            dQ[i] = 0.5 * Nb * Chord[i-1] * rho * Vr**2 * (Cl[i-1]*np.sin(phi[i-1]) + Cd[i-1]*np.cos(phi[i-1])) * R[i-1] * y[i-1]
        
        T += dT[i]
        Q += dQ[i]
        
        print(' ')
        print(f'section {i}번째 완료')
    RPM_Thrust.append(T)
    RPM_Torque.append(Q)
    
    stop_event.set()
    dot_thread.join()

    vi = np.sqrt((T)/(2 * A * rho))
    J = vi / (RPM * D)
    
    RPM += 500


for i in range(len(RPM_Thrust)):
    print(f'RPM = {3000+i*500}, Thrust = {RPM_Thrust[i]}, Torque = {RPM_Torque[i]}')
    
# print('phi = ', phi*rad2deg)
# print('alpha = ', alpha*rad2deg)
# print('Cl = ', Cl)
# print('Cd = ', Cd)
# print('vi = ', vi)
# print('J = ', J)
# print('dT = ', dT)
# print('dQ = ', dQ)
# print('Thrust = ', T)
# print('Torque = ', Q)


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