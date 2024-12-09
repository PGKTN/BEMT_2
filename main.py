import numpy as np
import configparser as cfp
import pandas as pd
import xlsxwriter as xw

import calculator as ct
import read_data as rd
import airfoil_storage as af

from xfoil import XFoil

rpm2omega = 2 * np.pi / 60
rad2deg = 180 / np.pi
deg2rad = np.pi / 180

rotor_path = r'C:\BEMT_2\rotor.ini'

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
section = config.get(ini_section_rotor, ini_option_section).split(' ')

# section의 개수.
Nsec = len(section)

# radius, chord, theta 값을 split으로 나눔과 동시에 float형으로 변환
radius_ = [float(value) for value in config.get(ini_section_rotor, ini_option_radius).split(' ')]
chord_ = [float(value) for value in config.get(ini_section_rotor, ini_option_chord).split(' ')]
theta_ = [float(value) for value in config.get(ini_section_rotor, ini_option_theta).split(' ')]

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

y = np.zeros(Nsec-1)
for i in range(len(y)):
    y[i] = R[i+1] - R[i]
print('y = ', y)

# Rsec = np.zeros(Nsec-1)
# for i in range(len(Rsec)):
#     Rsec[i] = Rhub + y[i]
# print('Rsec = ', Rsec)

dA = np.zeros(Nsec-1)
for i in range(len(dA)):
    dA[i] = Chord[i] * y[i]

Asec = np.zeros(Nsec-1)
for i in range(len(Asec)):
    for j in range(i+1):
        Asec[i] += dA[j]
A = Asec[Nsec-2]

Sigma = np.zeros(Nsec-1)
for i in range(len(Sigma)):
    Sigma[i] = (Nb * Asec[i]) / (np.pi * R[i+1]**2)

# Vrat = np.zeros(Nsec)
# for i in range(len(Vrat)):
#     Vrat[i] = Vinf / (omega * R[i])

Theta_star = np.zeros(Nsec)
for i in range(len(Theta_star)):
    Theta_star[i] = (beta + Theta[i] - alpha0)

Phi_result = np.zeros(Nsec)
Reynolds = np.zeros(Nsec)
VR = np.zeros(Nsec)
G = np.zeros(Nsec)
alpha = np.zeros(Nsec)
UP = np.zeros(Nsec)
UT = np.zeros(Nsec)
U = np.zeros(Nsec)
Cl = np.zeros(Nsec)
Cd = np.zeros(Nsec)
dL = np.zeros(Nsec)
dD = np.zeros(Nsec)
dT = np.zeros(Nsec)
dQ = np.zeros(Nsec)
phi = np.zeros(Nsec)
dFx = np.zeros(Nsec)
dFz = np.zeros(Nsec)

T_new = 0
T_old = 10000
T = 0

Vc = 0
vi = 0.0001

while np.abs(T_new - T_old) > 1e-2:
    T_old = T_new
    T_new = 0
    
    for i in range(Nsec-1):
        xf = XFoil()
        xf.airfoil = af.naca0012
        xf.max_iter = 100
        Reynolds[i] = rho * omega * R[i] * Chord[i] / mu
        xf.Re = Reynolds[i]
        
        UP[i] = Vc + vi
        UT[i] = omega * R[i]
        U[i] = np.sqrt(UT[i]**2 + UP[i]**2)
        
        phi[i] = np.arctan(UP[i]/UT[i])
        alpha[i] = Theta_star[i] - phi[i]
        
        Cl_, Cd_, Cm_, Cp_ = xf.a(alpha[i])
        while np.isnan(Cl_) or np.isnan(Cd_):
            Reynolds[i] += 1000
            xf.Re = Reynolds[i]
            Cl_, Cd_, Cm_, Cp_ = xf.a(alpha[i])
        Cl[i] = Cl_
        Cd[i] = Cd_
        
        dL[i] = 0.5 * rho * U[i]**2 * Chord[i] * Cl[i] * y[i]
        dD[i] = 0.5 * rho * U[i]**2 * Chord[i] * Cd[i] * y[i]
        
        dFx[i] = dL[i] * np.sin(phi[i]) + dD[i] * np.cos(phi[i])
        dFz[i] = dL[i] * np.cos(phi[i]) - dD[i] * np.sin(phi[i])
        
        dT[i] = Nb * dFz[i]
        dQ[i] = Nb * dFx[i] * R[i]
        
        T_new += dT[i]
    
    vi = np.sqrt(T_new/(2*rho*A))
    T = T_new
    print('T_new = ', T_new)
    print('T_old = ', T_old)

for i in range(Nsec-1):
    print('Reynolds[', i, '] = ', Reynolds[i])
    print('theta[', i, '] = ', Theta[i])
    print('phi[', i, '] = ', phi[i])
    print('alpha[', i, '] = ', alpha[i])
    print('UP[', i, '] = ', UP[i])
    print('UT[', i, '] = ', UT[i])
    print('U[', i, '] = ', U[i])
    print('Cl[', i, '] = ', Cl[i])
    print('Cd[', i, '] = ', Cd[i])
    print('dL[', i, '] = ', dL[i])
    print('dD[', i, '] = ', dD[i])
    print('dT[', i, '] = ', dT[i])
    print('dQ[', i, '] = ', dQ[i])

print('Thrust = ', T)
# for i in range(Nsec):
#     xf = XFoil()
#     xf.airfoil = af.naca0012
#     xf.max_iter = 100
#     xf.Re = rho * omega * R[i] * Chord[i] / mu
    
#     Cl, Cd, Cm, Cp = xf.a(Theta_star[i])
#     while(np.isnan(Cl) or np.isnan(Cd)):
#         xf.Re += xf.Re * 1000
#         Cl, Cd, Cm, Cp = xf.a(Theta_star[i])

#     SigCl = 0.25 * Sigma[i] * Cl
#     SigCd = 0.25 * Sigma[i] * Cd
    
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