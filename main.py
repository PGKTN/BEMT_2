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

y = np.zeros(Nsec)
for i in range(len(y)):
    if(i==0):
        y[i] = Rhub
    else:
        y[i] = R[i] - R[i-1]
print('y = ', y)

# Rsec = np.zeros(Nsec-1)
# for i in range(len(Rsec)):
#     Rsec[i] = Rhub + y[i]
# print('Rsec = ', Rsec)

dA = np.zeros(Nsec)
for i in range(len(dA)):
    if(i==0):
        dA[i] = 0
    else:
        dA[i] = Chord[i] * y[i]
print('dA = ', dA)

Asec = np.zeros(Nsec)
for i in range(len(Asec)):
    for j in range(i+1):
        Asec[i] += dA[j]
print('Asec = ', Asec)

Sigma = np.zeros(Nsec)
for i in range(Nsec):
    Sigma[i] = (Nb * Asec[i]) / (np.pi * R[i]**2)
print('Sigma = ', Sigma)

Vrat = np.zeros(Nsec)
for i in range(len(Vrat)):
    Vrat[i] = Vinf / (omega * R[i])
print('Vrat = ', Vrat)

Theta_star = np.zeros(Nsec)
for i in range(len(Theta_star)):
    Theta_star[i] = (beta + Theta[i] - alpha0)
print('Theta_star = ', Theta_star)

Phi_result = np.zeros(Nsec)
VR = np.zeros(Nsec)

for i in range(Nsec):
    xf = XFoil()
    xf.airfoil = af.naca0012
    xf.max_iter = 100
    xf.Re = rho * omega * R[i] * Chord[i] / mu
    print('xf.Re = ', xf.Re)
    
    Cl, Cd, Cm, Cp = xf.a(Theta_star[i])
    while(np.isnan(Cl) or np.isnan(Cd)):
        xf.Re += 1000
        Cl, Cd, Cm, Cp = xf.a(Theta_star[i])

    SigCl = 0.25 * Sigma[i] * Cl
    SigCd = 0.25 * Sigma[i] * Cd
    
    G_left = -(SigCl * (Theta_star[i]*deg2rad) + Vrat[i] * SigCd)
    G_right = 1 + SigCd + Vrat[i]*SigCl*(np.pi/2 - (Theta_star[i]*deg2rad))
    
    Phi_left = 0
    Phi_right = np.pi / 2
    
    G_new = 1
    while abs(G_new) > epsilon:
        DelPhi = Phi_right - Phi_left
        DelG = G_right - G_left
        Phi_new = Phi_left - (DelPhi / DelG) * G_left
        
        SinPhi = np.sin(Phi_new)
        CosPhi = np.cos(Phi_new)
        H_phi = SinPhi + SigCd
        E_phi = SigCl * ((Theta_star[i]*deg2rad) - Phi_new)
        F_SinPhi = H_phi * SinPhi - E_phi * CosPhi
        G_new = F_SinPhi - Vrat[i]*(H_phi * CosPhi + E_phi * SinPhi)
        
        if G_left * G_new < 0:
            Phi_right = Phi_new
            G_right = G_new
        else:
            Phi_left = Phi_new
            G_left = G_new
            
        Phi_result[i] = Phi_new * rad2deg
    VR[i] = omega * R[i] / G_new

alpha = np.zeros(Nsec)
for i in range(len(alpha)):
    alpha[i] = Theta_star[i] - Phi_result[i]
print('alpha = ', alpha)

dT = np.zeros(Nsec)
dQ = np.zeros(Nsec)
for i in range(Nsec):
    xf = XFoil()
    xf.airfoil = af.naca0012
    xf.max_iter = 100
    xf.Re = rho * omega * R[i] * Chord[i] / mu
    print('xf.Re = ', xf.Re)
    # breakpoint()
    
    Cl, Cd ,Cm, Cp = xf.a(alpha[i])
    while(np.isnan(Cl) or np.isnan(Cd)):
        xf.Re += 1000
        Cl, Cd, Cm, Cp = xf.a(alpha[i])
    print('Cl = ', Cl)
    print('Cd = ', Cd)
    # breakpoint()
    
    SinPhi = np.sin(Phi_result[i]*deg2rad)
    CosPhi = np.cos(Phi_result[i]*deg2rad)
    SigCl = 0.25 * Sigma[i] * Cl
    SigCd = 0.25 * Sigma[i] * Cd
    G = ((CosPhi * (SinPhi + SigCd)) + SigCl*SinPhi)/(SinPhi)
    
    VR = omega * R[i] / G
    if np.isnan(VR):
        VR = 0
    print('VR = ', VR)
    # breakpoint()
    
    dT[i] = 0.5 * Nb * Chord[i] * rho * VR**2 * (Cl*CosPhi - Cd*SinPhi) * y[i]
    dQ[i] = 0.5 * Nb * Chord[i] * rho * VR**2 * (Cl*SinPhi - Cd*CosPhi) * R[i] * y[i]
    print('Nb = ', Nb)
    print('Chord[', i, '] = ', Chord[i])
    print('rho = ', rho)
    print('dT[', i, '] = ', dT[i])
    print('dQ[', i, '] = ', dQ[i])
    # breakpoint()
    
Thrust = 0
Torque = 0
for i in range(Nsec):
    Thrust += 0.5 * dT[i]
    Torque += 0.5 * dQ[i]

print('Thrust = ', Thrust)
print('Torque = ', Torque)