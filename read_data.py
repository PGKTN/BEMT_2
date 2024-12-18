import numpy as np
import configparser as cfp
import pandas as pd

def get_rotor_data(rotor_path):
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
    ini_option_pitch = 'pitch'

    ini_section_fluid = 'fluid'
    ini_option_rho = 'rho'
    ini_option_mu = 'mu'

    # ini파일에서 값을 읽어와 string -> float로 변환.
    rpm = float(config.get(ini_section_case, ini_option_rpm))
    v_inf = float(config.get(ini_section_case, ini_option_v_inf))
    nblades = float(config.get(ini_section_rotor, ini_option_nblades))
    diameter = float(config.get(ini_section_rotor, ini_option_diameter))
    radius_hub = float(config.get(ini_section_rotor, ini_option_radius_hub))
    rho = float(config.get(ini_section_fluid, ini_option_rho))
    mu = float(config.get(ini_section_fluid, ini_option_mu))
    
    # ini파일에서 값을 읽어와 string -> list로 변환.
    section = config.get(ini_section_rotor, ini_option_section).split(' ')

    # section의 개수.
    sections_count = len(section)

    # radius, chord, pitch 값을 split으로 나눔과 동시에 float형으로 변환
    radius_ = [float(value) for value in config.get(ini_section_rotor, ini_option_radius).split(' ')]
    chord_ = [float(value) for value in config.get(ini_section_rotor, ini_option_chord).split(' ')]
    pitch_ = [float(value) for value in config.get(ini_section_rotor, ini_option_pitch).split(' ')]

    radius = np.array(radius_)
    chord = np.array(chord_)
    pitch = np.array(pitch_)

    # section, radius, chord, pitch의 개수가 일치하는지 확인.
    if sections_count != len(radius) or sections_count != len(chord) or sections_count != len(pitch):
        raise ValueError('section, radius, chord, pitch의 개수가 일치하지 않습니다.')
    else :
        print('sections_count 정상.')
        return rpm, v_inf, nblades, diameter, radius_hub, section, radius, chord, pitch, rho, mu

def Get_Airfoil_Xlsx(airfoil_path, airfoil, coefficient):
    sheet_name = airfoil + '_' + coefficient
    
    df = pd.read_excel(airfoil_path, sheet_name = sheet_name)
    Ncol = df.shape[1]
    
    start_num = 66
    
    Co_reynolds = []
    for i in range(Ncol-1):
        if i < 25:
            col_data = pd.read_excel(airfoil_path, sheet_name = sheet_name, usecols = chr(start_num+i))
        else:
            col_data = pd.read_excel(airfoil_path, sheet_name = sheet_name, usecols = 'A' + chr(start_num+i-26))
        Co_reynolds.append(col_data.values.flatten())
    Co_reynolds = np.array(Co_reynolds)
    return Co_reynolds

def Get_Coefficient(Co_reynolds, Reynolds, alpha):
    col_index = int(Reynolds / 10000) - 1
    col1 = Co_reynolds[col_index]
    col2 = Co_reynolds[col_index+1]
    
    row1 = int(alpha+180)
    row2 = int(alpha+180)+1
    
    Co_11 = col1[row1]
    Co_12 = col1[row2]
    Co_21 = col2[row1]
    Co_22 = col2[row2]
    
    Rey1 = int(Reynolds / 10000) * 10000
    Rey2 = int(Reynolds / 10000) * 10000 + 10000
    
    num1 = Co_11*(np.abs(Reynolds-Rey1)/np.abs(Rey2-Rey1)) + Co_21*(np.abs(Reynolds-Rey2)/np.abs(Rey2-Rey1))
    num2 = Co_12*(np.abs(Reynolds-Rey1)/np.abs(Rey2-Rey1)) + Co_22*(np.abs(Reynolds-Rey2)/np.abs(Rey2-Rey1))
    
    alpha1 = int(alpha)
    alpha2 = int(alpha)+1
    Co = num1*(np.abs(alpha-alpha1)/np.abs(alpha2-alpha1)) + num2*(np.abs(alpha-alpha2)/np.abs(alpha2-alpha1))
    
    return Co

import sys
import time

def print_dots(stop_event):
    dots = ['.', '..', '...', '....']
    while not stop_event.is_set():
        for dot in dots:
            if stop_event.is_set():
                break
            sys.stdout.write(f'\r{dot}')
            sys.stdout.flush()
            time.sleep(0.1)
        sys.stdout.write('\r    \r')
