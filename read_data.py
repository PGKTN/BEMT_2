import numpy as np
import configparser as cfp

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