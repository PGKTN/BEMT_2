import calculator as ct

rotor_path = r'C:\BEMT\rotor.ini'

phi_arr = ct.get_phi(rotor_path)
alpha_arr = ct.get_alpha(rotor_path, phi_arr)

print(phi_arr)
print(alpha_arr)