import pandas as pd
import numpy as np

airfoil = 'naca0012'
coefficient = 'Cl'
airfoil_path = r'C:\BEMT_2\airfoil.xlsx'
Reynolds = 15000
alpha = 9.5

sheet_name = airfoil + '_' + coefficient
df = pd.read_excel(airfoil_path, sheet_name = sheet_name)
Ncol = df.shape[1]

start_num = 66

i = 0

Cl_reynolds = []
for i in range(Ncol-1):
    col_data = pd.read_excel(airfoil_path, sheet_name = sheet_name, usecols = chr(start_num+i))
    Cl_reynolds.append(col_data.values.flatten())

Cl_reynolds = np.array(Cl_reynolds)

col_index = int(Reynolds / 10000) - 1
col1 = Cl_reynolds[col_index]
col2 = Cl_reynolds[col_index+1]

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

print('Co11 = ', Co_11)
print('Co12 = ', Co_12)
print('Co21 = ', Co_21)
print('Co22 = ', Co_22)
print('Co = ', Co)