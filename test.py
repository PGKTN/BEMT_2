import pandas as pd
import numpy as np

airfoil_path = r'C:\BEMT_2\airfoil.xlsx'
Reynolds = 15000
airfoil = 'naca0012'
coefficient = 'Cl'
alpha = 9.5

col_index = int(Reynolds / 10000) + 65

sheet_name = airfoil + '_' + coefficient
col1 = pd.read_excel(airfoil_path, sheet_name = sheet_name, usecols = chr(col_index))
col2 = pd.read_excel(airfoil_path, sheet_name = sheet_name, usecols = chr(col_index+1))

row1 = int(alpha+180)
row2 = int(alpha+180)+1

Co_11 = col1.iloc[row1, 0]
Co_12 = col1.iloc[row2, 0]
Co_21 = col2.iloc[row1, 0]
Co_22 = col2.iloc[row2, 0]

Rey1 = int(Reynolds / 10000) * 10000
Rey2 = int(Reynolds / 10000) * 10000 + 10000

num1 = Co_11*(np.abs(Reynolds-Rey1)/np.abs(Rey2-Rey1)) + Co_21*(np.abs(Reynolds-Rey2)/np.abs(Rey2-Rey1))
print('num1 = ', num1)
num2 = Co_12*(np.abs(Reynolds-Rey1)/np.abs(Rey2-Rey1)) + Co_22*(np.abs(Reynolds-Rey2)/np.abs(Rey2-Rey1))
print('num2 = ', num2)

alpha1 = int(alpha)
alpha2 = int(alpha)+1
Co = num1*(np.abs(alpha-alpha1)/np.abs(alpha2-alpha1)) + num2*(np.abs(alpha-alpha2)/np.abs(alpha2-alpha1))

print(Co)