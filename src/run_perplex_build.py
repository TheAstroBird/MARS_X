import subprocess as sub
import pandas as pd
import numpy as np

composition = 'T13'

mass_amount_data = pd.read_excel('../../chemical_models.xlsx', index_col=0)
therm_components = ['NA2O', 'MGO', 'AL2O3', 'SIO2', 'CAO', 'FEO']
if composition in ['BF97', 'MA79', 'LF97', 'S99', 'KC08', 'T13']:
    name = 'mars_' + composition.lower()
    i = 0
    while i < len(therm_components):
        if np.isnan(mass_amount_data[composition][therm_components[i]]):
            therm_components.pop(i)
        else:
            i += 1
    mass_amount = [str(mass_amount_data[composition][comp]) for comp in therm_components]
else:
    print("ERROR! There isn't such composition model")
    quit()
areotherm = [1.240095e3, 8.085953e-3, -1.774443e-9, -1.746166e-13, 4.278426e-19]
areotherm = [str(coeff) for coeff in areotherm]
stx21_solution_models = ['C2/c', 'Wus', 'Pv', 'Pl', 'Sp', 'O', 'Wad', 'Ring', 'Opx', 'Cpx',
                         'Aki', 'Gt', 'Ppv', 'CF', 'NaAl']

stherm_components = '\n'.join(therm_components) + '\n'
smass_amount = ' '.join(mass_amount)
sareotherm = '\n'.join(areotherm)
ssolution_models = '\n'.join(stx21_solution_models) + '\n'

input_data = {'name': name, 'therm data': 'stx21ver.dat', 'comput option': '', 'transform component': 'N',
              'specify comput mode': '', 'saturated components': 'N', 'independent vars': 'N',
              'therm components': stherm_components, 'dependent P/T': 'Y', 'dependent var': '2',
              'polyn degree': '4', 'polyn coeff': sareotherm, 'x-var': '1', 'min-max bar': '1 2.5e5',
              'min-max temp': '100 2500', 'mass?': 'Y', 'mass amounts': smass_amount, 'output print file': 'Y',
              'exclude pure phases': 'N', 'include solution models': 'Y', 'solution model': 'stx21_solution_model.dat',
              'select models': ssolution_models, 'calc_title': 'MARS_X'}

s = '\n'.join(input_data.values()) + '\n'

sub.run('./build', input=s, shell=True, text=True, cwd='../Perple_X')
sub.run('./vertex', input=name+'\n', shell=True, text=True, cwd='../Perple_X')