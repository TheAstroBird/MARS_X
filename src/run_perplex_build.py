import subprocess as sub
import pandas as pd
import numpy as np

# composition = 'T13'
# areo = 'ATL'
for composition in ['BF97', 'MA79', 'LF97', 'S99', 'KC08', 'T13']:
    for areo in ['ATL', 'ATM', 'ATH']:
        mass_amount_data = pd.read_excel('../../chemical_models.xlsx', index_col=0)
        therm_components = ['NA2O', 'MGO', 'AL2O3', 'SIO2', 'CAO', 'FEO']
        if composition in ['BF97', 'MA79', 'LF97', 'S99', 'KC08', 'T13']:
            name = 'mars_' + composition.lower() + '_' + areo.lower()
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
        areotherm = {'AT0': [1.240095e3, 8.085953e-3, -1.774443e-9, -1.746166e-13, 4.278426e-19],
                     'ATL': [2.871026e2, 3.238349e-2, -2.671839e-7, 9.589441e-13, -1.228137e-18],
                     'ATM': [3.508448e2, 3.797712e-2, -3.620745e-7, 1.460068e-12, -2.051437e-18],
                     'ATH': [4.145869e2, 4.357076e-2, -4.569651e-7, 1.961191e-12, -2.874737e-18]}
        areotherm = [str(coeff) for coeff in areotherm[areo]]
        stx21_solution_models = ['C2/c', 'Wus', 'Pv', 'Pl', 'Sp', 'O', 'Wad', 'Ring', 'Opx', 'Cpx',
                                 'Aki', 'Gt', 'Ppv', 'CF', 'NaAl']

        stherm_components = '\n'.join(therm_components) + '\n'
        smass_amount = ' '.join(mass_amount)
        sareotherm = '\n'.join(areotherm)
        ssolution_models = '\n'.join(stx21_solution_models) + '\n'

        input_data = {'name': name, 'therm data': 'stx21ver.dat', 'comput option': '', 'transform component': 'N',
                      'specify comput mode': '', 'saturated components': 'N', 'independent vars': 'N',
                      'therm components': stherm_components, 'dependent P/T': 'N', 'x-var': '1', 'min-max bar': '1 2.5e5',
                      'min-max temp': '100 2500', 'mass?': 'Y', 'mass amounts': smass_amount, 'output print file': 'Y',
                      'exclude pure phases': 'N', 'include solution models': 'Y', 'solution model': 'stx21_solution_model.dat',
                      'select models': ssolution_models, 'calc_title': 'MARS_X'}

        s = '\n'.join(input_data.values()) + '\n'

        areoterm_name = 'mars_areoterm_dot_' + areo + '.dat'
        werami_data = {'name': name,
                       'comp_mode1': '4', 'path_descr1': '2', 'areoterm_file1': areoterm_name, 'every nth1': '1',
                       'property1': '2', 'individual properties1': 'N',
                       'property2': '13', 'individual properties2': 'N',
                       'property3': '14', 'individual properties3': 'N', 'finish': '0',
                       'comp_mode2': '4', 'path_descr2': '2', 'areoterm_file2': areoterm_name, 'every nth2': '1',
                       'property4': '25', 'cumulative1': 'N',
                       'comp_mode3': '4', 'path_descr3': '2', 'areoterm_file3': areoterm_name, 'every nth3': '1',
                       'property5': '25', 'cumulative2': 'Y',
                       'exit': '0'}

        swerami = '\n'.join(werami_data.values()) + '\n'

        sub.run('./build', input=s, shell=True, text=True, cwd='../Perple_X')
        sub.run('./vertex', input=name+'\n', shell=True, text=True, cwd='../Perple_X')
        sub.run('./werami', input=swerami, shell=True, text=True, cwd='../Perple_X')