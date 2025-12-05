import subprocess
import pandas as pd


composition = 'BF97'
areoterm = 'ATH'
core_s = 0.1
crust_den = 2.7
crust_depth = 72
rewr = 'N'
visc = 1e18
# Не забудь удалить строчку из model_integral.py
for core_s in [0.7, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]:
    for visc in [1e21, 1e20]:#, 1e19, 1e18]:
        param_distr = {'AUTO': 'Y', 'CONSTCORE': 'N', 'GRAPH': 'N', 'SAVEPIC': 'N', 'save_out': 'distr_test.png',
                       'composition': composition, 'areoterm': areoterm}
        param_integral = {'AUTO': 'Y', 'VISCOSITY': 'Y', 'MELTLAYER': 'Y', 'CREEPFUNCTION': 'N', 'REWRITEFILE': rewr,
                          'composition': composition, 'areoterm': areoterm}
        integral_param = {'crust density': [crust_den], 'crust depth': [crust_depth],
                          'core density': [6.1], 'core hydro': [0], 'core sulfur': [core_s],
                          'viscosity': [visc], 'andrade': [0.3], 'visc_melt': [1e11]}

        if param_distr['GRAPH'] == 'N' or param_distr['SAVEPIC'] == 'N':
            param_distr.pop('save_out')
            if param_distr['GRAPH'] == 'N':
                param_distr.pop('SAVEPIC')
        if param_integral['VISCOSITY'] == 'N':
            param_integral.pop('model_name')
            param_integral.pop('MELTLAYER')
            param_integral.pop('CREEPFUNCTION')

        s = '\n'.join(param_distr.values()) + '\n'

        DATAs = pd.DataFrame(integral_param)
        DATAs.to_excel("../data/dynamic/input_param.xlsx", index=False)

        subprocess.run('python3 model_distribution.py', input=s, shell=True, text=True)

        s = '\n'.join(param_integral.values()) + '\n'

        subprocess.run('python3 model_integral.py', input=s, shell=True, text=True)

subprocess.run('python3 plots/plot_k2_MOI.py', shell=True)