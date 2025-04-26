import subprocess
import pandas as pd

param_distr = {'AUTO': 'Y', 'CONSTCORE': 'N', 'GRAPH': 'N', 'SAVEPIC': 'N', 'save_out': 'distr_test.png'}
param_integral = {'AUTO': 'Y', 'VISCOSITY': 'Y', 'model_name': 'DW85', 'MELTLAYER': 'Y', 'CREEPFUNCTION': 'N',
                  'REWRITEFILE': 'N'}
integral_param = {'crust density': [2.7], 'crust depth': [70],
                  'core density': [6.1], 'core hydro': [0], 'core sulfur': [0.8],
                  'viscosity': [1e21], 'andrade': [0.3], 'visc_melt': [1e9]}

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