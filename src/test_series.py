import subprocess
import pandas as pd

param_distr = {'AUTO': 'Y', 'CONSTCORE': 'N', 'GRAPH': 'Y', 'SAVEPIC': 'Y', 'save_out': 'smth.png'}
integral_param = {'crust density': [2.9], 'crust depth': [50],
                  'core density': [6.1], 'core hydro': [0], 'core sulfur': [0.8],
                  'viscosity': [12], 'andrade': [0.3], 'visc_melt': [0]}

if param_distr['GRAPH'] == 'N' or param_distr['SAVEPIC'] == 'N':
    param_distr.pop('save_out')
    if param_distr['GRAPH'] == 'N':
        param_distr.pop('SAVEPIC')
s = ''
for k in param_distr.keys():
    s += param_distr[k] + '\n'

DATAs = pd.DataFrame(integral_param)
DATAs.to_excel("../data/dynamic/integral_param.xlsx", index=False)

subprocess.run('python3 model_distribution.py', input=s, shell=True, text=True)