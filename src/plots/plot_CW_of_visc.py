import pandas as pd
import matplotlib.pyplot as plt

s = '../../data/archive/PhD_2nd_year/comparing_chem_models_'
df = pd.read_excel(s+'ma79_ath_2.xlsx')
for ds in ['bf97_atm_2.xlsx', 'bf97_ath.xlsx', 't13_atm.xlsx', 't13_ath.xlsx']:
    df = pd.concat([df, pd.read_excel(s+ds)], ignore_index=True)
for areo in ['atl', 'atm', 'ath']:
    for comp in ['lf97', 's99', 'kc08']:
        df = pd.concat([df, pd.read_excel(s+comp+'_'+areo+'.xlsx')], ignore_index=True)

plt.figure()
plt.scatter(df['viscosity'], df['Chandler_period'], c='k')
plt.xlabel('$\eta_0$, Pa$\cdot$s')
plt.ylabel('T$_W$, day')
plt.ylim(bottom=180)
plt.axhline(y=206.9, c='k')
plt.grid()
plt.xscale('log')
plt.show()
