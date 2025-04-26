import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

datas = pd.read_excel('../../data/dynamic/mineral_distributions_cumulative.xlsx')
for k in datas.keys():
    for i in range(datas.shape[0]):
        if np.isnan(datas[k][i]):
            datas.loc[i, k] = 0

plt.figure()
for k in datas:
    if k != 'pressure' and k != 'temperature':
        plt.plot(datas['pressure']/1e4, datas[k], c='k')
plt.xlim((0, 25))
plt.ylim((0, 100))
plt.xlabel('P, GPa')
plt.ylabel('wt, %')
plt.savefig('../../data/png/minerals_lf97.png', dpi=600)
plt.show()