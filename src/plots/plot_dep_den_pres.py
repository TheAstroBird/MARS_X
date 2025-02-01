import matplotlib.pyplot as plt
import pandas as pd

distr = pd.read_excel('../../data/dynamic/mantle_distr.xlsx')
plt.figure(figsize=(5, 4), layout='tight', dpi=600)
plt.plot(distr['pressure']/1e4, distr['density']/1e3, c='b')
plt.grid(True)
plt.xlabel('Давление, $ГПа$')
plt.ylabel('Плотность, $г/см^3$')
plt.xlim(left=0)
plt.savefig('../../data/archive/dep_den_pres_perplex_test.jpg')
