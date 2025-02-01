import matplotlib.pyplot as plt
import pandas as pd

distr = pd.read_excel('../../data/dynamic/mantle_distr.xlsx')
plt.figure(figsize=(5, 4), layout='tight', dpi=600)
plt.plot(distr['pressure']/1e4, distr['compressional velocity'], 'r', distr['pressure']/1e4, distr['shear velocity'], 'm')
plt.grid(True)
plt.xlabel('Давление, $ГПа$')
plt.ylabel('Скорость, $км/с$')
plt.xlim(left=0)
plt.savefig('../../data/archive/dep_seis_pres_perplex_test.jpg')