import pandas as pd
import matplotlib.pyplot as plt

DATAs_perplex = pd.read_excel('../data/dynamic/mantle_distr.xlsx')

plt.figure()
plt.plot(DATAs_perplex['pressure'], DATAs_perplex['temperature'])
plt.show()