import pandas as pd
import matplotlib.pyplot as plt

d = pd.read_excel('../data/dynamic/mantle_distributions_lf97_atl.xlsx')
print(d)
plt.figure()
plt.plot(d['pressure'], d['density'], c='k')
plt.show()