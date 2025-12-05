import pandas as pd
import matplotlib.pyplot as plt

d = pd.read_excel('../data/dynamic/mantle_distributions_ma79_ath.xlsx')
print(d)
plt.figure()
plt.plot(d['pressure'], d['temperature'], c='k')
plt.show()