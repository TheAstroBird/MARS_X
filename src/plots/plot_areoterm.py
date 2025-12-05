import pandas as pd
import matplotlib.pyplot as plt

d0 = pd.read_excel('../../data/archive/areoterms/areoterm_dot_AT0.xlsx')
dL = pd.read_excel('../../data/archive/areoterms/areoterm_dot_ATL.xlsx')
dM = pd.read_excel('../../data/archive/areoterms/areoterm_dot_ATM.xlsx')
dH = pd.read_excel('../../data/archive/areoterms/areoterm_dot_ATH.xlsx')
plt.figure()
# plt.plot(d0['pressure']/1e4, d0['temperature'])
plt.plot(dL['pressure']/1e4, dL['temperature'])
plt.plot(dM['pressure']/1e4, dM['temperature'])
plt.plot(dM['pressure']/1e4, dH['temperature'])
plt.xlim(0, 25)
plt.ylim(bottom=300)
plt.xlabel('P, GPa')
plt.ylabel('T, K')
plt.grid()
plt.show()
positive = True
d = dL
for i in range(1, d.shape[0]):
    if d['temperature'][i] < d['temperature'][i-1]:
        positive = False
        break
print('Good') if positive else print('Bad')