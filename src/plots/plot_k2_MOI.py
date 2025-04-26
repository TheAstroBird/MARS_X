import matplotlib.pyplot as plt
import pandas as pd

datas = pd.read_excel('../data/dynamic/integral_parameters.xlsx')
# datas = pd.read_excel('../../data/archive/PhD_2nd_year/comparing_chem_models_bf97.xlsx')

p, (sp1, sp2) = plt.subplots(1, 2, figsize=(10, 5))

for sp in (sp1, sp2):
    sp.axhline(0.3634, c='k')
    sp.axhline(0.3646, c='k')
    sp.axvline(0.166, c='k')
    sp.axvline(0.182, c='k')

# plt.scatter(datas['Love number (real)'][0], datas['inertia'][0], c='b')
# if datas.shape[-1] > 1:
#     plt.scatter(datas['Love number (real)'][1:], datas['inertia'][1:], c='k')
sc1 = sp1.scatter(datas['Love number (real)'], datas['inertia'], c=datas['Chandler period'], cmap='plasma')
sc2 = sp2.scatter(datas['Love number (real)'], datas['inertia'], c=datas['Chandler period elastic'], cmap='plasma')
sp1.set_title('Inelastic CW')
sp2.set_title('Elastic CW')
p.colorbar(sc1, ax=sp1)
p.colorbar(sc2, ax=sp2)
plt.show()