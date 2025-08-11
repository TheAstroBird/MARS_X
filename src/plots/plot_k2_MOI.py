import matplotlib.pyplot as plt
import pandas as pd

# datas = pd.read_excel('../../data/dynamic/integral_parameters.xlsx')
datas = pd.read_excel('../../data/archive/PhD_2nd_year/comparing_chem_models_bf97_atl.xlsx')

p, (sp1, sp2) = plt.subplots(1, 2, figsize=(14, 5))

for sp in (sp1, sp2):
    sp.axhline(0.3634, c='k')
    sp.axhline(0.3646, c='k')
    sp.axvline(0.166, c='k')
    sp.axvline(0.182, c='k')

# plt.scatter(datas['Love number (real)'][0], datas['inertia'][0], c='b')
# if datas.shape[-1] > 1:
#     plt.scatter(datas['Love number (real)'][1:], datas['inertia'][1:], c='k')
sc1 = sp1.scatter(datas['Love_number'].apply(lambda x: complex(x).real), datas['inertia'], c=datas['Chandler_period'], cmap='plasma')
sc2 = sp2.scatter(datas['Love_number'].apply(lambda x: complex(x).real), datas['inertia'], c=datas['Chandler_period_elastic'], cmap='plasma')
#sp1.set_title('Неупру CW')
#sp2.set_title('Elastic CW')
sp1.set_ylabel("I/(MR$^2$)")
sp2.set_ylabel("I/(MR$^2$)")
sp1.set_xlabel('k$_2$')
sp2.set_xlabel('k$_2$')
sp1.set_xticks((0.166, 0.168, 0.170, 0.172, 0.174, 0.176, 0.178, 0.180, 0.182))
sp2.set_xticks((0.166, 0.168, 0.170, 0.172, 0.174, 0.176, 0.178, 0.180, 0.182))
sp1.tick_params(axis='x', labelrotation=45)
sp2.tick_params(axis='x', labelrotation=45)
sp1.grid()
sp2.grid()
p.colorbar(sc1, ax=sp1)
p.colorbar(sc2, ax=sp2)
#plt.savefig('../../data/png/selection_bf97_bw', dpi=600, bbox_inches='tight')
plt.show()