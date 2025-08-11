import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams.update({'font.size': 13})

# datas = pd.read_excel('../../data/dynamic/integral_parameters.xlsx')
models = ['bf97', 'ma79', 't13', 'lf97', 's99', 'kc08']
D = []
for model in models:
    D.append(pd.read_excel('../../data/archive/PhD_2nd_year/comparing_chem_models_' + model + '.xlsx'))

p, sps = plt.subplots(3, 2, figsize=(18, 22))
t = ["a)", "b)", "c)", "d)", "e)", "f)"]
x = [0.164, 0.1635, 0.164, 0.164, 0.164, 0.164]
y = [0.3647, 0.3650, 0.3647, 0.3647, 0.3647, 0.3647]
ind = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]

for i in range(6):
    datas = D[i]
    sp = sps[ind[i]]
    sp.axhline(0.3634, c='k')
    sp.axhline(0.3646, c='k')
    sp.axvline(0.166, c='k')
    sp.axvline(0.182, c='k')

    sc = sp.scatter(datas['Love number (real)'], datas['inertia'], c=datas['Chandler period'], cmap='plasma',
                      vmin=209.4, vmax=213.1)
    sp.set_ylabel("I/(MR$^2$)")
    sp.set_xlabel('k$_2$')
    sp.set_xticks((0.166, 0.168, 0.170, 0.172, 0.174, 0.176, 0.178, 0.180, 0.182))
    sp.tick_params(axis='x', labelrotation=45)
    sp.grid()
    p.colorbar(sc, ax=sp)
    sp.text(x[i], y[i], t[i], ha="center", va="center", size=20, weight="bold")
plt.savefig('../../data/png/selection_chemmodels_en', dpi=600, bbox_inches="tight")
plt.show()