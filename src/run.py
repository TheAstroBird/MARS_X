from model import Mars
from plot import plot_2d_models

comp = 'BF97'
areo = 'ATL'

t = 'archive/PhD_2nd_year/comparing_chem_models'
t_ca = t + '_' + comp.lower() + '_' + areo.lower()
for xs in []:
    c = Mars(composition='BF97',
             areoterm='ATL',
             density_crust=2.7, # 2.7-3.1
             depth_crust=72, # 24-72
             sulfur_core=xs)
    c.ced()
    for visc in [10**21, 10**20, 10**19, 10**18]:
        print(f'\nCalculation of model with {xs:.1f}FeS and viscosity of {visc:.0e} Pa*s')
        c.viscosity = visc
        c.cid()
        c.cip()
        c.save_integrals(target=t,
                         rewrite=False,
                         add_composition_name=True)
plot_2d_models(path=t_ca)