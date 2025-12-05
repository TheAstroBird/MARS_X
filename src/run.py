from model import Mars
from plot import plot_2d_models
from datasets import refine
import pandas as pd

comp = 'LF97'
areo = 'ATM'

t = 'archive/PhD_2nd_year/comparing_chem_models'
t_ca = t + '_' + comp.lower() + '_' + areo.lower()
# tar = 'png/comparing_chem_models_' + comp.lower() + '_' + areo.lower() + '.png'
if False:
    c = Mars(composition=comp,
             areoterm=areo,
             density_crust=3.1,  # 2.7-3.1
             depth_crust=24,  # 24-72
             sulfur_core=0.0)
    c.calculate()
    print(c.MOI)
    # c.save_integrals(target=t,
    #                  rewrite=False,
    #                  add_composition_name=True)
if False:
    for xs in [0.3]:
        c = Mars(composition=comp,
                 areoterm=areo,
                 density_crust=2.7, # 2.7-3.1
                 depth_crust=72, # 24-72, 2.7: -, 2.8: -, 2.9: -, 3.0: -
                 sulfur_core=xs)
        c.ced()
        for visc_exp in range(21, 17, -1): #21-18
            print(f'\nCalculation of model with {xs:.1f} FeS and viscosity of exponent of {visc_exp}')
            c.viscosity = 10**visc_exp
            c.cid()
            c.cip()
            c.save_integrals(target=t,
                             rewrite=False,
                             add_composition_name=True,
                             suffix='')
            if c.Love_sun.real > 0.182 or c.MOI > 0.3646:
                break

#plotting
# plot_2d_models(path=t_ca + '', color_axis=True, compare=True, save=True, target=tar)
plot_2d_models(path=t_ca + '', color_axis=False, compare=True, save=True, target='png/comparing_chem_models_lf97_atm_bw.png')

# refining
if False:
    refine(pd.read_excel('../data/' + t_ca + '.xlsx')).to_excel('../data/' + t_ca + '.xlsx', index=False)