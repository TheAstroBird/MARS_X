import pandas as pd
import numpy as np

d = pd.read_excel('../data/archive/areoterms/areoterm_ATH.xlsx')

# pr = [1.032776e-24, -1.978521e-18, 1.148117e-12, -2.864355e-7, 3.305711e-2, 2.835077e2] # ATL
# pr = [1.193499e-23, -1.072304e-17, 3.646187e-12, -5.845511e-7, 4.576164e-2, 3.093017e2] # ATM
pr = [1.55991752e-33, -1.73074711e-27, 7.82149469e-22, -1.85986645e-16,
      2.50733884e-11, -1.91681834e-6, 7.93348343e-2, 3.01694513e2] # ATH
x = np.linspace(1, 249901, 2500)
y1 = np.polyval(pr, x[:1000])
y2 = np.linspace(y1[-1], 2070, 1500) # stop: ATL - 1860; ATM - 1965; ATH - 2070
y = np.concatenate((y1, y2))
dic = {'pressure': x, 'temperature': y}
newd = pd.DataFrame(dic)
newd.to_excel('../data/archive/areoterms/areoterm_dot_ATH.xlsx', index=False)