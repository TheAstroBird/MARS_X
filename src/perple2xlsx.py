import pandas as pd

s = input('\nВведите название файла с распределениями в формате .tab (без расширения):\n')
DATAs_perplex = pd.read_csv('../Perple_X/'+s+'.tab', engine='python', sep='      ', header=8)
cols = []
for i in range(len(DATAs_perplex.columns)):
    if 'rho' in DATAs_perplex.columns[i]:
        cols.append('density')
    elif 'P' in DATAs_perplex.columns[i]:
        cols.append('pressure')
    elif 'T' in DATAs_perplex.columns[i]:
        cols.append('temperature')
    elif 'vp' in DATAs_perplex.columns[i]:
        cols.append('compressional velocity')
    elif 'vs' in DATAs_perplex.columns[i]:
        cols.append('shear velocity')
    else:
        cols.append(str(DATAs_perplex.columns[i]))
DATAs_perplex.to_excel('../data/dynamic/mantle_distr.xlsx', header=cols, index=False)