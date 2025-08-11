import numpy as np
import pandas as pd

comp = input('\nВведите название химической модели:\n')
areo = input('\nВведите название ареотермы:\n')
print('\nВыберите номер файла:')
print('1 - mantle_distributions.xlsx (распределения параметров модели):')
print('2 - mineral_distributions.xlsx (минералогический состав мантии вдоль ареотермы):')
print('3 - mineral_distributions_cumulative.xlsx (как 2, только в кумулятивном режиме)')
while True:
    n = input()
    if n in ['1', '2', '3']:
        break
    else:
        print("\nНЕИЗВЕСТНЫЙ ФОРМАТ ВВОДА! Введите еще раз")
s = 'mars_' + comp.lower() + '_' + areo.lower() + '_' + n
r_limit = {
    'ma79': 12,
    'bf97': 9.9,
    'lf97': 9.7,
    's99': 9.5,
    'kc08': 10,
    't13': 9.8
}
left_limit = {
    'atl': 10,
    'atm': 8,
    'ath': 5
}
r_i = round(r_limit[comp.lower()]*100)
l_i = round(min(r_limit[comp.lower()] - 1, left_limit[areo.lower()])*100)
with open('../Perple_X/'+s+'.tab') as f:
    for i in range(8):
        f.readline()
    keys = f.readline().split()
    datas = {}
    for i in range(len(keys)):
        if keys[i] == 'P(bar)':
            keys[i] = 'pressure'
        elif keys[i] == 'T(K)':
            keys[i] = 'temperature'
        elif keys[i] == 'rho,kg/m3':
            keys[i] = 'density'
        elif keys[i] == 'vp,km/s':
            keys[i] = 'compressional velocity'
        elif keys[i] == 'vs,km/s':
            keys[i] = 'shear velocity'
        j = 2
        while keys[i] in datas:
            keys[i] += '_' + str(j)
            j += 1
        datas[keys[i]] = []
    while True:
        s = f.readline().split()
        if not s:
            break
        else:
            for i in range(len(keys)):
                datas[keys[i]].append(float(s[i]))
if n == '1':
    print(l_i, r_i)
    print(datas['pressure'][r_i])
    for i in range(l_i):
        for k in datas:
            if k != 'pressure':
                datas[k][i] = datas[k][r_i] - (r_i - i)*(datas[k][r_i]-datas[k][l_i])/(r_i-l_i)
pd_datas = pd.DataFrame(datas)
drop = []
for k in datas:
    for k1 in datas:
        if k1 == k + '_2':
            pd_datas[k] = (pd_datas[k1].fillna(0) + pd_datas[k].fillna(0)).replace(0, np.nan)
            drop.append(k1)
            print(k, k1)
pd_datas.drop(columns=drop, inplace=True)
if n == '1':
    s = '../data/dynamic/mantle_distributions_' + comp.lower() + '_' + areo.lower() + '.xlsx'
elif n == '2':
    s = '../data/dynamic/mineral_distributions_' + comp.lower() + '_' + areo.lower() + '.xlsx'
elif n == '3':
    s = '../data/dynamic/mineral_distributions_cumulative_' + comp.lower() + '_' + areo.lower() + '.xlsx'
pd_datas.to_excel(s, index=False)