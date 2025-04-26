import pandas as pd

s = input('\nВведите название файла с распределениями в формате .tab (без расширения):\n')
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
        datas[keys[i]] = []
    while True:
        s = f.readline().split()
        if not s:
            break
        else:
            for i in range(len(keys)):
                datas[keys[i]].append(float(s[i]))
pd_datas = pd.DataFrame(datas)
print('\nВыберите корректное название файла для сохранения:')
print('1 - mantle_distributions.xlsx (распределения параметров модели):')
print('2 - mineral_distributions.xlsx (минералогический состав мантии вдоль ареотермы):')
print('3 - mineral_distributions_cumulative.xlsx (как 2, только в кумулятивном режиме)')
while True:
    n = input()
    if n == '1':
        s = '../data/dynamic/mantle_distributions.xlsx'
        break
    elif n == '2':
        s = '../data/dynamic/mineral_distributions.xlsx'
        break
    elif n == '3':
        s = '../data/dynamic/mineral_distributions_cumulative.xlsx'
        break
    else:
        print("\nНЕИЗВЕСТНЫЙ ФОРМАТ ВВОДА! Введите еще раз")
pd_datas.to_excel(s, index=False)