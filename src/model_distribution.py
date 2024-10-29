# Данная программа строит распределения для всей планеты для моделей внутреннего строения Марса
import numpy as np
import scipy.integrate as ing
import scipy.interpolate as inp
import scipy.optimize as opt
import math as m
import matplotlib.pyplot as plt
import pandas as pd

# ----------------------------------------------------------------------------------------

# Определяемые функции

# ----------------------------------------------------------------------------------------

# основное тело программы

# ----------------------------------------------------------------------------------------

# Режим работы программы

AUTO = False
if input('\nДля запуска кода нажмите ввод:\n') == 'Y':
    AUTO = True
    print('Включен режим задвоения ввода (необходим для автоматической обратки программы)\n')

CONSTCORE = False
print('Использовать приближение постоянной плотности ядра? (Y или N)')
while True:
    s = input()
    if AUTO:
        print(s)
    if s == "Y" or s == "y":
        CONSTCORE = True
        break
    elif s == "N" or s == "n":
        break
    else:
        print("\nНЕИЗВЕСТНЫЙ ФОРМАТ ВВОДА! Введите еще раз")
GRAPH = False
print('\nПостроить график полученных распределений? (Y или N)')
while True:
    s = input()
    if AUTO:
        print(s)
    if s == "Y" or s == "y":
        GRAPH = True
        break
    elif s == "N" or s == "n":
        break
    else:
        print("\nНЕИЗВЕСТНЫЙ ФОРМАТ ВВОДА! Введите еще раз")
if GRAPH:
    SAVEPIC = False
    print('\nСохранить построенный график в виде изображения? (Y или N)')
    while True:
        s = input()
        if AUTO:
            print(s)
        if s == "Y" or s == "y":
            SAVEPIC = True
            break
        elif s == "N" or s == "n":
            break
        else:
            print("\nНЕИЗВЕСТНЫЙ ФОРМАТ ВВОДА! Введите еще раз")
    if SAVEPIC:
        s_out = '../data/png/'
        s = input('\nВведите название и формат для сохранения изображения (default "test.png"):\n')
        if s == '':
            s = 'test.png'
        if AUTO:
            print(s)
        s_out += s

# внешние данные

DATAsIN = pd.read_excel("../data/dynamic/integral_param.xlsx")
den_crust_const = DATAsIN['crust density'][0] # г/см3
depth_crust_const = DATAsIN['crust depth'][0] # км
sulf = DATAsIN['core sulfur'][0] # молекулярная концентрация, максимально 1
hydr = DATAsIN['core hydro'][0]
if (sulf + hydr) > 1:
    print('ОШИБКА! Суммарное содержание серы и водорода превышает 100%')
    exit()

DATAs_perplex = pd.read_excel('../data/dynamic/mantle_distr.xlsx')

# Строим распределение с ГЛУБИНОЙ

R_Mars = 3389.92 # км
M_Mars = 6.4185e23 # кг
Gravity_const = 6.6743e-11 # в СИ
rad_core_const = 1650 # км
den_av = 3*M_Mars/(4*m.pi*R_Mars**3) * 1e-12 # г/см3
den_crust_const /= den_av
depth_crust_const /= R_Mars
rad_core_const /= R_Mars
P_const = 3*Gravity_const*M_Mars**2/(4*m.pi*R_Mars**4) * 1e-21
n_grid = 10000 # деление сетки в мантии (определяет максимальный шаг)
mantle_step = 1 / n_grid
distr_numeric = {'density': [], 'radius': [1], 'mass': [1], 'pressure': [0], 'temperature': [0]}

# 1 - кора
distr_numeric['density'].append(den_crust_const)

class system15:
    def __init__(self, dens, temp):
        self.dens = dens
        self.temp = temp

    def __call__(self, t, y):
        der1 = 3 * t**2 * self.dens(t, y)
        der2 = - P_const * self.dens(t, y) * y[0]/t**2
        der3 = self.temp(t, y)
        return [der1, der2, der3]

def den_crust_func(t, y):
    return den_crust_const
def temp_crust_func(t, y):
    return 0

sys15_crust = system15(den_crust_func, lambda t, y: 0)
crust_distr = ing.RK45(sys15_crust, 1, [1, 0, 0], 1-depth_crust_const)
while crust_distr.status == 'running':
    crust_distr.step()
    distr_numeric['radius'].append(crust_distr.t)
    distr_numeric['density'].append(den_crust_const)
    distr_numeric['mass'].append(crust_distr.y[0])
    distr_numeric['pressure'].append(crust_distr.y[1])
    distr_numeric['temperature'].append(crust_distr.y[2])
i_crust = len(distr_numeric['radius']) - 1

# 2 - мантия
def den_mantle_func(t, y):
    fun = inp.interp1d(DATAs_perplex['pressure']/1e4, DATAs_perplex['density']/(den_av*1e3), fill_value='extrapolate')
    return fun(y[1])
def temp_mantle_func(t, y):
    fun = inp.interp1d(DATAs_perplex['pressure']/1e4, DATAs_perplex['temperature'], fill_value='extrapolate')
    return fun(y[1])

sys15_mantle = system15(den_mantle_func, lambda x: 0) # температура считается интерполяцией, а не из диф. уравнения
mantle_distr = ing.RK45(sys15_mantle, distr_numeric['radius'][-1],
                        [distr_numeric['mass'][-1], distr_numeric['pressure'][-1], distr_numeric['temperature'][-1]],
                        rad_core_const, max_step=mantle_step)
while mantle_distr.status == 'running':
    mantle_distr.step()
    distr_numeric['radius'].append(mantle_distr.t)
    distr_numeric['density'].append(den_mantle_func(mantle_distr.t, mantle_distr.y))
    distr_numeric['mass'].append(mantle_distr.y[0])
    distr_numeric['pressure'].append(mantle_distr.y[1])
    distr_numeric['temperature'].append(temp_mantle_func(mantle_distr.t, mantle_distr.y))
print(distr_numeric)
quit()
# 3 - ядро (первый прогон)
# class Core_state:
#     def __init__(self, temp, mode):
#         self.temp = temp
#         self.mode = mode
#
#     def __call__(self, den):
#         if self.mode == 'FeS':
#
#         elif self.mode == 'FeH':
#             return pressure?
# def den_core_func(t, y):
#
#     return opt.fsolve(Vine, 10)
# def temp_core_func(t, y):


P_core = numpy.linspace(P_b, P_c, n_core)
den_core = numpy.zeros(n_core)
ferrum = 1 - sulf - hydr
ferrum *= 55.85
sulf *= 87.92
hydr *= 56.86
mixt = ferrum + sulf + hydr
ferrum /= mixt
sulf /= mixt
hydr /= mixt
T_grid[i] = numpy.interp(P_grid[i], P_den, T_den)
if CONSTCORE:
    den_core_ = DATAsIN['core density']/den_av
    for j in range(n_core):
        den_core[j] = den_core_
else:
    for j in range(n_core):
        f = dihotomia(1, 300, T_grid[i], P_core[j], parameter=param[21])/den_av
        s = dihotomia(1, 800, T_grid[i], P_core[j], parameter=param[22])/den_av
        h = dihotomia(2, 2100, T_grid[i], P_core[j], parameter=FEH)/den_av
        den_core[j] = 1/(ferrum/f + sulf/s + hydr/h)
den_grid[i] = numpy.interp(P_grid[i], P_core, den_core)
while i > 0:
    mass_grid[i - 1] = mass_grid[i] + RK4(0, i, den_mode='linear approx', rad=rad_grid,
                                          mass=mass_grid, P = P_grid, den_app=den_core, P_app=P_core)
    if rad_grid[i] > 0.02:
        T_grid[i - 1] = T_grid[i] + (T_grid[i] * (rad_grid[i] - rad_grid[i - 1]) * Gravity_const * mass_grid[i] *
                                     M_Mars / rad_grid[i] ** 2 / (R_Mars*1e3) * 3e-5 * 51 / 5e4)
        P_grid[i - 1] = P_grid[i] + P_const * RK4(1, i, den_mode='linear approx', rad=rad_grid,
                                          mass=mass_grid, P = P_grid, den_app=den_core, P_app=P_core)
        den_grid[i-1] = numpy.interp(P_grid[i-1], P_core, den_core)
    else:
        T_grid[i-1] = T_grid[i]
        P_grid[i-1] = P_grid[i]
        den_grid[i-1] = den_grid[i]
    i -= 1
print(i_core, ' ', mass_grid[0])
# 4 - jump - сдвиг границы исходя из примерной разности масс
i = i_core
r3_c = rad_grid[i_core]**3 + mass_grid[0]/(den_grid[i_core] - den_grid[i_core+1])
if r3_c > 0:
    r_c = r3_c**(1/3)
else:
    r_c = -(-r3_c) ** (1 / 3)
i_core = int(numpy.interp(r_c, rad_grid, range(n_grid)))
while i > i_core:
        T_grid[i] = numpy.interp(P_grid[i], P_den, T_den)
        den_grid[i] = numpy.interp(P_grid[i], P_den, den) / den_av
        mass_grid[i - 1] = mass_grid[i] + RK4(0, i, den_mode='linear approx', rad=rad_grid, mass=mass_grid,
                                              P=P_grid, den_app=den / den_av, P_app=P_den)
        P_grid[i - 1] = P_grid[i] + P_const * RK4(1, i, den_mode='linear approx', rad=rad_grid, mass=mass_grid,
                                                  P=P_grid, den_app=den / den_av, P_app=P_den)
        i -= 1
i = i_core
den_grid[i] = numpy.interp(P_grid[i], P_core, den_core)
T_grid[i] = numpy.interp(P_grid[i], P_den, T_den)
while i > 0:
    mass_grid[i - 1] = mass_grid[i] + RK4(0, i, den_mode='linear approx', rad=rad_grid,
                                          mass=mass_grid, P = P_grid, den_app=den_core, P_app=P_core)
    if rad_grid[i] > 0.02:
        T_grid[i - 1] = T_grid[i] + (T_grid[i] * (rad_grid[i] - rad_grid[i - 1]) * Gravity_const * mass_grid[i] *
                                     M_Mars / rad_grid[i] ** 2 / (R_Mars*1e3) * 3e-5 * 51 / 5e4)
        P_grid[i - 1] = P_grid[i] + P_const * RK4(1, i, den_mode='linear approx', rad=rad_grid,
                                          mass=mass_grid, P = P_grid, den_app=den_core, P_app=P_core)
        den_grid[i-1] = numpy.interp(P_grid[i-1], P_core, den_core)
    else:
        T_grid[i-1] = T_grid[i]
        P_grid[i-1] = P_grid[i]
        den_grid[i-1] = den_grid[i]
    i -= 1
print(i_core, ' ', mass_grid[0])
# 5 - сдвиг границы ядра, пока масса в центре не обнулится

mass_center_count = 0
if mass_grid[0] > 0:
    mass_center_count = 1
elif mass_grid[0] < 0:
    mass_center_count = -1
while mass_center_count:
    i_core += mass_center_count
    i = i_core
    T_grid[i] = numpy.interp(P_grid[i], P_den, T_den)
    if mass_center_count < 0:
        i += 1
        den_grid[i] = numpy.interp(P_grid[i], P_den, den)/den_av
        mass_grid[i - 1] = mass_grid[i] + RK4(0, i, den_mode='linear approx', rad=rad_grid, mass=mass_grid,
                                              P=P_grid, den_app=den/den_av, P_app=P_den)
        P_grid[i - 1] = P_grid[i] + P_const * RK4(1, i, den_mode='linear approx', rad=rad_grid,
                                                  mass=mass_grid,
                                                  P=P_grid, den_app=den/den_av, P_app=P_den)
        i -= 1
    else:
        den_grid[i] = numpy.interp(P_grid[i], P_core, den_core)
    while i > 0:
        mass_grid[i - 1] = mass_grid[i] + RK4(0, i, den_mode='linear approx', rad=rad_grid,
                                              mass=mass_grid, P=P_grid, den_app=den_core, P_app=P_core)
        if rad_grid[i] > 0.02:
            T_grid[i - 1] = T_grid[i] + (T_grid[i] * (rad_grid[i] - rad_grid[i - 1]) * Gravity_const * mass_grid[i] *
                                         M_Mars / rad_grid[i] ** 2 / (R_Mars * 1e3) * 3e-5 * 51 / 5e4)
            P_grid[i - 1] = P_grid[i] + P_const * RK4(1, i, den_mode='linear approx', rad=rad_grid,
                                                      mass=mass_grid, P=P_grid, den_app=den_core, P_app=P_core)
            den_grid[i - 1] = numpy.interp(P_grid[i - 1], P_core, den_core)
        else:
            T_grid[i - 1] = T_grid[i]
            P_grid[i - 1] = P_grid[i]
            den_grid[i - 1] = den_grid[i]
        i -= 1
    if mass_center_count*mass_grid[0] <= 0:
        mass_center_count = 0
    print(i_core, ' ', mass_grid[0])
