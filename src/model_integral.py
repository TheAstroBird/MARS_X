# Данная программа считает интегральные параметры для построенной модели внутреннего строения Марса
import math as m
import pandas as pd
import numpy as np
import scipy.integrate as ing
import scipy.interpolate as inp
import scipy.linalg as lin

# --------------------------------------------------------------------------------

# Главное тело программы

# --------------------------------------------------------------------------------

# Режим работы программы

AUTO = False
if input('\nДля запуска кода нажмите ввод:\n') == 'Y':
    AUTO = True
    print('Включен режим задвоения ввода (необходим для автоматической обратки программы)\n')

VISCOSITY = False
print('При расчете учитывать неупругость? (Y или N)')
while True:
    s = input()
    if AUTO:
        print(s)
    if s == "Y" or s == "y":
        VISCOSITY = True
        break
    elif s == "N" or s == "n":
        break
    else:
        print("\nНЕИЗВЕСТНЫЙ ФОРМАТ ВВОДА! Введите еще раз")
if VISCOSITY:
    DW85 = False
    print("\nВведите обозначение рассматриваемой химической модели (напр., DW85):")
    s = input()
    if AUTO:
        print(s)
    if 'DW' in s:
        DW85 = True
    MELTLAYER = False
    print('\nУчитьвать наличие расплавленного слоя? (Y или N)')
    while True:
        s = input()
        if AUTO:
            print(s)
        if s == "Y" or s == "y":
            MELTLAYER = True
            break
        elif s == "N" or s == "n":
            break
        else:
            print("\nНЕИЗВЕСТНЫЙ ФОРМАТ ВВОДА! Введите еще раз")
    CREEPFUNCTION = False
    print('\nДополнительно посчитать неупругий вклад в чандлеровской период с использованием функции крипа? (Y или N)')
    while True:
        s = input()
        if AUTO:
            print(s)
        if s == "Y" or s == "y":
            CREEPFUNCTION = True
            break
        elif s == "N" or s == "n":
            break
        else:
            print("\nНЕИЗВЕСТНЫЙ ФОРМАТ ВВОДА! Введите еще раз")
REWRITEFILE = False
print('\nПерезаписать файл integral_parameters.xlsx? (Y или N)')
while True:
    s = input()
    if AUTO:
        print(s)
    if s == "Y" or s == "y":
        REWRITEFILE = True
        break
    elif s == "N" or s == "n":
        break
    else:
        print("\nНЕИЗВЕСТНЫЙ ФОРМАТ ВВОДА! Введите еще раз")

# Чтение внешних данных
print('\nЧтение внешних данных...')
distr_numeric = pd.read_excel("../data/dynamic/model_distributions.xlsx")
param_integral = pd.read_excel("../data/dynamic/input_param.xlsx")
if VISCOSITY and DW85:
    distr_mineral = pd.read_excel("../data/dynamic/mineral_distributions.xlsx")

n = 2  # порядок определяемого числа Лява
R_Mars = 3389.92  # км
depth_crust_const = param_integral['crust depth'].iloc[-1]  # км
H_melt = 200  # км - толщина расплавленного слоя
rad_core = param_integral['core radius'].iloc[-1]
M_Mars = 6.4185e23  # кг
den_av = 3 * M_Mars / (4 * m.pi * R_Mars ** 3) * 1e-12  # г/см3
grav_av = 3.7279  # м/с^2
l_planet = len(distr_numeric['radius'])
for i in range(l_planet):
    if distr_numeric['radius'][i] <= R_Mars - depth_crust_const:
        l_crust = i + 1
        break
for i in range(l_crust, l_planet):
    if distr_numeric['radius'][i] <= rad_core:
        l_to_core = i
        break
den_crust = param_integral['crust density'].iloc[-1]
den_core = distr_numeric['density'].iloc[-1]
den_core_bound = distr_numeric['density'][l_to_core]
t_core_bound = 0.02  # точка остановки в ядре (определяет предел интегрирования, дальше плотность полагается постоянной)
n_grid = 10000  # деление сетки в мантии (определяет максимальный шаг)
grid_step = 1 / n_grid

distr_elast = {'viscosity': [], 'lame1_elast': [], 'lame2_elast': [], 'lame1': [], 'lame2': [],
               'lame1_cw': [], 'lame2_cw': [], 'Q': []}

# Введение неупругости
if VISCOSITY:
    print('Введение неупругости...')
    # 1 - Инициализация параметров
    print('Инициализация параметров неупругости...')
    andrade = param_integral['andrade'].iloc[-1]
    visc_0_exp = m.log10(param_integral['viscosity'].iloc[-1]) - 9  # вязкость коры, степень в ГПа-с
    freq_cw = 2 * m.pi / 206.9 / 86400  # частота прилива - чандлеровская частота
    freq_sun = 2 * m.pi / 44340  # частота солнечного прилива 1/с
    visc_melt_exp = m.log10(param_integral['visc_melt'].iloc[-1]) - 9

    # 2 - Расчет слоев
    print('Расчет слоев вязкости...')
    rad_trans = [R_Mars, R_Mars - depth_crust_const, R_Mars - depth_crust_const - grid_step**2]
    visc_exp = [visc_0_exp, visc_0_exp, visc_0_exp - 2]
    if DW85:
        pres_trans = []
        for i in range(1, distr_mineral.shape[0]):
            if len(pres_trans) == 0 and not m.isnan(distr_mineral['Wad'][i]) and m.isnan(distr_mineral['Wad'][i - 1]):
                pres_trans += [distr_mineral['pressure'][i] / 1e4]
            if len(pres_trans) == 1 and m.isnan(distr_mineral['O'][i]) and not m.isnan(distr_mineral['O'][i - 1]):
                pres_trans += [distr_mineral['pressure'][i] * (1+grid_step**2) / 1e4]
            if len(pres_trans) == 2 and not m.isnan(distr_mineral['Ring'][i]) and m.isnan(distr_mineral['Ring'][i - 1]):
                pres_trans += [distr_mineral['pressure'][i] / 1e4]
            if len(pres_trans) == 3 and m.isnan(distr_mineral['Wad'][i]) and not m.isnan(distr_mineral['Wad'][i - 1]):
                pres_trans += [distr_mineral['pressure'][i] * (1+grid_step**2) / 1e4]
                break
        pres_to_rad = inp.interp1d(distr_numeric['pressure'], distr_numeric['radius'])
        visc_exp += [visc_0_exp - 2, visc_0_exp - 1, visc_0_exp - 1, visc_0_exp]
        for p in pres_trans:
            r = pres_to_rad(p)
            if r < rad_core + H_melt * MELTLAYER:
                break
            else:
                rad_trans += [r]
        visc_exp = visc_exp[:len(rad_trans)]

    visc_exp += [visc_exp[-1]]
    if MELTLAYER:
        rad_trans += [rad_core + H_melt, rad_core + H_melt - grid_step**2]
        visc_exp += [visc_melt_exp] * 2
    rad_trans += [rad_core]

    def visc(r):
        fun = inp.interp1d(rad_trans, visc_exp)
        return 10 ** (fun(r))

    # 3 - расчет коэффициентов Ламе
    print('Расчет коэффициентов Ламе...')
    for i in range(l_to_core):
        lame2_elast = (distr_numeric['shear velocity'][i]**2) * distr_numeric['density'][i]
        distr_elast['lame2_elast'].append(lame2_elast)
        J = ((1 + (1j*visc(distr_numeric['radius'][i])/lame2_elast*freq_sun)**(-andrade)*m.gamma(1 + andrade))
             /lame2_elast
             - 1j/(visc(distr_numeric['radius'][i]) * freq_sun))
        lame2_inelast = 1 / J
        distr_elast['lame2'].append(lame2_inelast)
        bulk_elast = (distr_numeric['compr velocity'][i]**2)*distr_numeric['density'][i] - 4/3*distr_elast['lame2_elast'][i]
        bulk_inelast = (distr_numeric['compr velocity'][i]**2)*distr_numeric['density'][i] - 4/3*distr_elast['lame2'][i]
        distr_elast['lame1_elast'].append(bulk_elast - 2 / 3 * lame2_elast)
        distr_elast['lame1'].append(bulk_inelast - 2/3*lame2_inelast)

        J_cw = ((1 + (1j*visc(distr_numeric['radius'][i])/lame2_elast*freq_cw)**(-andrade)*m.gamma(1 + andrade))
                /lame2_elast
                - 1j/(visc(distr_numeric['radius'][i])*freq_cw))
        distr_elast['lame2_cw'].append(1 / J_cw)
        bulk_cw = ((distr_numeric['compr velocity'][i]**2)*distr_numeric['density'][i] -
                   4/3*distr_elast['lame2_cw'][i])
        distr_elast['lame1_cw'].append(bulk_cw - 2/3/J_cw)

        distr_elast['Q'].append(lame2_inelast.real / lame2_inelast.imag)
        distr_elast['viscosity'].append(visc(distr_numeric['radius'][i])*1e9)

# Расчет коэффициентов Ламе для упругой модели
else:
    print('Расчет коэффициентов Ламе...')
    for i in range(l_to_core):
        distr_elast['lame2'].append((distr_numeric['shear velocity'][i] ** 2) * distr_numeric['density'][i])
        distr_elast['lame2_elast'].append(distr_elast['lame2'][i])
        distr_elast['lame1'].append((distr_numeric['compr velocity'][i] ** 2) * distr_numeric['density'][i] -
                                    4 / 3 * distr_elast['lame2'][i])
        distr_elast['lame1_elast'].append(distr_elast['lame1'][i])

        distr_elast['lame2_cw'].append(distr_elast['lame2'][i])
        distr_elast['lame1_cw'].append(distr_elast['lame1'][i])

        distr_elast['Q'].append(m.inf)
        distr_elast['viscosity'].append(0)

for i in range(l_to_core, l_planet):
    distr_elast['lame2'].append(0)
    distr_elast['lame2_elast'].append(distr_elast['lame2'][i])
    distr_elast['lame1'].append((distr_numeric['compr velocity'][i] ** 2) * distr_numeric['density'][i])
    distr_elast['lame1_elast'].append(distr_elast['lame1'][i])

    distr_elast['lame2_cw'].append(0)
    distr_elast['lame1_cw'].append(distr_elast['lame1'][i])

    distr_elast['Q'].append(m.inf)
    distr_elast['viscosity'].append(0)

# Обезразмеривание
print('Обезразмеривание данных...')
distr_numeric['radius'] /= R_Mars
distr_numeric['density'] /= den_av
distr_numeric['gravity'] /= grav_av
for i in range(l_planet):
    for k in distr_elast:
        if k != 'Q' and k != 'viscosity':
            distr_elast[k][i] /= den_av * R_Mars * grav_av / 1e3
rad_core /= R_Mars

# Расчет модельных массы и момента инерции формулой Симпсона
print('Расчет модельного значения массы...')
integral_dict = {'mass': [], 'inertia': []}
for i in range(l_planet):
    integral_dict['mass'] += [4 * m.pi * distr_numeric['density'][i] * distr_numeric['radius'][i] ** 2]
    integral_dict['inertia'] += [8 / 3 * m.pi * distr_numeric['density'][i] * distr_numeric['radius'][i] ** 4]
M_Mars_model = - ing.trapezoid(integral_dict['mass'], x=distr_numeric['radius'])
M_Mars_model *= 3/4/m.pi
print('M =', M_Mars_model)
print('Расчет модельного значения момента инерции...')
I_Mars_model = - ing.trapezoid(integral_dict['inertia'], x=distr_numeric['radius'])
I_Mars_model *= 3/4/m.pi
print('I =', I_Mars_model)
I_core_model = - ing.trapezoid(integral_dict['inertia'][l_to_core:], x=distr_numeric['radius'][l_to_core:])
I_core_model *= 3/4/m.pi

# Вычисление числа Молоденского
print('Вычисление числа Молоденского...')
den_func = inp.interp1d(distr_numeric['radius'], distr_numeric['density'])
grav_func = inp.interp1d(distr_numeric['radius'], distr_numeric['gravity'])

def Molodensky_func(t, y):
    return - y[0] ** 2 - 2 * (n + 1) / t * y[0] - 3 * (
                den_func(t + 0.0001) - den_func(t - 0.0001)) / 0.0002 / grav_func(t)

Molodensky_solution = ing.RK45(Molodensky_func, t_core_bound, [0], distr_numeric['radius'][l_to_core] - 0.00011,
                               max_step=grid_step)
while Molodensky_solution.status == 'running':
    Molodensky_solution.step()
Molodensky_number = Molodensky_solution.y[0]
print('gamma =', Molodensky_number)

# Вычисление числа Лява k2
print('Расчет числа Лява k2...')
lame1_func = inp.interp1d(distr_numeric['radius'], distr_elast['lame1'])
lame2_func = inp.interp1d(distr_numeric['radius'], distr_elast['lame2'])
lame1_cw_func = inp.interp1d(distr_numeric['radius'], distr_elast['lame1_cw'])
lame2_cw_func = inp.interp1d(distr_numeric['radius'], distr_elast['lame2_cw'])
lame1_elst_func = inp.interp1d(distr_numeric['radius'], distr_elast['lame1_elast'])
lame2_elst_func = inp.interp1d(distr_numeric['radius'], distr_elast['lame2_elast'])

class Love_RHS:
    def __init__(self, freq_mode):
        if freq_mode == 'sun':
            self.lame1 = lame1_func
            self.lame2 = lame2_func
        elif freq_mode == 'cw':
            self.lame1 = lame1_cw_func
            self.lame2 = lame2_cw_func
        elif freq_mode == 'elast':
            self.lame1 = lame1_elst_func
            self.lame2 = lame2_elst_func

    def __call__(self, x, y):
        lame1 = self.lame1(x)
        lame2 = self.lame2(x)
        den = den_func(x)
        grav = grav_func(x)
        eq = [0] * 6
        eq[0] = -2*lame1/(lame1 + 2*lame2)*y[0]/x + y[1]/(lame1 + 2*lame2) + lame1*n*(n + 1)/(lame1 + 2*lame2)*y[4]/x
        eq[1] = ((-4*den*grav*x + 4*lame2*(3*lame1 + 2*lame2)/(lame1 + 2*lame2))*y[0]/x**2 - 4*lame2/(lame1 + 2*lame2)*y[1]/x +
                 (n*(n + 1)*den*grav*x - 2*lame2*(3*lame1 + 2*lame2)*n*(n + 1)/(lame1 + 2*lame2))*y[4]/x**2 + n*(n + 1)/x*y[5] -
                 den*y[3])
        eq[2] = 3*den*y[0] + y[3]
        eq[3] = -3*den*(n + 1)*n/x*y[4] + n*(n + 1)/x**2*y[2] - 2/x*y[3]
        eq[4] = -y[0]/x + y[4]/x + y[5]/lame2
        eq[5] = ((den*grav/x - 2*lame2*(3*lame1 + 2*lame2)/(lame1 + 2*lame2)/x**2)*y[0] -
                 lame1/(lame1 + 2*lame2)*y[1]/x +
                 2*lame2/(lame1 + 2*lame2)*(lame1*(2*n**2 + 2*n - 1) + 2*lame2*(n**2 + n - 1))*y[4]/x**2 - 3*y[5]/x -
                 den*y[2]/x)
        return eq

print('Расчет упругого числа Лява k2...')
Love_RHS_elast = Love_RHS('elast')
Love_elast_solution = [ing.RK45(Love_RHS_elast, distr_numeric['radius'][l_to_core - 1],
                              np.array([1,
                                        distr_numeric['density'][l_to_core - 1]*distr_numeric['gravity'][l_to_core - 1],
                                        0, -3*distr_numeric['density'][l_to_core - 1], 0, 0], dtype="complex_"),
                              1, max_step=grid_step),
                     ing.RK45(Love_RHS_elast, distr_numeric['radius'][l_to_core - 1],
                              np.array([0, -distr_numeric['density'][l_to_core - 1], 1,
                               n/distr_numeric['radius'][l_to_core - 1] + Molodensky_number, 0, 0], dtype="complex_"),
                              1, max_step=grid_step),
                     ing.RK45(Love_RHS_elast, distr_numeric['radius'][l_to_core - 1],
                              np.array([0, 0, 0, 0, 1, 0], dtype="complex_"), 1, max_step=grid_step)]
for sol in Love_elast_solution:
    while sol.status == 'running':
        sol.step()

bound_cond_matrix = [[Love_elast_solution[0].y[1], Love_elast_solution[1].y[1], Love_elast_solution[2].y[1]],
                     [Love_elast_solution[0].y[3] + (n+1) * Love_elast_solution[0].y[2],
                      Love_elast_solution[1].y[3] + (n + 1) * Love_elast_solution[1].y[2],
                      Love_elast_solution[2].y[3] + (n + 1) * Love_elast_solution[2].y[2]],
                     [Love_elast_solution[0].y[5], Love_elast_solution[1].y[5], Love_elast_solution[2].y[5]]]
bound_cond_RHS = [0, 2*n + 1, 0]
solutions_coeffs = lin.solve(bound_cond_matrix, bound_cond_RHS)
Love_elast = (solutions_coeffs[0] * Love_elast_solution[0].y[2] + solutions_coeffs[1] * Love_elast_solution[1].y[2] +
            solutions_coeffs[2] * Love_elast_solution[2].y[2] - 1)
print('k2_e =', Love_elast)

print('Расчет числа Лява k2 для чандлеровского колебания...')
Love_RHS_CW = Love_RHS('cw')
Love_CW_solution = [ing.RK45(Love_RHS_CW, distr_numeric['radius'][l_to_core - 1],
                             np.array([1,
                                       distr_numeric['density'][l_to_core - 1]*distr_numeric['gravity'][l_to_core - 1],
                                       0, -3*distr_numeric['density'][l_to_core - 1], 0, 0], dtype="complex_"),
                             1, max_step=grid_step),
                    ing.RK45(Love_RHS_CW, distr_numeric['radius'][l_to_core - 1],
                             np.array([0, -distr_numeric['density'][l_to_core - 1], 1,
                                       n/distr_numeric['radius'][l_to_core - 1] + Molodensky_number, 0, 0],
                                      dtype="complex_"), 1, max_step=grid_step),
                    ing.RK45(Love_RHS_CW, distr_numeric['radius'][l_to_core - 1], np.array([0, 0, 0, 0, 1, 0],
                                                                                           dtype="complex_"),
                             1, max_step=grid_step)]
for sol in Love_CW_solution:
    while sol.status == 'running':
        sol.step()

bound_cond_matrix = [[Love_CW_solution[0].y[1], Love_CW_solution[1].y[1], Love_CW_solution[2].y[1]],
                     [Love_CW_solution[0].y[3] + (n+1) * Love_CW_solution[0].y[2],
                      Love_CW_solution[1].y[3] + (n+1) * Love_CW_solution[1].y[2],
                      Love_CW_solution[2].y[3] + (n+1) * Love_CW_solution[2].y[2]],
                     [Love_CW_solution[0].y[5], Love_CW_solution[1].y[5], Love_CW_solution[2].y[5]]]
bound_cond_RHS = [0, 2*n + 1, 0]
solutions_coeffs = lin.solve(bound_cond_matrix, bound_cond_RHS)
Love_CW = (solutions_coeffs[0]*Love_CW_solution[0].y[2] + solutions_coeffs[1]*Love_CW_solution[1].y[2] +
           solutions_coeffs[2]*Love_CW_solution[2].y[2] - 1)
print('k2_cw =', Love_CW)

print('Расчет числа Лява k2 для Солнечного прилива...')
Love_RHS_Sun = Love_RHS('sun')
Love_Sun_solution = [ing.RK45(Love_RHS_Sun, distr_numeric['radius'][l_to_core - 1],
                              np.array([1,
                                        distr_numeric['density'][l_to_core - 1]*distr_numeric['gravity'][l_to_core - 1],
                                        0, -3*distr_numeric['density'][l_to_core - 1], 0, 0], dtype="complex_"),
                              1, max_step=grid_step),
                     ing.RK45(Love_RHS_Sun, distr_numeric['radius'][l_to_core - 1],
                              np.array([0, -distr_numeric['density'][l_to_core - 1], 1,
                               n/distr_numeric['radius'][l_to_core - 1] + Molodensky_number, 0, 0], dtype="complex_"),
                              1, max_step=grid_step),
                     ing.RK45(Love_RHS_Sun, distr_numeric['radius'][l_to_core - 1],
                              np.array([0, 0, 0, 0, 1, 0], dtype="complex_"), 1, max_step=grid_step)]
for sol in Love_Sun_solution:
    while sol.status == 'running':
        sol.step()

bound_cond_matrix = [[Love_Sun_solution[0].y[1], Love_Sun_solution[1].y[1], Love_Sun_solution[2].y[1]],
                     [Love_Sun_solution[0].y[3] + (n+1) * Love_Sun_solution[0].y[2],
                      Love_Sun_solution[1].y[3] + (n + 1) * Love_Sun_solution[1].y[2],
                      Love_Sun_solution[2].y[3] + (n + 1) * Love_Sun_solution[2].y[2]],
                     [Love_Sun_solution[0].y[5], Love_Sun_solution[1].y[5], Love_Sun_solution[2].y[5]]]
bound_cond_RHS = [0, 2*n + 1, 0]
solutions_coeffs = lin.solve(bound_cond_matrix, bound_cond_RHS)
Love_Sun = (solutions_coeffs[0] * Love_Sun_solution[0].y[2] + solutions_coeffs[1] * Love_Sun_solution[1].y[2] +
            solutions_coeffs[2] * Love_Sun_solution[2].y[2] - 1)
print('k2_sun =', Love_Sun)

# Восстановление размерности
print('Восстановление размерности...')
distr_numeric['radius'] *= R_Mars
distr_numeric['density'] *= den_av
distr_numeric['gravity'] *= grav_av
for i in range(l_planet):
    for k in distr_elast:
        if k != 'Q' and k != 'viscosity':
            distr_elast[k][i] *= den_av * R_Mars * grav_av / 1e3

# Расчет чандлеровского периода
print('Расчет чандлеровского периода...')
A = 0.362976
B = 0.363229
C = 0.365067
period_rot = 24.6229 * 3600
grav_const = 6.6743e-11

A_red = (C - B) / A
B_red = (C - A) / B
A_aver = (A + B) / 2
freq_rot = 2 * m.pi / period_rot
period_Euler = period_rot / m.sqrt(A_red * B_red)
Love_sec = 3 * grav_const * (C - A_aver) * M_Mars / freq_rot ** 2 / (R_Mars * 1e3) ** 3
print('R_core =', rad_core * R_Mars)

period_CW = period_Euler * (1 - I_core_model / m.sqrt(A * B)) / (1 - Love_CW.real / Love_sec) / 86400
print('Tw0 = ', period_CW)

period_CW_e = period_Euler * (1 - I_core_model / m.sqrt(A * B)) / (1 - Love_elast.real / Love_sec) / 86400
print('Tw0e = ', period_CW_e)

# Расчет чандлеровского периода через функцию крипа
if VISCOSITY and CREEPFUNCTION:
    print('Расчет чандлеровского периода с использованием функции крипа...')
    period_CW_creep = period_Euler * (1 - I_core_model / m.sqrt(A * B)) / (1 - Love_Sun.real / Love_sec) / 86400
    n = 0.4
    dT = period_CW_creep / 1118 * ((freq_sun / freq_cw) ** n - 1) / (Love_sec - Love_CW.real) / m.tan(n * m.pi / 2)
    print(f'Tw_cr = {period_CW_creep} + {dT}')

# Запись данных в файлы
print('Запись данных в файлы...')
param_integral['core density'] = [den_core] # in the center of the core
param_integral['inertia'] = [I_Mars_model]
param_integral['Love number (real)'] = [Love_Sun.real]
param_integral['Love number (imaginary)'] = [Love_Sun.imag]
param_integral['Chandler period'] = [period_CW]
param_integral['Chandler period elastic'] = [period_CW_e]
if not REWRITEFILE:
    param_integral_old = pd.read_excel("../data/dynamic/integral_parameters.xlsx")
    param_integral = pd.concat([param_integral_old, param_integral], ignore_index=True)
param_integral.to_excel("../data/dynamic/integral_parameters.xlsx", index=False)

for k in distr_elast:
    distr_numeric[k] = distr_elast[k]
distr_numeric.to_excel("../data/dynamic/model_distributions.xlsx", index=False)
print('Программа выполнена')
