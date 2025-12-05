# Данная программа строит распределения для всей планеты для моделей внутреннего строения Марса
import math as m
import numpy as np
import scipy.integrate as ing
import scipy.interpolate as inp
import scipy.optimize as opt
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------

# Главное тело программы

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

composition = input('\nВведите название используемой химической модели:\n').lower()
areoterm = input('\nВведите название используемой ареотермы:\n').lower()

# Чтение внешних данных

print('\nПолучение доступа к внешним данным...')
DATAsIN = pd.read_excel("../data/dynamic/input_param.xlsx")
den_crust_const = DATAsIN['crust density'].iloc[-1] # г/см3
den_core_const = DATAsIN['core density'].iloc[-1]
depth_crust_const = DATAsIN['crust depth'].iloc[-1] # км
sulf = DATAsIN['core sulfur'].iloc[-1] # молекулярная концентрация, максимально 1
hydr = DATAsIN['core hydro'].iloc[-1]
if (sulf + hydr) > 1:
    print('ОШИБКА! Суммарное содержание серы и водорода превышает 100%')
    exit()
ferrum = 1 - sulf - hydr
mix = ferrum*55.85 + sulf*87.92 + hydr*56.86
x_fe = ferrum*55.85/mix
x_fes = sulf*87.92/mix
x_feh = hydr*56.86/mix

DATAs_perplex = pd.read_excel('../data/dynamic/mantle_distributions_' + composition + '_' + areoterm + '.xlsx')

# Построение распределений с глубиной

print('Определение первоначальных констант...')
R_Mars = 3389.92 # км
M_Mars = 6.4185e23 # кг
Gravity_const = 6.6743e-11 # в СИ
rad_core_const = 1650 # км
den_av = 3*M_Mars/(4*m.pi*R_Mars**3) * 1e-12 # г/см3
den_crust_const /= den_av
den_core_const /= den_av
depth_crust_const /= R_Mars
rad_core_const /= R_Mars
P_const = 3*Gravity_const*M_Mars**2/(4*m.pi*R_Mars**4) * 1e-21
T_exp_const = 3e-5*Gravity_const*M_Mars/5e2/R_Mars/1e3
n_grid = 10000 # деление сетки (определяет максимальный шаг)
t_core_bound = 0.02 # точка остановки в ядре (определяет предел интегрирования, дальше плотность полагается постоянной)
grid_step = 1 / n_grid
distr_numeric = {'density': [], 'radius': [1], 'mass': [1], 'pressure': [0], 'temperature': [300]}
DATAs_perplex['density'] /= den_av*1e3
DATAs_perplex['pressure'] /= 1e4
class system15:
    def __init__(self, dens, temp):
        self.dens = dens
        self.temp = temp

    def __call__(self, t, y):
        der1 = 3 * t**2 * self.dens(t, y)
        der2 = - P_const * self.dens(t, y) * y[0]/t**2
        der3 = self.temp(t, y)
        return [der1, der2, der3]

# 1 - в коре
print('Расчет коры...')
distr_numeric['density'].append(den_crust_const)
def den_crust_func(t, y):
    return den_crust_const
def temp_crust_func(t, y):
    fun0 = inp.interp1d(DATAs_perplex['pressure'], DATAs_perplex['temperature'],
                        fill_value='extrapolate')
    fun = inp.interp1d([1-depth_crust_const,1], [float(fun0(y[1])), 300])
    return float(fun(t))

sys15_crust = system15(den_crust_func, lambda t, y: 0) # температура считается интерполяцией, а не из диф. уравнения
crust_distr = ing.RK45(sys15_crust, 1, [1, 0, 300], 1-depth_crust_const, max_step=grid_step)
while crust_distr.status == 'running':
    crust_distr.step()
    distr_numeric['radius'].append(crust_distr.t)
    distr_numeric['density'].append(den_crust_const)
    distr_numeric['mass'].append(crust_distr.y[0])
    distr_numeric['pressure'].append(crust_distr.y[1])
    distr_numeric['temperature'].append(temp_crust_func(crust_distr.t, crust_distr.y))
l_crust = len(distr_numeric['radius'])

# 2 - в мантии
print('Расчет мантии...')
def den_mantle_func(t, y):
    fun = inp.interp1d(DATAs_perplex['pressure'], DATAs_perplex['density'], fill_value='extrapolate')
    return float(fun(y[1]))
def temp_mantle_func(t, y):
    fun = inp.interp1d(DATAs_perplex['pressure'], DATAs_perplex['temperature'], fill_value='extrapolate')
    return float(fun(y[1]))

distr_numeric['radius'].append(distr_numeric['radius'][-1] - grid_step**2)
distr_numeric['density'].append(den_mantle_func(distr_numeric['radius'][-1],
                                                [distr_numeric['mass'][-1], distr_numeric['pressure'][-1] * (1+grid_step**2),
                                                    distr_numeric['temperature'][-1]]))
distr_numeric['mass'].append(distr_numeric['mass'][-1])
distr_numeric['pressure'].append(distr_numeric['pressure'][-1] * (1+grid_step**2))
distr_numeric['temperature'].append(distr_numeric['temperature'][-1])

sys15_mantle = system15(den_mantle_func, lambda t, y: 0) # температура считается интерполяцией, а не из диф. уравнения
mantle_distr = ing.RK45(sys15_mantle, distr_numeric['radius'][-1],
                        [distr_numeric['mass'][-1], distr_numeric['pressure'][-1], distr_numeric['temperature'][-1]],
                        rad_core_const, max_step=grid_step)

while mantle_distr.status == 'running':
    mantle_distr.step()
    distr_numeric['radius'].append(mantle_distr.t)
    distr_numeric['density'].append(den_mantle_func(mantle_distr.t, mantle_distr.y))
    distr_numeric['mass'].append(mantle_distr.y[0])
    distr_numeric['pressure'].append(mantle_distr.y[1])
    distr_numeric['temperature'].append(temp_mantle_func(mantle_distr.t, mantle_distr.y))

# 3 - в ядре (первый прогон)
print('Первоначальный расчет ядра...')
class Core_state:
    def __init__(self, temp, mode):
        self.temp = temp
        self.mode = mode
        # создадим таблицы параметров в виде [den_0, alpha_0, K_T_0, K'_T_0, Kdot_T_0]
        if mode == 'Fe':
            self.param = [8/den_av, 75e-6, 183, 4, 0]
        elif mode == 'FeS':
            self.param = [4.94/den_av, 68.52e-6, 54, 4, 0.02]
        elif mode == 'FeH':
            self.param = [6.7/den_av, 0, 121, 5.31, 0]

    def __call__(self, den):
        if self.mode == 'Fe' or self.mode == 'FeS':
            f = 1/2 * ((den/self.param[0])**(2/3)-1)
            p0 = 3 * f * (1+2*f)**(5/2) * self.param[2] * (1+3/2*(self.param[3]-4)*f)
            K = self.param[2] + self.param[4] * (self.temp - 800)
            p = p0 + self.param[1]*K*(self.temp-800)
            return p
        elif self.mode == 'FeH':
            f = 1 - (self.param[0]/den)**(1/3)
            p = 3 * self.param[2] * (den/self.param[0])**(2/3) * f * m.exp(3/2*(self.param[3]-1)*f)
            return p

def den_core_func(t, y):
    if CONSTCORE:
        return den_core_const
    else:
        pfe = Core_state(y[2], 'Fe')
        pfes = Core_state(y[2], 'FeS')
        pfeh = Core_state(y[2], 'FeH')
        den_fe = opt.fsolve(lambda k: pfe(k)-y[1], 8.5/den_av)
        den_fes = opt.fsolve(lambda k: pfes(k)-y[1], 6/den_av)
        den_feh = opt.fsolve(lambda k: pfeh(k)-y[1], 7/den_av)
        return 1/(x_fe/den_fe[0] + x_fes/den_fes[0] + x_feh/den_feh[0])
def temp_core_func(t, y):
    return -T_exp_const*y[0]*y[2]/t**2

distr_core = {'density': [den_core_func(distr_numeric['radius'][-1] - grid_step**2,
                                        [distr_numeric['mass'][-1], distr_numeric['pressure'][-1] * (1+grid_step**2),
                                            distr_numeric['temperature'][-1]])],
              'radius': [distr_numeric['radius'][-1] - grid_step**2],
              'mass': [distr_numeric['mass'][-1]],
              'pressure': [distr_numeric['pressure'][-1] * (1+grid_step**2)],
              'temperature': [distr_numeric['temperature'][-1]]}
sys15_core = system15(den_core_func, temp_core_func)
core_distr = ing.RK45(sys15_core, distr_core['radius'][-1],
                        [distr_core['mass'][-1], distr_core['pressure'][-1], distr_core['temperature'][-1]],
                        t_core_bound)

while core_distr.status == 'running':
    core_distr.step()
    distr_core['radius'].append(core_distr.t)
    distr_core['density'].append(den_core_func(core_distr.t, core_distr.y))
    distr_core['mass'].append(core_distr.y[0])
    distr_core['pressure'].append(core_distr.y[1])
    distr_core['temperature'].append(core_distr.y[2])

distr_core['radius'].append(0)
distr_core['density'].append(distr_core['density'][-1])
distr_core['mass'].append(distr_core['mass'][-1] - distr_core['density'][-1]*t_core_bound**3)
distr_core['pressure'].append(distr_core['pressure'][-1] + P_const*(distr_core['density'][-1]*t_core_bound)**2/2)
distr_core['temperature'].append(distr_core['temperature'][-1]*m.exp(T_exp_const*distr_core['density'][-1]*t_core_bound**2/2))

l_to_core = len(distr_numeric['radius'])
print(l_to_core, ' ', distr_core['mass'][-1])

# 4 - jump - примерный сдвиг границы исходя из остатка масс в центре планеты
print('Вычисление "jump"-а...')
r3_c = distr_numeric['radius'][-1]**3 + distr_core['mass'][-1]/(distr_core['density'][0] - distr_numeric['density'][-1])
if r3_c > 0:
    r_c = r3_c**(1/3)
else:
    r_c = -(-r3_c) ** (1 / 3)
while distr_numeric['radius'][-1] < r_c:
    for k in distr_numeric:
        distr_numeric[k].pop()

mantle_distr = ing.RK45(sys15_mantle, distr_numeric['radius'][-1],
                        [distr_numeric['mass'][-1], distr_numeric['pressure'][-1], distr_numeric['temperature'][-1]],
                        r_c, max_step=grid_step)
while mantle_distr.status == 'running':
    mantle_distr.step()
    distr_numeric['radius'].append(mantle_distr.t)
    distr_numeric['density'].append(den_mantle_func(mantle_distr.t, mantle_distr.y))
    distr_numeric['mass'].append(mantle_distr.y[0])
    distr_numeric['pressure'].append(mantle_distr.y[1])
    distr_numeric['temperature'].append(temp_mantle_func(mantle_distr.t, mantle_distr.y))

# 5 - сдвиг границы ядра, пока масса в центре не обнулится, в два шага
# шаг 1 - интегрирование на крупной сетке
# шаг 2 - интегрирование на мелкой сетке
steps_in_core = [np.inf, grid_step]
for i in range(len(steps_in_core)):
    print(f'Подбор радиуса ядра - шаг {i+1}/2...')
    distr_core = {'density': [den_core_func(distr_numeric['radius'][-1] - grid_step**2,
                                            [distr_numeric['mass'][-1], distr_numeric['pressure'][-1] * (1+grid_step**2),
                                             distr_numeric['temperature'][-1]])],
                  'radius': [distr_numeric['radius'][-1] - grid_step**2],
                  'mass': [distr_numeric['mass'][-1]],
                  'pressure': [distr_numeric['pressure'][-1] * (1+grid_step**2)],
                  'temperature': [distr_numeric['temperature'][-1]]}
    core_distr = ing.RK45(sys15_core, distr_core['radius'][-1],
                            [distr_core['mass'][-1], distr_core['pressure'][-1], distr_core['temperature'][-1]],
                            t_core_bound, max_step=steps_in_core[i])

    while core_distr.status == 'running':
        core_distr.step()
        distr_core['radius'].append(core_distr.t)
        distr_core['density'].append(den_core_func(core_distr.t, core_distr.y))
        distr_core['mass'].append(core_distr.y[0])
        distr_core['pressure'].append(core_distr.y[1])
        distr_core['temperature'].append(core_distr.y[2])

    distr_core['radius'].append(0)
    distr_core['density'].append(distr_core['density'][-1])
    distr_core['mass'].append(distr_core['mass'][-1] - distr_core['density'][-1]*t_core_bound**3)
    distr_core['pressure'].append(distr_core['pressure'][-1] + P_const*(distr_core['density'][-1]*t_core_bound)**2/2)
    distr_core['temperature'].append(distr_core['temperature'][-1]*m.exp(T_exp_const*distr_core['density'][-1]*t_core_bound**2/2))

    l_to_core = len(distr_numeric['radius'])
    print(l_to_core, ' ', distr_core['mass'][-1])

    mass_center_count = 0
    if distr_core['mass'][-1] > 0:
        mass_center_count = 1
    elif distr_core['mass'][-1] < 0:
        mass_center_count = -1

    distr_numeric_last = {}
    while mass_center_count:
        if mass_center_count > 0:
            for k in distr_numeric:
                distr_numeric_last[k] = distr_numeric[k][-1]
                distr_numeric[k].pop()
        else:
            mantle_distr.t_bound -= grid_step
            mantle_distr.status = 'running'
            while mantle_distr.status == 'running':
                mantle_distr.step()
                distr_numeric['radius'].append(mantle_distr.t)
                distr_numeric['density'].append(den_mantle_func(mantle_distr.t, mantle_distr.y))
                distr_numeric['mass'].append(mantle_distr.y[0])
                distr_numeric['pressure'].append(mantle_distr.y[1])
                distr_numeric['temperature'].append(temp_mantle_func(mantle_distr.t, mantle_distr.y))

        distr_core_prev = distr_core
        distr_core = {'density': [den_core_func(distr_numeric['radius'][-1] - grid_step**2,
                                                [distr_numeric['mass'][-1], distr_numeric['pressure'][-1] * (1+grid_step**2),
                                                 distr_numeric['temperature'][-1]])],
                      'radius': [distr_numeric['radius'][-1] - grid_step**2],
                      'mass': [distr_numeric['mass'][-1]],
                      'pressure': [distr_numeric['pressure'][-1] * (1+grid_step**2)],
                      'temperature': [distr_numeric['temperature'][-1]]}
        core_distr = ing.RK45(sys15_core, distr_core['radius'][-1],
                              [distr_core['mass'][-1], distr_core['pressure'][-1], distr_core['temperature'][-1]],
                              t_core_bound, max_step=steps_in_core[i])

        while core_distr.status == 'running':
            core_distr.step()
            distr_core['radius'].append(core_distr.t)
            distr_core['density'].append(den_core_func(core_distr.t, core_distr.y))
            distr_core['mass'].append(core_distr.y[0])
            distr_core['pressure'].append(core_distr.y[1])
            distr_core['temperature'].append(core_distr.y[2])

        distr_core['radius'].append(0)
        distr_core['density'].append(distr_core['density'][-1])
        distr_core['mass'].append(distr_core['mass'][-1] - distr_core['density'][-1] * t_core_bound ** 3)
        distr_core['pressure'].append(distr_core['pressure'][-1] +
                                      P_const * (distr_core['density'][-1] * t_core_bound) ** 2 / 2)
        distr_core['temperature'].append(distr_core['temperature'][-1] *
                                         m.exp(T_exp_const * distr_core['density'][-1] * t_core_bound ** 2 / 2))

        l_to_core = len(distr_numeric['radius'])
        print(l_to_core, ' ', distr_core['mass'][-1])

        if mass_center_count * distr_core['mass'][-1] <= 0:
            if abs(distr_core_prev['mass'][-1]) < abs(distr_core['mass'][-1]):
                distr_core = distr_core_prev
                if mass_center_count > 0:
                    for k in distr_numeric:
                        distr_numeric[k] += [distr_numeric_last[k]]
                else:
                    while distr_numeric['radius'][-1] < mantle_distr.t_bound + grid_step:
                        for k in distr_numeric:
                            distr_numeric[k].pop()
                l_to_core = len(distr_numeric['radius'])
                print(l_to_core, ' ', distr_core['mass'][-1])
            mass_center_count = 0

for i in range(len(distr_core['mass'])):
    distr_core['mass'][i] -= distr_core['mass'][-1]

# Расчет скорости волн

# 1 - в коре
print('Расчет волн в коре...')
distr_numeric['compr velocity'] = [7.0] * l_crust
distr_numeric['shear velocity'] = [4.0] * l_crust

# 2 - в мантии
print('Расчет волн в мантии...')
vel_p_mantle_func = inp.interp1d(DATAs_perplex['pressure'], DATAs_perplex['compressional velocity'],
                                 fill_value='extrapolate')
vel_s_mantle_func = inp.interp1d(DATAs_perplex['pressure'], DATAs_perplex['shear velocity'],
                                 fill_value='extrapolate')
for i in range(l_crust, l_to_core):
    distr_numeric['compr velocity'].append(vel_p_mantle_func(distr_numeric['pressure'][i]))
    distr_numeric['shear velocity'].append(vel_s_mantle_func(distr_numeric['pressure'][i]))

# 3 - в ядре
print('Расчет волн в ядре...')
distr_core['compr velocity'] = []
distr_core['shear velocity'] = []
class Core_elasticity:
    def __init__(self, mode):
        self.mode = mode
        # таблица параметров в виде [den_0, alpha_0, K_T_0, K'_T_0, Kdot_T_0, T_0]
        if mode == 'Fe':
            self.param = [7.03/den_av, 75e-6, 105, 4.5, 0.025, 2100]
        elif mode == 'FeS':
            self.param = [4.94/den_av, 68.52e-6, 54, 4, 0.02, 1100]
        elif mode == 'FeH':
            self.param = [6.7/den_av, 1, 121, 5.31, 0, 0]

    def __call__(self, T):
        if self.mode == 'FeH':
            self.param[5] = T
        den_star = self.param[0]*m.exp(-self.param[1]*(T-self.param[5]))
        expon = self.param[4] / self.param[1] / self.param[2]
        K = self.param[2] * (den_star / self.param[0]) ** expon
        K_der = self.param[3] * m.exp(self.param[1] * (T-self.param[5]))
        return den_star, K, K_der

Fe_vel = Core_elasticity('Fe')
FeS_vel = Core_elasticity('FeS')
FeH_vel = Core_elasticity('FeH')

for i in range(len(distr_core['radius'])):
    T = distr_core['temperature'][i]

    den_fe, K_fe, K_der_fe = Fe_vel(T)
    den_fes, K_fes, K_der_fes = FeS_vel(T)
    den_feh, K_feh, K_der_feh = FeH_vel(T)

    denv = x_fe*den_fe + x_fes*den_fes + x_feh*den_feh
    denr = x_fe/den_fe + x_fes/den_fes + x_feh/den_feh
    Kv = x_fe*K_fe + x_fes*K_fes + x_feh*K_feh
    Kr = x_fe/K_fe + x_fes/K_fes + x_feh/K_feh
    Kv_der = x_fe*K_der_fe + x_fes*K_der_fes + x_feh*K_der_feh
    Kr_der = x_fe/K_der_fe + x_fes/K_der_fes + x_feh/K_der_feh

    den_star = (denv + 1/denr) / 2
    K = (Kv + 1/Kr) / 2
    K_der = (Kv_der + 1/Kr_der) / 2

    L1 = K
    L2 = 5*K - 3*K*K_der
    epsil = (1 - (distr_core['density'][i]/den_star)**(2/3)) / 2

    distr_core['compr velocity'].append(m.sqrt((1-2*epsil)**(5/2) * (L1+L2*epsil) / distr_core['density'][i] / den_av))
    distr_core['shear velocity'].append(0)

for k in distr_numeric:
    distr_numeric[k] += distr_core[k]

l_planet = len(distr_numeric['radius'])
for i in range(l_planet):
    distr_numeric['radius'][i] *= R_Mars
    distr_numeric['density'][i] *= den_av
    distr_numeric['mass'][i] *= M_Mars

if GRAPH:
    print('Построение графика...')
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.25, right=0.75)

    pres = ax.twinx()
    temp = ax.twinx()
    velo = ax.twinx()

    velo.yaxis.tick_left()
    velo.yaxis.set_label_position('left')
    temp.spines.right.set_position(("axes", 1.2))
    velo.spines.left.set_position(("axes", -0.2))

    pr, = ax.plot(distr_numeric['radius'], distr_numeric['density'], "C0")
    ax.text(700, distr_numeric['density'][l_planet - l_planet//10] + 0.1, r"$\mathrm{\rho}$")
    pp, = pres.plot(distr_numeric['radius'], distr_numeric['pressure'], "C1")
    pres.text(1000, distr_numeric['pressure'][l_planet - 3*l_planet//10] , "$P$",
              horizontalalignment='right', verticalalignment='top')
    pt, = temp.plot(distr_numeric['radius'], distr_numeric['temperature'], "C2")
    temp.text(1500, distr_numeric['temperature'][l_planet - 4*l_planet//10] + 30, '$T$')
    pvp, = velo.plot(distr_numeric['radius'], distr_numeric['compr velocity'], "C3")
    velo.text(700, distr_numeric['compr velocity'][l_planet - l_planet//10] + 0.3, "$v_\mathrm{p}$")
    velo.text(3000, distr_numeric['compr velocity'][l_planet - 9*l_planet//10] - 0.1, "$v_\mathrm{p}$",
              horizontalalignment='right', verticalalignment='top')
    pvs, = velo.plot(distr_numeric['radius'], distr_numeric['shear velocity'], "C4")
    velo.text(2200, distr_numeric['shear velocity'][l_planet - 7*l_planet//10], "$v_\mathrm{s}$",
              horizontalalignment='right', verticalalignment='top')

    ax.set(xlim=(0, R_Mars), ylim=(0, 10), xlabel="Радиус, $км$", ylabel="Плотность, $г/см^3$")
    pres.set(ylim=(0, 50), ylabel="Давление, $ГПа$")
    temp.set(ylim=(0, 2500), ylabel="Температура, $К$")
    velo.set(ylim=(0, 25), ylabel="Скорость, $км/с$")

    if SAVEPIC:
        print('Сохранение графика...')
        plt.savefig(fname=s_out, dpi=600, bbox_inches='tight', pad_inches=0.1)
    plt.show()

# Запись данных в файл
print('Запись данных в файл...')

distr_numeric['gravity'] = []
for i in range(l_planet-1):
    distr_numeric['gravity'].append(Gravity_const*distr_numeric['mass'][i] / (distr_numeric['radius'][i]*1e3)**2)
distr_numeric['gravity'].append(0)

DATAs = pd.DataFrame(distr_numeric)
DATAs.to_excel("../data/dynamic/model_distributions.xlsx", index=False)

DATAsIN['core radius'] = [distr_numeric['radius'][l_to_core]]
DATAsIN.to_excel("../data/dynamic/input_param.xlsx", index=False)

print('Программа выполнена')