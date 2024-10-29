import numpy
import math
import matplotlib.pyplot as plt
from pylatex import Math, NoEscape
import os
import pandas as pd
import scipy as sp

def dihotomia(i0, T_norm, T, P, parameter=[]):
    den_1 = 0
    den_2 = 20
    den_0 = (den_1 + den_2) / 2
    P_0 = 0
    while (den_2 - den_1) > .0005:
        P_0 = formulae(i0, den_0, parameter[0], T, T_norm, parameter[1], parameter[2], parameter[3], parameter[4],
                       parameter[5], parameter[6])
        if P_0 < P:
            den_1 = den_0
        else:
            den_2 = den_0
        den_0 = (den_1 + den_2) / 2
    return den_0

def RHS(i0, rad, den, mass):
    if i0 == 0:
        return 3*rad**2*den
    elif i0 == 1:
        return -den*mass/rad**2

def linear_approx(i0, x0, x=[], y=[]):
    i = i0
    if i == 0:
        i = 1
    elif i == len(x):
        i -= 1
    elif x0 > x[i]:
        while i < len(x) and x0 > x[i]:
            i += 1
    return y[i] - (y[i]-y[i-1])*(x[i]-x0)/(x[i]-x[i-1])

def RK4(i0, i_current, T_0 = 0, T = 0, den_0 = 0, den_mode=[], rad=[], mass=[], den=[0], P=[0], den_app=[0], P_app=[0], parameter=[0]):
    k1 = numpy.zeros(2)
    if den_mode == 'given':
        den_current = den[i_current]
    elif den_mode == 'linear approx':
        den_current = numpy.interp(P[i_current], P_app, den_app,
                                   left=den_app[0] - (P_app[0]-P[i_current])/(P_app[1]-P_app[0])*(den_app[1]-den_app[0]))
    elif den_mode == 'function':
        den_current = den[i_current]
    for i in range(2):
        k1[i] = RHS(i, rad[i_current], den_current, mass[i_current])

    k2 = numpy.zeros(2)
    if den_mode == 'given':
        den_current = (den[i_current]+den[i_current-1])/2
    elif den_mode == 'linear approx':
        den_current = numpy.interp(P[i_current]-k1[1]/2*(rad[i_current]-rad[i_current-1]), P_app, den_app,
                                   left=den_app[0] - (P_app[0]-(P[i_current]-k1[1]/2*(rad[i_current]-rad[i_current-1])))/(P_app[1]-P_app[0])*(den_app[1]-den_app[0]))
    elif den_mode == 'function':
        den_current = dihotomia(2, T_0, T, P[i_current]-k1[1]/2*(rad[i_current]-rad[i_current-1]), parameter)/den_0
    for i in range(2):
        k2[i] = RHS(i, (rad[i_current]+rad[i_current-1])/2, den_current,
                    mass[i_current]-k1[0]/2*(rad[i_current]-rad[i_current-1]))

    k3 = numpy.zeros(2)
    if den_mode == 'given':
        den_current = (den[i_current]+den[i_current-1])/2
    elif den_mode == 'linear approx':
        den_current = numpy.interp(P[i_current]-k2[1]/2*(rad[i_current]-rad[i_current-1]), P_app, den_app,
                                   left=den_app[0] - (P_app[0]-(P[i_current]-k2[1]/2*(rad[i_current]-rad[i_current-1])))/(P_app[1]-P_app[0])*(den_app[1]-den_app[0]))
    elif den_mode == 'function':
        den_current = dihotomia(2, T_0, T, P[i_current] - k2[1] / 2 * (rad[i_current] - rad[i_current - 1]), parameter)/den_0
    for i in range(2):
        k3[i] = RHS(i, (rad[i_current]+rad[i_current-1])/2, den_current,
                    mass[i_current]-k2[0]/2*(rad[i_current]-rad[i_current-1]))

    k4 = numpy.zeros(2)
    if den_mode == 'given':
        den_current = den[i_current-1]
    elif den_mode == 'linear approx':
        den_current = numpy.interp(P[i_current] - k3[1]* (rad[i_current] - rad[i_current-1]), P_app, den_app,
                                   left=den_app[0] - (P_app[0]-(P[i_current]-k3[1]* (rad[i_current]-rad[i_current-1])))/(P_app[1]-P_app[0])*(den_app[1]-den_app[0]))
    elif den_mode == 'function':
        den_current = dihotomia(2, T_0, T, P[i_current] - k3[1] * (rad[i_current] - rad[i_current - 1]), parameter)/den_0
    k4[0] = RHS(0, rad[i_current - 1], den_current,
                mass[i_current] - k3[0] * (rad[i_current] - rad[i_current - 1]))
    if rad[i_current-1] == 0:
        k4[1] = 0
    else:
        k4[1] = RHS(1, rad[i_current - 1], den_current,
                    mass[i_current] - k3[0]* (rad[i_current] - rad[i_current-1]))
    return -(rad[i_current] - rad[i_current-1])/6*(k1[i0]+2*k2[i0]+2*k3[i0]+k4[i0])

# ----------------------------------------------------------------------------------------

# основное тело программы

# ----------------------------------------------------------------------------------------

# Режим работы программы

CONSTCORE = True
GRAPH = False
SAVEPIC = False
if SAVEPIC:
    s_out = '../data/png/'
    inp = input('\nВведите название и формат для сохранения изображения (default "test.png"):\n')
    if inp == '':
        inp = 'test.png'
    print(inp)
    s_out += inp
    quit()

# внешние изменяемые данные

DATAsIN = pd.read_excel("datas/Input_data.xlsx")
den_crust_const = DATAsIN['crust density'] # г/см3
depth_crust_const = 50 # км
fe_s_float = DATAsIN['Fe_s'][0]/100
sulf = DATAsIN['sulfur'][0] # молекулярная концентрация, максимально 1
hydr = DATAsIN['hydro'][0]
if (sulf + hydr) > 1:
    print('ERROR')
    exit()

# ввод данных

mineralogy = numpy.zeros((16, 12))
mineralogy[0] = [0.58, 0, 0, 0.2, 0.062, 0.058, 0.088, 0.012, 0, 0, 0, 0]  # 2 GPa
mineralogy[1] = [0.58, 0, 0, 0.2, 0.076, 0.044, 0.088, 0.012, 0, 0, 0, 0]  # 3 GPa
mineralogy[2] = [0.58, 0, 0, 0.19, 0.073, 0.057, 0.088, 0.012, 0, 0, 0, 0]  # 4 GPa
mineralogy[3] = [0.58, 0, 0, 0.18, 0.072, 0.068, 0.088, 0.012, 0, 0, 0, 0]  # 5 GPa
mineralogy[4] = [0.58, 0, 0, 0.12, 0.123, 0.077, 0.088, 0.012, 0, 0, 0, 0]  # 7 GPa
mineralogy[5] = [0.58, 0, 0, 0.08, 0.156, 0.084, 0.088, 0.012, 0, 0, 0, 0]  # 8
mineralogy[6] = [0.58, 0, 0, 0.03, 0.200, 0.090, 0.088, 0.012, 0, 0, 0, 0]  # 9
mineralogy[7] = [0.58, 0, 0, 0, 0.23, 0.09, 0.088, 0.012, 0, 0, 0, 0]  # 10
mineralogy[8] = [0.58, 0, 0, 0, 0.22, 0.10, 0.088, 0.012, 0, 0, 0, 0]  # 11
mineralogy[9] = [0.58, 0, 0, 0, 0.185, 0.115, 0.088, 0.012, 0, 0, 0, 0]  # 13
mineralogy[10] = [0, 0.51, 0.16, 0, 0.08, 0.10, 0.08, 0.023, 0.047, 0, 0, 0]  # 14
mineralogy[11] = [0, 0.42, 0.16, 0, 0.06, 0.11, 0.16, 0.038, 0.052, 0, 0, 0]  # 15
mineralogy[12] = [0, 0, 0.60, 0, 0, 0, 0.08, 0.06, 0.26, 0, 0, 0]  # 19
mineralogy[13] = [0, 0, 0.58, 0, 0, 0, 0.09, 0.05, 0.28, 0, 0, 0]  # 21
mineralogy[14] = [0, 0, 0, 0, 0, 0, 0.021, 0.019, 0.06, 0.743, 0.057, 0.1]  # 22.5
mineralogy[15] = [0, 0, 0, 0, 0, 0, 0.021, 0.019, 0.06, 0.743, 0.057, 0.1]  # 23.5
P = [2, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 19, 21, 22.5, 23.5]
T_p = [1400, 1475, 1550, 1625, 1750, 1800, 1850, 1900, 1940, 2000, 2025, 2040, 2070, 2075, 2079, 2080]

P_Fe_s = numpy.zeros((5, 9))
P_Fe_s[0] = [16, 18, 20, 22, 24, 25, 26, 28, 30] # Fe# numbers
P_Fe_s[1] = [14, 13.8, 13.6, 13.36, 13.12, 13, 12.9, 12.7, 12.5] # P_1 - alpha -> alpha+beta(alpha+gamma)
P_Fe_s[2] = [14.45, 14.35, 14.25, 14.15, 14.05, 14, 13.94, 13.82, 13.7] # P_2 - alpha+beta(alpha+gamma) -> beta
P_Fe_s[3] = [17.85, 17.55, 17.25, 16.95, 16.65, 16.5, 16.3, 15.9, 15.5] # P_3 - beta -> (beta+gamma)
P_Fe_s[4] = [18.35, 18.05, 17.75, 17.45, 17.15, 17, 16.8, 16.4, 16] # P_4 - (beta+gamma) -> gamma
P_1 = numpy.interp(fe_s_float*100, P_Fe_s[0], P_Fe_s[1])
P_2 = numpy.interp(fe_s_float*100, P_Fe_s[0], P_Fe_s[2])
P_3 = numpy.interp(fe_s_float*100, P_Fe_s[0], P_Fe_s[3])
P_4 = numpy.interp(fe_s_float*100, P_Fe_s[0], P_Fe_s[4])
T_1 = numpy.interp(P_1, P, T_p)
T_2 = numpy.interp(P_2, P, T_p)
T_3 = numpy.interp(P_3, P, T_p)
T_4 = numpy.interp(P_4, P, T_p)

# подсчет зависимости плотности от давления в мантии (блок 1)

mineralogy_den = numpy.zeros((11, 12))  # для подсчетов плотности в мантии
mineralogy_den[0] = [0.58, 0, 0, 0.20, 0.062, 0.058, 0.088, 0.012, 0, 0, 0, 0]  # ~2
mineralogy_den[1] = [0.58, 0, 0, 0, 0.185, 0.115, 0.089, 0.013, 0.018, 0, 0, 0]  # ~13
mineralogy_den[2] = [0, 0.51, 0.16, 0, 0.08, 0.10, 0.08, 0.023, 0.047, 0, 0, 0]  # ~14
mineralogy_den[3] = [0, 0.42, 0.16, 0, 0.06, 0.11, 0.16, 0.038, 0.052, 0, 0, 0]  # ~15
mineralogy_den[4] = mineralogy_den[3]
mineralogy_den[6] = [0, 0, 0.60, 0, 0, 0, 0.08, 0.06, 0.26, 0, 0, 0]  # ~19
mineralogy_den[5] = mineralogy_den[6]
mineralogy_den[7] = [0, 0, 0.58, 0, 0, 0, 0.09, 0.05, 0.28, 0, 0, 0]  # ~21-0
mineralogy_den[8] = mineralogy_den[7]
mineralogy_den[9] = [0, 0, 0, 0, 0, 0, 0.021, 0.019, 0.06, 0.743, 0.057, 0.1]  # ~22.5+0
mineralogy_den[10] = mineralogy_den[9]  # ~23.5

given = [0, 1, 2, 3, 6, 7, 9, 10]

dfe = 0.25-fe_s_float
fe_s = numpy.zeros((11, 12))
fe_s[0] = [0.242, 0, 0, .22, .247, .247, .341, 0, 0, 0, 0, 0]
fe_s[1] = [0.215, 0, 0, 0, 0.195, 0.195, 0.62, 0, 0, 0, 0, 0]
fe_s[2] = [0, 0.243, .33, 0, .165, .165, 0.54, 0, 0, 0, 0, .543]
fe_s[3] = [0, 0.249, .33, 0, 0.16, 0.16, 0.54, 0, 0, 0, 0, .543]
fe_s[4] = fe_s[3]
fe_s[6] = [0, 0, 0.281, 0, 0, 0, 1.0, 0, 0, 0, 0, 0.5]
fe_s[5] = fe_s[6]
fe_s[7] = [0, 0, 0.266, 0, 0, 0, 1.0, 0, 0, 0, 0, 0.5]
fe_s[8] = fe_s[7]
fe_s[9] = [0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0.16, 0, 0.493]
fe_s[10] = fe_s[9]
if dfe != 0:
    for i in range(11):
        for j in range(12):
            if fe_s[i][j] > 0:
                fe_s[i][j] -= dfe
den_miner = numpy.zeros((11, 21))
T_den = [1400, T_1, T_2, 2040, T_3, T_4, 2070, 2075, 2079, 2079, 2080]
P_den = [2, P_1, P_2, 15, P_3, P_4, 19, 21, 22.5, 22.5, 23.5]
T_0 = 300
param = numpy.zeros((23, 7))
# Olivine Mg and Fe
param[0] = [3.228, 0.3034, 0.7422, -0.5381, 129.0, 5.37, -.0224]
param[1] = [4.402, 0.2386, 1.1530, -0.0518, 137.9, 4.00, -.0258]
# β-Phase Mg and Fe
param[2] = [3.473, 0.2893, 0.5772, -0.8903, 174.0, 4.00, -.0323]
param[3] = [4.723, 0.2319, 0.7117, -0.2430, 166.0, 4.00, -.0215]
# Spinel Mg and Fe
param[4] = [3.564, 0.2497, 0.3639, -0.6531, 183.0, 4.30, -.0348]
param[5] = [4.849, 0.2697, 0.0000, -0.0000, 197.0, 4.00, -.0375]
# Pyroxene Mg-opx, Fe-opx, Mg-cpx, Fe-cpx, Ca-Mg-cpx, Ca-Fe-cpx
param[6] = [3.204, 0.2947, 0.2694, -0.5588, 107.0, 4.20, -.0200]
param[7] = [4.002, 0.3930, 0.0000, -0.0000, 101.0, 4.20, -.0200]
param[8] = [3.188, 0.2947, 0.2694, -0.5588, 107.0, 4.20, -.0200]
param[9] = [4.005, 0.3930, 0.0000, -0.0000, 101.0, 4.20, -.0200]
param[10] = [3.280, .3330, 0.0000, -0.0000, 113.0, 4.80, -.0200]
param[11] = [3.656, .2980, 0.0000, -0.0000, 119.0, 4.00, -.0200]
# Garnet-majorite Mg, Fe, Ca, Mj
param[12] = [3.566, .2311, 0.5956, -0.4538, 179.0, 4.00, -.0220]
param[13] = [4.312, .1776, 1.2140, -0.5071, 175.0, 4.00, -.0220]
param[14] = [3.600, .1951, 0.8089, -0.4972, 168.0, 6.20, -.0220]
param[15] = [3.518, .2311, 0.5956, -0.4538, 161.0, 4.00, -.0220]
# Perovskite Mg, Fe, Ca
param[16] = [4.107, .1982, 0.8180, -0.4740, 262.0, 4.00, -.0550]
param[17] = [5.154, .1982, 0.8180, -0.4740, 287.2, 4.00, -.0596]
param[18] = [4.252, .1982, 0.8180, -0.4740, 281.0, 4.00, -.0220]
# Magnesiowüstite Mg and Fe
param[19] = [3.584, .3768, 0.7404, -0.7446, 160.3, 4.13, -.0272]
param[20] = [5.865, .3203, 0.6293, -0.0000, 146.0, 4.00, -.0200]
# Iron-iron sulfide
param[21] = [7.88, 0.7700, 0.0000, -.0000, 170.0, 4.00, -.0200] # 300 К
param[22] = [4.94, 0.6852, .0000, -0.0000, 54.0, 4.00, -.0200] # 800 К
FEH = [6.7, 0, 0, 0, 121, 5.31, 0] # ？2100 ? К
# вектор молярных масс
mol = [140.73, 203.79, 140.73, 203.79, 140.73, 203.79, 200.79, 263.88, 200.79, 263.88, 216.55, 248.11, 403.19, 497.78,
       450.47, 401.64, 100.41, 131.94, 116.17, 40.32, 71.85, 55.85, 87.92]
# пересчитываем первый столбец под плотности в г/см^3
for i in range(23):
    #param[i][0] = mol[i] / param[i][0]
    param[i][1] *= 1e-4
    param[i][2] *= 1e-8

den = numpy.zeros(11)
den_min = numpy.zeros((11, 12))  # комбинирует fe- и mg- минералы
for i in range(11):
    for j in range(7):
        if mineralogy_den[i][j] != 0:
            #den_eq = param[2*j][0] + (param[2*j+1][0]-param[2*j][0])*fe_s[i][j]
            #a1_eq = param[2*j+1][1]*(1-fe_s[i][j]) + param[2*j][1]*fe_s[i][j]
            #a2_eq = param[2 * j + 1][2] * (1 - fe_s[i][j]) + param[2 * j][2] * fe_s[i][j]
            #a3_eq = param[2 * j + 1][3] * (1 - fe_s[i][j]) + param[2 * j][3] * fe_s[i][j]
            #K_eq = param[2*j+1][4] + (param[2*j][4] - param[2*j+1][4])*fe_s[i][j]
            #Kp_eq = param[2 * j + 1][5] + (param[2 * j][5] - param[2 * j + 1][5]) * fe_s[i][j]
            #Kt_eq = param[2 * j + 1][6] + (param[2 * j][6] - param[2 * j + 1][6]) * fe_s[i][j]
            #param_eq = [den_eq, a1_eq, a2_eq, a3_eq, K_eq, Kp_eq, Kt_eq]
            den_miner[i][2*j] = dihotomia(0, T_0, T_den[i], P_den[i], param[2*j])
            den_miner[i][2 * j+1] = dihotomia(0, T_0, T_den[i], P_den[i], param[2 * j+1])
            fe = fe_s[i][j] * mol[2*j+1]
            mg = (1-fe_s[i][j])*mol[2*j]
            femg = fe + mg
            fe /= femg
            den_min[i][j] = (1 - fe) / den_miner[i][2*j] + fe / den_miner[i][2*j+1]
            den_min[i][j] = 1 / den_min[i][j]
            den[i] += mineralogy_den[i][j] / den_min[i][j]
    for j in [7, 8]:
        if mineralogy_den[i][j] != 0:
            den_miner[i][7 + j] = dihotomia(0, T_0, T_den[i], P_den[i], param[7 + j])
            den_min[i][j] = den_miner[i][7 + j]
            den[i] += mineralogy_den[i][j] / den_min[i][j]
    if mineralogy_den[i][9] != 0:
        den_miner[i][16] = dihotomia(0, T_0, T_den[i], P_den[i], param[16])
        den_miner[i][17] = dihotomia(0, T_0, T_den[i], P_den[i], param[17])
        fe = fe_s[i][9] * mol[17]
        mg = (1 - fe_s[i][9]) * mol[16]
        femg = fe + mg
        fe /= femg
        den_min[i][9] = (1 - fe) / den_miner[i][16] + fe / den_miner[i][17]
        den_min[i][9] = 1 / den_min[i][9]
        den[i] += mineralogy_den[i][9] / den_min[i][9]
    if mineralogy_den[i][10] != 0:
        den_miner[i][18] = dihotomia(0, T_0, T_den[i], P_den[i], param[18])
        den_min[i][10] = den_miner[i][18]
        den[i] += mineralogy_den[i][10] / den_min[i][10]
    if mineralogy_den[i][11] != 0:
        den_miner[i][19] = dihotomia(0, T_0, T_den[i], P_den[i], param[19])
        den_miner[i][20] = dihotomia(0, T_0, T_den[i], P_den[i], param[20])
        fe = fe_s[i][11] * mol[20]
        mg = (1 - fe_s[i][11]) * mol[19]
        femg = fe + mg
        fe /= femg
        den_min[i][11] = (1 - fe) / den_miner[i][19] + fe / den_miner[i][20]
        den_min[i][11] = 1 / den_min[i][11]
        den[i] += mineralogy_den[i][11] / den_min[i][11]
    if den[i] != 0:
        den[i] = 1 / den[i]
    else:
        print("ERROR")
        exit()
    print(den[i])

# система 15: находим связь давления и глубины (блок 2)

R_Mars = 3389.92 # км
M_Mars = 6.4185e23 # кг
Gravity_const = 6.6743e-11 # в СИ
rad_core_const = 1750 #км
den_av = 3*M_Mars/(4*math.pi*R_Mars**3) * 1e-12 # г/см3
den_crust_const /= den_av
depth_crust_const /= R_Mars
rad_core_const /= R_Mars
P_const = 3*Gravity_const*M_Mars**2/(4*math.pi*R_Mars**4) * 1e-21
n_grid = 10000 # количество шагов
den_grid = numpy.zeros(n_grid)
rad_grid = numpy.linspace(0, 1, n_grid)
mass_grid = numpy.linspace(0, 1, n_grid)
P_grid = numpy.zeros(n_grid)
T_grid = numpy.zeros(n_grid)
i = n_grid-1
# 1 - кора
den_grid[i] = den_crust_const
while rad_grid[i-1] >= 1 - depth_crust_const:
    den_grid[i-1] = den_crust_const
    mass_grid[i-1] = mass_grid[i] + RK4(0, i, den_mode='given', rad=rad_grid, mass=mass_grid, den=den_grid)
    P_grid[i-1] = P_grid[i] + P_const*RK4(1, i, den_mode='given', rad=rad_grid, mass=mass_grid, den=den_grid)
    i -= 1
i_crust = i
# 2 - мантия
while rad_grid[i-1] >= rad_core_const:
    T_grid[i] = numpy.interp(P_grid[i], P, T_p,
                             left=T_p[0] - (P[0]-P_grid[i])/(P[1]-P[0])*(T_p[1]-T_p[0]))
    den_grid[i] = numpy.interp(P_grid[i], P_den, den,
                               left=den[0] - (P_den[0]-P_grid[i])/(P_den[1]-P_den[0])*(den[1]-den[0]))/den_av
    mass_grid[i-1] = mass_grid[i] + RK4(0, i, den_mode='linear approx', rad=rad_grid, mass=mass_grid,
                                        P=P_grid, den_app=den/den_av, P_app=P_den)
    P_grid[i-1] = P_grid[i] + P_const*RK4(1, i, den_mode='linear approx', rad=rad_grid, mass=mass_grid,
                                        P=P_grid, den_app=den/den_av, P_app=P_den)
    i -= 1

# 3 - ядро
n_core = 40  # точность деления ядра
i_core = i
P_b = 0
P_c = 50  # ГПа - предполагаемое максимальное давление (не строго)
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

# построим рис 2 из Жарков, Гудкова 2005

depth_grid = numpy.zeros(n_grid)
P_grid_depth = numpy.zeros(n_grid)
for i in range(n_grid):
    depth_grid[i] = (1-rad_grid[n_grid-1-i])*R_Mars
    P_grid_depth[i] = P_grid[n_grid-1-i]

def depth_to_P(x):
    return numpy.interp(x, depth_grid, P_grid_depth)

def P_to_depth(x):
    return numpy.interp(x, P_grid_depth, depth_grid)

print(rad_grid[i_core]*R_Mars)


# Скорости волн в мантии

param_speed = numpy.zeros((2, 14, 8))

param_speed[0][0] = [3.222, 129, 82, 4.2, 1.4, 0.017, 0.014, 26.6e-6]
param_speed[0][1] = [3.472, 170, 114, 4.3, 1.4, 0.018, 0.014, 22.0e-6]
param_speed[0][2] = [3.564, 186, 124, 4.1, 1.3, 0.021, 0.016, 21.0e-6]
param_speed[0][3] = [3.204, 104, 77, 5.0, 2.0, 0.012, 0.011, 27.0e-6]
param_speed[0][4] = [3.188, 114, 77, 5.0, 2.0, 0.012, 0.011, 27.0e-6]
param_speed[0][5] = [3.280, 113, 67, 4.5, 1.7, 0.013, 0.010, 27.0e-6]
param_speed[0][6] = [3.566, 175, 90, 4.9, 1.4, 0.021, 0.010, 18.0e-6]
param_speed[0][7] = [3.600, 169, 104, 4.9, 1.6, 0.016, 0.015, 16.0e-6]
param_speed[0][8] = [3.518, 175, 90, 4.9, 1.4, 0.021, 0.010, 20.0e-6]
param_speed[0][9] = [4.107, 266, 153, 3.9, 2.0, 0.031, 0.028, 17.0e-6]
param_speed[0][10] = [4.252, 227, 125, 3.9, 1.9, 0.027, 0.023, 17.0e-6]
param_speed[0][11] = [3.584, 163, 131, 4.2, 2.5, 0.016, 0.024, 32.5e-6]
param_speed[0][12] = [7.03, 105, 0, 4.5, 0, 0.025, 0, 75e-6] # fe T = 2100 K
param_speed[0][13] = [4.94, 54, 0, 4.0, 0, 0.020, 0, 68.52e-6] # FeS T = 1100 K
FEH_speed = [6.7, 5.31, 0, 121, 0, 0, 0, 0]

param_speed[1][0] = [1.18, 9, -31, 0, 0, 0, 0, 0]
param_speed[1][1] = [1.24, 15, -41, 0, 0, 0, 0, 0]
param_speed[1][2] = [1.285, 15, -41, 0, 0, 0, 0, 0]
param_speed[1][3] = [0.798, 0, -24, 0, 0, 0, 0, 0]
param_speed[1][4] = [0.817, 0, -24, 0, 0, 0, 0, 0]
param_speed[1][5] = [0.376, 7, -6, 0, 0, 0, 0, 0]
param_speed[1][6] = [0.746, 1, 8, 0, 0, 0, 0, 0]
param_speed[1][9] = [1.047, 0, 0, 0, 0, 0, 0, 0]
param_speed[1][11] = [2.281, -8, -77, 0, 0, 0, 0, 0]

v_p = numpy.zeros(n_grid)
v_s = numpy.zeros(n_grid)

for i in range(n_grid):
    den_grid[i] *= den_av

# Далее считаем скорости, в ядре:

i = 0
v_p_f = 0
v_p_s = 0
while i < i_core + 1:
    Kv = 0
    Kr = 0
    Kv_der = 0
    Kr_der = 0
    denv = 0
    denr = 0
    if ferrum:
        den_star = param_speed[0][12][0] * math.exp(-param_speed[0][12][7] * (T_grid[i]-2100))
        expon = param_speed[0][12][5]/param_speed[0][12][7]/param_speed[0][12][1]
        K = param_speed[0][12][1]*(den_star/param_speed[0][12][0])**expon
        K_der = param_speed[0][12][3] * math.exp(param_speed[0][12][7] * (T_grid[i]-2100))
        Kv += ferrum*K
        Kr += ferrum/K
        Kv_der += ferrum*K_der
        Kr_der += ferrum/K_der
        denv += ferrum*den_star
        denr += ferrum/den_star
    if sulf:
        den_star = param_speed[0][13][0] * math.exp(-param_speed[0][13][7] * (T_grid[i]-1100))
        expon = param_speed[0][13][5]/param_speed[0][13][7]/param_speed[0][13][1]
        K = param_speed[0][13][1]*(den_star/param_speed[0][13][0])**expon
        K_der = param_speed[0][13][3] * math.exp(param_speed[0][13][7] * (T_grid[i]-1100))
        Kv += sulf*K
        Kr += sulf/K
        Kv_der += sulf*K_der
        Kr_der += sulf/K_der
        denv += sulf*den_star
        denr += sulf/den_star
    if hydr:
        den_star = FEH_speed[0]
        K = FEH_speed[1]
        K_der = FEH_speed[3]
        Kv += hydr*K
        Kr += hydr/K
        Kv_der += hydr*K_der
        Kr_der += hydr/K_der
        denv += hydr*den_star
        denr += hydr/den_star
    K = (Kv + 1/Kr)/2
    K_der = (Kv_der + 1/Kr_der)/2
    den_star = (denv + 1/denr)/2
    L1 = K
    L2 = 5*K - 3*K*K_der
    epsil = (1 - (den_grid[i]/den_star)**(2/3))/2
    v_p[i] = math.sqrt( (1 - 2*epsil)**(5/2) * (L1 + L2*epsil) / den_grid[i] )
    v_s[i] = 0
    i += 1
# в мантии:
v_p_den = numpy.zeros(11)
v_s_den = numpy.zeros(11)
param_speed_cur = numpy.zeros(8)
T_0 = 300
for i in range(11):
    Kv = 0
    Kr = 0
    Gv = 0
    Gr = 0
    Kv_der = 0
    Kr_der = 0
    Gv_der = 0
    Gr_der = 0
    denv = 0
    denr = 0
    for j in range(12):
        if mineralogy_den[i][j] != 0:
            T = T_den[i]
            for k in range(8):
                param_speed_cur[k] = param_speed[0][j][k] + param_speed[1][j][k]*fe_s[i][j]
            den_star = param_speed_cur[0] * math.exp(-param_speed_cur[7] * (T-T_0))
            expon = param_speed_cur[5]/param_speed_cur[7]/param_speed_cur[1]
            K = param_speed_cur[1]*(den_star/param_speed_cur[0])**expon
            expon = param_speed_cur[6]/param_speed_cur[7]/param_speed_cur[2]
            G = param_speed_cur[2]*(den_star/param_speed_cur[0])**expon
            K_der = param_speed_cur[3] * math.exp(param_speed_cur[7] * (T-T_0))
            G_der = param_speed_cur[4] * math.exp(param_speed_cur[7] * (T-T_0))
            Kv += mineralogy_den[i][j]*K
            Kr += mineralogy_den[i][j]/K
            Gv += mineralogy_den[i][j]*G
            Gr += mineralogy_den[i][j]/G
            Kv_der += mineralogy_den[i][j]*K_der
            Kr_der += mineralogy_den[i][j]/K_der
            Gv_der += mineralogy_den[i][j]*G_der
            Gr_der += mineralogy_den[i][j]/G_der
            denv += mineralogy_den[i][j]*den_star
            denr += mineralogy_den[i][j]/den_star
    K = (Kv + 1/Kr)/2
    G = (Gv + 1/Gr)/2
    K_der = (Kv_der + 1/Kr_der)/2
    G_der = (Gv_der + 1/Gr_der)/2
    den_star = (denv + 1/denr)/2
    M1 = G
    L1 = K + 4*G/3
    M2 = 5*G - 3*K*G_der
    L2 = 5*(K + 4*G/3) - 3*K*(K_der + 4*G_der/3)
    epsil = (1 - (den[i]/den_star)**(2/3))/2
    v_p_den[i] = math.sqrt( (1 - 2*epsil)**(5/2) * (L1 + L2*epsil) / den[i] )
    v_s_den[i] = math.sqrt( (1 - 2*epsil)**(5/2) * (M1 + M2*epsil) / den[i] )
for i in range(i_core+1, i_crust):
    v_p[i] = numpy.interp(P_grid[i], P_den, v_p_den,
                          left=v_p_den[0]-(P_den[0]-P_grid[i])/(P_den[1]-P_den[0])*(v_p_den[1]-v_p_den[0]))
    v_s[i] = numpy.interp(P_grid[i], P_den, v_s_den,
                          left=v_s_den[0]-(P_den[0]-P_grid[i])/(P_den[1]-P_den[0])*(v_s_den[1]-v_s_den[0]))

# Кора
i = i_crust
while i < n_grid:
    v_p[i] = 7
    v_s[i] = 4
    i += 1

if GRAPH:
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.25, right=0.75)
    pres = ax.twinx()
    temp = ax.twinx()
    velo = ax.twinx()
    velo.yaxis.tick_left()
    velo.yaxis.set_label_position('left')
    temp.spines.right.set_position(("axes", 1.2))
    velo.spines.left.set_position(("axes", -0.2))
    pr, = ax.plot(rad_grid*R_Mars, den_grid, "C0")
    ax.text(700, den_grid[n_grid//10]+0.1, r"$\rho$")
    pp, = pres.plot(rad_grid*R_Mars, P_grid, "C1")
    pres.text(1000, P_grid[3*n_grid//10]-1, "$P$", horizontalalignment='right', verticalalignment='top')
    pt, = temp.plot(rad_grid*R_Mars, T_grid, "C2")
    temp.text(1500, T_grid[4*n_grid//10], '$T$')
    pvp, = velo.plot(rad_grid*R_Mars, v_p, "C3")
    velo.text(700, v_p[n_grid//10]+0.3, "$v_p$")
    velo.text(3000, v_p[9*n_grid // 10]-0.1, "$v_p$", horizontalalignment='right', verticalalignment='top')
    pvs, = velo.plot(rad_grid*R_Mars, v_s, "C4")
    velo.text(2200, v_s[7*n_grid//10], "$v_s$", horizontalalignment='right', verticalalignment='top')
    ax.set(xlim=(0, R_Mars), ylim=(0, 9), xlabel="Радиус, $км$", ylabel="Плотность, $10^3$ $кг/м^3$")
    pres.set(ylim=(0, 50), ylabel="Давление, $ГПа$")
    temp.set(ylim=(0, 2500), ylabel="Температура, $К$")
    velo.set(ylim=(0, 25), ylabel="Скорость, $км/с$")
    if SAVEPIC:
        plt.savefig(fname=s_out, dpi=600, bbox_inches='tight', pad_inches=0.1)
    plt.show()

# Запись данных в файл
DATA_grid = []
i_const = int(numpy.interp(0.04, rad_grid, range(n_grid)))
for i in range(i_const):
    DATA_grid.append([rad_grid[i]*R_Mars, den_grid[i], v_p[i], v_s[i], 4/3*math.pi*Gravity_const*R_Mars*1e6*rad_grid[i]*den_grid[i]])
for i in range(i_const, n_grid):
    DATA_grid.append([rad_grid[i]*R_Mars, den_grid[i], v_p[i], v_s[i], Gravity_const*M_Mars*mass_grid[i]/(R_Mars*1e3*rad_grid[i])**2])
DATA_grid.append([round(numpy.interp(P_grid[0]-P_4, P_grid[0]-P_grid[:], range(n_grid))),
                  round(numpy.interp(P_grid[0]-P_3, P_grid[0]-P_grid[:], range(n_grid))),
                  round(numpy.interp(P_grid[0]-P_2, P_grid[0]-P_grid[:], range(n_grid))),
                  round(numpy.interp(P_grid[0]-P_1, P_grid[0]-P_grid[:], range(n_grid)))])
DATA_grid.append([i_core, i_crust])
DATAs = pd.DataFrame(DATA_grid, columns=['radius', 'density', 'v_p', 'v_s', 'gravity'])

DATAs.to_excel("datas/Mars_datas.xlsx", index=False)