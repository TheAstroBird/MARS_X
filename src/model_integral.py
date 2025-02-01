# Данная программа считает интегральные параметры для построенной модели внутреннего строения Марса
"""import math
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Нахождение массы планеты методом трапеций
def mass(den=[], rad=[], dif=[]):
    m = 0
    for i in range(len(dif)):
        m += (den[i]*rad[i]**2 + den[i+1]*rad[i+1]**2) * 2 * math.pi * dif[i]
    return m

# Нахождение момента инерции планеты методом трапеций
def inertia(den=[], rad=[], dif=[]):
    I = 0
    for i in range(len(dif)):
        I += (den[i]*rad[i]**4 + den[i+1]*rad[i+1]**4) * 4 / 3 * math.pi * dif[i]
    return I

# Нахождение числа Лява k2 - правая часть СДУ
def RHS(n, ie, den, rad, lame1, lame2, grav, zs=[]):
    f = [0, 0, 0, 0, 0, 0]
    for d in range(6):
        f[d] = zs[d]
    f[0] = - 2*lame1/(lame1+2*lame2)*zs[0]/rad + zs[1]/(lame1+2*lame2) + lame1*n*(n+1)/(lame1+2*lame2)*zs[4]/rad
    f[1] = ( (-4*den*grav*rad+4*lame2*(3*lame1+2*lame2)/(lame1+2*lame2))*zs[0]/(rad**2)
             - 4*lame2/(lame1+2*lame2)*zs[1]/rad
             + (n*(n+1)*den*grav*rad-2*n*(n+1)*lame2*(3*lame1+2*lame2)/(lame1+2*lame2))*zs[4]/(rad**2) + n*(n+1)/rad*zs[5]
             - den*zs[3] )
    f[2] = 3*den*zs[0] + zs[3]
    f[3] = - 3*n*(n+1)*den/rad*zs[4] + n*(n+1)/(rad**2)*zs[2] - 2/rad*zs[3]
    f[4] = - zs[0]/rad + zs[4]/rad + zs[5]/lame2
    f[5] = ( (den*grav/rad-2*lame2*(3*lame1+2*lame2)/(lame1+2*lame2)/(rad**2))*zs[0] - lame1/(lame1+2*lame2)*zs[1]/rad
             + 2*lame2/(lame1+2*lame2)*((2*n**2+2*n-1)*lame1+2*(n**2+n-1)*lame2)*zs[4]/(rad**2) - 3*zs[5]/rad - den*zs[2]/rad )
    return f[ie]


# Нахождение числа Лява k2 - СДУ
def SDE(n, i_core, z0=[], den=[], rad=[], dif=[], lame1=[], lame2=[], grav=[]):
    z = [0, 0, 0, 0, 0, 0]
    zs = [0, 0, 0, 0, 0, 0]
    zsint = [0, 0, 0, 0, 0, 0]
    k1 = [0, 0, 0, 0, 0, 0]
    k2 = [0, 0, 0, 0, 0, 0]
    k3 = [0, 0, 0, 0, 0, 0]
    k4 = [0, 0, 0, 0, 0, 0]
    for d in range(6):
        z[d] = z0[d]
    for i in range(i_core + 1, len(rad)-1):
        for d in range(6):
            zs[d] = z[d]
            zsint[d] = z[d]
        for d1 in range(6):
            k1[d1] = RHS(n, d1, den[i], rad[i], lame1[i], lame2[i], grav[i], zsint)
        for j in range(6):
            zsint[j] = zs[j]+dif[i]/2*k1[j]
        for d1 in range(6):
            k2[d1] = RHS(n, d1, (den[i]+den[i+1])/2, (rad[i]+rad[i+1])/2, (lame1[i]+lame1[i+1])/2, (lame2[i]+lame2[i+1])/2,
                 (grav[i]+grav[i+1])/2, zsint)
        for j in range(6):
            zsint[j] = zs[j]+dif[i]/2*k2[j]
        for d1 in range(6):
            k3[d1] = RHS(n, d1, (den[i]+den[i+1])/2, (rad[i]+rad[i+1])/2, (lame1[i]+lame1[i+1])/2, (lame2[i]+lame2[i+1])/2,
                 (grav[i]+grav[i+1])/2, zsint)
        for j in range(6):
            zsint[j] = zs[j]+dif[i]*k3[j]
        for d1 in range(6):
            k4[d1] = RHS(n, d1, den[i+1], rad[i+1], lame1[i+1], lame2[i+1], grav[i+1], zsint)
        for j in range(6):
            z[j] += dif[i]*(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6
    return z


# Нахождение числа Лява k2
def Love(load, n, i_core, gamma, love=[], den=[], rad=[], dif=[], lame1=[], lame2=[], grav=[]):
    # инициируем начальные векторы
    z1 = [1, den[i_core] * grav[i_core], 0, -3 * den[i_core], 0, 0]
    z2 = [0, -den[i_core], 1, n / rad[i_core] + gamma, 0, 0]
    z3 = [0, 0, 0, 0, 1, 0]

    # прогоним векторы через СДУ
    z1 = SDE(n, i_core, z1, den, rad, dif, lame1, lame2, grav)
    z2 = SDE(n, i_core, z2, den, rad, dif, lame1, lame2, grav)
    z3 = SDE(n, i_core, z3, den, rad, dif, lame1, lame2, grav)

    # решим уравнения на поверхности методом Гаусса
    A = [[z1[5], z2[5], z3[5]], [z1[1], z2[1], z3[1]], [z1[3]+(n+1)*z1[2], z2[3]+(n+1)*z2[2], z3[3]+(n+1)*z3[2]]]
    b = [0, 0, 2*n+1]
    if load:
        b[1] = -(2*n+1)/3
    for i in range(1, 3):
        A[1][i] -= A[0][i] * A[1][0] / A[0][0]
        A[2][i] -= A[0][i] * A[2][0] / A[0][0]
    A[2][2] -= A[1][2] * A[2][1] / A[1][1]
    b[2] -= b[1] * A[2][1] / A[1][1]
    b[2] /= A[2][2]
    b[1] -= b[2] * A[1][2]
    b[1] /= A[1][1]
    b[0] -= b[1] * A[0][1] + b[2] * A[0][2]
    b[0] /= A[0][0]

    # и само число Лява k2
    if love == 'k':
        lov = b[0]*z1[2] + b[1]*z2[2] + b[2]*z3[2] - 1
    elif love == 'h':
        lov = b[0] * z1[0] + b[1] * z2[0] + b[2] * z3[0]
    elif love == 'l':
        lov = b[0] * z1[4] + b[1] * z2[4] + b[2] * z3[4]
    return lov

# --------------------------------------------------------------------------------

# Главное тело программы

# --------------------------------------------------------------------------------

# Режим работы программы

VISCOSITY = True
MELTLAYER = False
CREEPFUNCTION = False
REWRITEFILE = False

# Вводим данные

n = 2
DATAs = pd.read_excel("datas/Mars_datas.xlsx")
n_grid = DATAs.shape[0] - 2
rad = np.zeros(n_grid)
den = np.zeros(n_grid)
v_p = np.zeros(n_grid)
v_s = np.zeros(n_grid)
grav = np.zeros(n_grid)
for i in range(n_grid):
    rad[i] = DATAs['radius'][i]
    den[i] = DATAs['density'][i]
    v_p[i] = DATAs['v_p'][i]
    v_s[i] = DATAs['v_s'][i]
    grav[i] = DATAs['gravity'][i]

DATAadd = pd.read_excel("datas/Input_data.xlsx")
sulf = DATAadd['sulfur']
hydr = DATAadd['hydro']

rad_av = 3389.92    # км
h_melt = 200 # км - толщина расплавленного слоя
den_av = 3.935      # г/см^3
mass_av = 4/3*math.pi*(1000*den_av) * (1000*rad_av)**3    # кг
grav_av = 3.7279        # м/с^2
transit4 = int(DATAs['radius'][n_grid])
transit3 = int(DATAs['density'][n_grid])
transit2 = int(DATAs['v_p'][n_grid])
transit1 = int(DATAs['v_s'][n_grid])
i_core = int(DATAs['radius'][n_grid+1])
i_crust = int(DATAs['density'][n_grid+1])
den_crust = den[n_grid-1]
den_core = den[0]
if MELTLAYER:
    di_melt = int(h_melt/rad_av*n_grid)
else:
    di_melt = 0

alpha = DATAadd['andrade'][0]
eta_0 = DATAadd['viscosity'][0]/1e9 # коры, ГПа-с
sigmaw = 2*math.pi/206.9/86400
chi = 2*math.pi/44340 # 1/с
eta_melt = DATAadd['eta melt layer'][0]
eta_layers = [eta_melt, math.log(eta_0, 10), math.log(eta_0, 10), math.log(eta_0, 10)-1, math.log(eta_0, 10)-1, math.log(eta_0, 10)-2,
              math.log(eta_0, 10)-2, math.log(eta_0, 10)]
i_layers = [i_core+di_melt, i_core+di_melt, transit4, transit3, transit2, transit1, i_crust, i_crust]

lame1 = []
lame2 = []
Q = []
if VISCOSITY:
    for i in range(len(v_p)):
        if v_s[i] == 0:
            lame2.append(0)
            lame1.append((v_p[i]**2)*den[i])
            Q.append(math.inf)
        else:
            eta = 10**np.interp(i, i_layers, eta_layers)
            mu = (v_s[i]**2)*den[i]
            J = (1+(1j*eta/mu*chi)**(-alpha)*math.gamma(1+alpha))/mu-1j/(eta*chi)
            lame2.append(1/J)
            K = (v_p[i]**2)*den[i] - 4/3*lame2[i]
            lame1.append(K-2/3/J)
            mu = 1/J
            Q.append(mu.real/mu.imag)
else:
    for i in range(len(v_p)):
        if v_s[i] == 0:
            lame2.append(0)
            lame1.append((v_p[i]**2)*den[i])
            Q.append(math.inf)
        else:
            lame2.append((v_s[i]**2)*den[i])
            lame1.append((v_p[i]**2)*den[i] - 4/3*lame2[i])
            Q.append(math.inf)

lame1_ch = []
lame2_ch = []
if VISCOSITY:
    for i in range(len(v_p)):
        if v_s[i] == 0:
            lame2_ch.append(0)
            lame1_ch.append((v_p[i]**2)*den[i])
        else:
            eta = 10**np.interp(i, i_layers, eta_layers)
            mu = (v_s[i]**2)*den[i]
            J = (1+(1j*eta/mu*sigmaw)**(-alpha)*math.gamma(1+alpha))/mu-1j/(eta*sigmaw)
            lame2_ch.append(1/J)
            K = (v_p[i]**2)*den[i] - 4/3*lame2_ch[i]
            lame1_ch.append(K-2/3/J)
else:
    for i in range(len(v_p)):
        if v_s[i] == 0:
            lame2_ch.append(0)
            lame1_ch.append((v_p[i] ** 2) * den[i])
        else:
            lame2_ch.append((v_s[i] ** 2) * den[i])
            lame1_ch.append((v_p[i] ** 2) * den[i] - 4 / 3 * lame2_ch[i])

# обезразмериваем
for i in range(len(den)):
    rad[i] /= rad_av
    den[i] /= den_av
    lame1[i] /= den_av*rad_av*grav_av/1000
    lame2[i] /= den_av*rad_av*grav_av/1000
    grav[i] /= grav_av
    lame1_ch[i] /= den_av*rad_av*grav_av/1000
    lame2_ch[i] /= den_av*rad_av*grav_av/1000

dif = []
for i in range(len(rad)-1):
    dif.append(rad[i+1]-rad[i])

def ro(x):
    return np.interp(x, rad, den)
def gravit(x):
    return np.interp(x, rad, grav)
def Molod(x, y):
    return -y**2-2*(n+1)/x*y-3*(ro(x+0.0001)-ro(x-0.0001))/0.0002/gravit(x)
solution = sp.integrate.RK45(Molod, 0.02, [0], rad[i_core]-0.00011)
while solution.status != 'finished':
    solution.step()
gamm = solution.y[0]
print('gamma =', gamm)

Ma = mass(den, rad, dif)
In = inertia(den, rad, dif)
K2 = Love(False, n, i_core, gamm, 'k', den, rad, dif, lame1, lame2, grav)

# возвращаем размерность
Ma *= 3/4/math.pi
In *= 3/4/math.pi
print('M =', Ma)
print('I =', In)
print('k2 =', K2)

# Расчет чандлеровского периода
A = 0.362976
B = 0.363229
C = 0.365067
T = 24.6229*3600
G = 6.6743e-11

K2_ch = Love(False, 2, i_core, gamm, 'k', den, rad, dif, lame1_ch, lame2_ch, grav)
print('K2_ch =', K2_ch)
alp = (C-B)/A
bet = (C-A)/B
AB = (A+B)/2
omega = 2*math.pi/T
Te = T/math.sqrt(alp*bet)
k0 = 3*G*(C-AB)*mass_av/omega**2/(rad_av*1e3)**3
Ac = 8/15*math.pi*den_core*1000*(rad[i_core]*rad_av*1000)**5/mass_av/(rad_av*1000)**2
print('R_core =', rad[i_core]*rad_av)

Tw = Te*(1-Ac/math.sqrt(A*B))/(1-K2_ch.real/k0)/86400
print('Tw0 = ', Tw)

if CREEPFUNCTION:
    Tw_cr = Te * (1 - Ac / math.sqrt(A * B)) / (1 - K2.real / k0) / 86400
    n = 0.4
    dT = Tw_cr/1118*((chi/sigmaw)**n-1)/(k0-K2_ch.real)/math.tan(n*math.pi/2)

Data_grid = [[In, K2.real, K2.imag, den_crust, den_core, float(sulf), float(hydr), rad[i_core]*rad_av, DATAadd['viscosity'][0],
              alpha, float(10)**(eta_melt+9), Tw]]
if REWRITEFILE:
    DATAsOUT = pd.DataFrame(Data_grid, columns=['inertia', 'Re k2', 'Im k2', 'crust density', 'core density', 'sulfur',
                                                'hydro', 'core radius', 'viscosity', 'andrade', 'eta melt layer', 'Tw'])
else:
    DATAsOUTold = pd.read_excel("datas/Probe_models.xlsx")
    DATAsOUTnew = pd.DataFrame(Data_grid, columns=['inertia', 'Re k2', 'Im k2', 'crust density', 'core density',
                                                   'sulfur', 'hydro', 'core radius', 'viscosity', 'andrade',
                                                   'eta melt layer', 'Tw'])
    DATAsOUT = pd.concat([DATAsOUTold, DATAsOUTnew], ignore_index=True)
DATAsOUT.to_excel("datas/Probe_models.xlsx", index=False)

Q_data = []
for i in range(n_grid):
    Q_data.append([rad[i]*rad_av, Q[i]])
DATAsOUT2 = pd.DataFrame(Q_data, columns=['radius', 'Q_mu'])
DATAsOUT2.to_excel("datas/Q_fourlayer_data.xlsx", index=False)"""