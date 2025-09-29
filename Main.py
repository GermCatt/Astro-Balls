import numpy as np
import matplotlib.pyplot as plt
from math import pi
import json


with open('atmosphere.json', 'r', encoding='utf-8') as f:
    atmo = json.load(f)


def rad(pres, temp, nu):        # Функция возвращает радиус шара
    global E, R0, d0, pi, R
    if nu < 4 / 3 * pi * R0 ** 3 * pres / (temp * R):
        raise ValueError("Слишком маленькое значение количества гелия для заданных параметров, отрицательное давление")

    p = [4 * pi * pres, 0, 8 * pi * E * d0 * R0, -8 * pi * E * d0 * R0 ** 2 - 3 * nu * R * temp]
    rs = np.roots(p)
    real = [x.real for x in rs if abs(x.imag) < 1e-6 and x.real > R0]
    return float(real[0])


def d(r):       # Функция возвращает толщину стенки шара
    global R0, d0
    return R0 ** 2 / r ** 2 * d0


def T(pres, temp, nu, R1, Q):       # Функция возвращает количество температуру, которая будет внутри шара после нагрева / охлаждения
    global E, R0, d0, R
    p = [4 * pi * pres * R1 ** 3 + 2 * nu * R * temp, 0, 8 * pi * E * d0 * R1 ** 3, -8 * pi * E * d0 * R1 ** 3 * R0 ** 2 - 2 * Q * R1 ** 3 - 5 * nu * R * temp * R1 ** 3]
    rs = np.roots(p)
    real = [x.real for x in rs if abs(x.imag) < 1e-6 and x.real > R0]
    if not real:
        raise ValueError("Неверный ввод, хз чё не так, пойти не так могло что угодно")
    real = float(real[0])
    T2 = 2 * Q / (3 * nu * R) + 5 / 3 * temp - 2 / 3 * temp * (real / R1) ** 3
    return T2


def P(T_atmo, temp, R1):        # Функция возвращает мощность нагрева/охлаждения газа внутри шара
    global kapa
    power = (T_atmo - temp) * 4 * pi * kapa * R1 ** 2 / d(R1)
    return power


def a(velo, R1, mu, nu, rho, g):
    global m0, M
    Re = rho * velo * R1 * 2 / mu
    if Re < 1:
        Cd = 24 / (Re + 0.001)
    elif Re > 2 * 10 ** 5:
        Cd = 0.1
    else:
        Cd = 0.47
    Fs = 0.5 * pi * R1 ** 2 * Cd * rho * velo ** 2
    Farh = 4 / 3 * pi * R1 ** 3 * g * rho
    Fg = (M + m0 + nu * 0.004) * g
    return (Farh - Fs - Fg) / (M + m0 + nu * 0.004)


h0 = 0
temp0 = atmo['0']['T']
pres0 = atmo['0']['P']
Rcrit = 2.9
velo = 7
epsV = 2
H = 40000
epsH = 500
dt = 1
M = 10
m0 = 0.5
R = 8.314462618
E = 1 * 10 ** 6
R0 = 0.3
d0 = 1/1000
kapa = 1
R1 = R0 + 0
nu = 8 / 3 * pi * R0 ** 3 * pres0 / (temp0 * R)
res = {}
for i in range(50000):
    velos = []
    atts = []
    times = []
    velo = 0
    temp = temp0 + 0
    h = 1
    t = 0
    while R1 < Rcrit:
        pres = (atmo[str(int(h // 100 * 100))]['P'] * (h % 100) + atmo[str(int(h // 100 * 100 + 100))]['P'] * (100 - h % 100)) / 100
        T_atmo = (atmo[str(int(h // 100 * 100))]['T'] * (h % 100) + atmo[str(int(h // 100 * 100 + 100))]['T'] * (100 - h % 100)) / 100
        rho = (atmo[str(int(h // 100 * 100))]['rho'] * (h % 100) + atmo[str(int(h // 100 * 100 + 100))]['rho'] * (100 - h % 100)) / 100
        g = (atmo[str(int(h // 100 * 100))]['g'] * (h % 100) + atmo[str(int(h // 100 * 100 + 100))]['g'] * (100 - h % 100)) / 100
        mu = (atmo[str(int(h // 100 * 100))]['mu'] * (h % 100) + atmo[str(int(h // 100 * 100 + 100))]['mu'] * (100 - h % 100)) / 100
        velos.append(velo)
        atts.append(h)
        times.append(t)
        t += dt
        h += velo * dt
        if h < 0:
            break
        velo += a(velo, R1, mu, nu, rho, g) * dt
        R1 = rad(pres, temp, nu)
        Q = P(T_atmo, temp, R1) * dt
        try:
            temp = T(pres, temp, nu, R1, Q)
        except ValueError:
            break
    if abs(h - H) < epsH:
        res[nu] = (times, atts, velos)
    nu += 0.2

print(res)