import numpy as np
import matplotlib.pyplot as plt
import math
import json


with open('atmosphere.json', 'r', encoding='utf-8') as f:
    atmo = json.load(f)
pi = math.pi
R = 8.314462618
E = 1 * 10 ** 6
R0 = 0.3
d0 = 1/1000


def rad(pres, temp, nu):
    global E, R0, d0, pi, R
    if nu < 4 / 3 * pi * R0 ** 3 * pres / (temp * R):
        raise ValueError("Слишком маленькое значение количества гелия для заданных параметров, отрицательное давление")

    p = [4 * pi * pres, 0, 8 * pi * E * d0 * R0, -8 * pi * E * d0 * R0 ** 2 - 3 * nu * R * temp]
    rs = np.roots(p)
    real = [x.real for x in rs if abs(x.imag)<1e-6 and x.real>R0]
    return float(real[0])


def d(r):
    global R0, d0
    return R0 ** 2 / r ** 2 * d0
