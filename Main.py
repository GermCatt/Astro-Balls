import numpy as np
import matplotlib.pyplot as plt
import math
import json


with open('atmosphere.json', 'r', encoding='utf-8') as f:
    atmo = json.load(f)
pi = math.pi
R = 8.314462618

def radius(temp, pressure, nu, sigma):
    q = nu * R * temp
    Ds = np.sqrt(q * (81 * pressure ** 2 * q - 1024 * pi * sigma ** 3) / (576 * pi ** 2 * pressure ** 4))
    alpha = np.cbrt(3 * q / (8 * pi * pressure) - 64 * sigma ** 3 / (27 * pressure ** 3) + Ds)
    beta = np.cbrt(3 * q / (8 * pi * pressure) - 64 * sigma ** 3 / (27 * pressure ** 3) - Ds)
    return alpha + beta - 4 * sigma / (3 * pressure)