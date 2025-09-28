import numpy as np
import matplotlib.pyplot as plt
import math
import json


with open('atmosphere.json', 'r', encoding='utf-8') as f:
    atmo = json.load(f)