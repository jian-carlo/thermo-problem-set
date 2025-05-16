import numpy as np
from scipy.constants import R

R_bar = R*10**(-5)

x = np.array([
    0.7, 0.3
])

v = np.array([
    0.7, 0.3
])

i = 1, 2

P = 1 #bar

T = 298.15 #K

n = 0.0
for i_ in i:
    j = i_ - 1
    n += P * v[j] / (R_bar * T)

print(f"n = {n:.3f}")

ds = 0.0
for i_ in i:
    j = i_ - 1
    ds += -1 * n * R * x[j] * np.log(x[j]) 

print(f"Delta S = {ds:.3f}")
