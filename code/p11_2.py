import numpy as np
from scipy.constants import R
import sympy as sp
from helpers.solvers.numerical import newtons_method as nm

R_bar = R*10**(-6)

n = np.array([
    4, 2.5
])

t = np.array([
    75+273.15, 130+273.15 #K
])

p = np.array([
    30, 20 #bar
])

i = 1, 2

ntotal = np.sum(n)
x = np.array([
    n[0]/ntotal, 1 - n[0]/ntotal
])

cv = np.array([
    2.5*R, 1.5*R
])

cv = np.array([
    3.5*R, 2.5*R
])

T = sp.symbols("T")
terms = np.array([])
for i_ in i:
    j = i_ - 1
    n_, t_, cv_ = n[j], t[j], cv[j]
    ith_term = n_ * cv_ * (T - t_)
    terms = np.append(terms, ith_term)
expr = terms[0] + terms[1]
guess = (t[0] + t[1]) / 2
solution = nm(expr, T, 100, 1e-10)
print(solution)
