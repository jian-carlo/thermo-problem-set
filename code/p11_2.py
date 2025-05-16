import numpy as np
from scipy.constants import R
import sympy as sp
from helpers.solvers.numerical import newtons_method as nm
from helpers.solvers.numerical import newtons_method_debug as nmd

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

cp = np.array([
    3.5*R, 2.5*R
])

T = sp.symbols("T")
T_terms = np.array([])
for i_ in i:
    j = i_ - 1
    n_, t_, cv_ = n[j], t[j], cv[j]
    ith_term = n_ * cv_ * (T - t_)
    T_terms = np.append(T_terms, ith_term)
T_expr = T_terms[0] + T_terms[1]
T_guess = (t[0] + t[1]) / 2
T_value = nm(T_expr, T, T_guess)

P = sp.symbols("P")
sum_term = 0
for i_ in i :
    j = i_ - 1
    n_, t_, p_ = n[j], t[j], p[j]
    ith_term = (n_*R_bar*t_)/p_
    sum_term += ith_term
P_expr = (np.sum(n)*R_bar*T_value)/P - sum_term
P_guess = (p[0] + p[1]) / 2
P_value = nm(P_expr, P, P_guess) #bar

ds = np.array([])
for i_ in i:
    j = i_ - 1
    n_, cp_, t_, p_ = n[j], cp[j], t[j], p[j]
    dsi = n_*(cp_*np.log(T_value/t_)-R*np.log(P_value/p_))
    ds = np.append(ds, dsi)
