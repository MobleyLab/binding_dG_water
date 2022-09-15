import sys, os

from scipy.integrate import quad
from scipy.integrate import simps, trapz, romb
import numpy as np
import matplotlib.pyplot as plt
from math import erf

plt.rcParams['text.usetex'] = True


#RT =0.5924849497137634 kcal/mol
#k = 20 kcal/mol*A^2
#C = 55M = 55 mol/L = 55 *6.023*10^23/10^27 A^3 =  0.033121 A^-3

def numerical(x, r=1.0):
    return (x**2.)*np.exp(-20./(2*0.5924849497137634)*((x-r)**2.))

def analytic(a, r):
    r2 = r*r
    A = 0.5*r/a*np.exp(-a*r2)
    B = (r2 + 0.5/a)
    C = np.sqrt(np.pi/4./a) * (1 + erf(np.sqrt(a)*r))
    return A + B*C

RT = 0.5924849497137634
k = 20.0
a = k/2.0/RT
C = 0.033121

numer = []
exact = []

r0_arr = np.linspace(0, 3.0, 100)
for r0 in r0_arr:

    y = RT*np.log(C*4*np.pi*quad(numerical, 0.0, np.inf, r0)[0])
    numer.append(y)

    y = RT*np.log(C*4*np.pi*analytic(a, r0))
    exact.append(y)
plt.plot(r0_arr, numer, label="numerical")
plt.plot(r0_arr, exact, label="analytic")
plt.xlabel(r"$r_0$ / \AA")
plt.ylabel(r"$dG_{res}$")
plt.legend()
plt.show()

