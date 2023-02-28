import math
import numpy as np

from num_functions import Edo


def f(x, y):
    y = y[0]
    return [-1.2*y + 7*math.exp(-0.3*x)]
    # return x+y**2


fun = Edo(f)
fun.IC = [3]
fun.N = 10
fun.intervalo = (0, 1.5)
print(fun.IC)
print(fun.N)
print(fun.intervalo)
print(fun.h)
fun.rungekutta('rk4', 'classico')
print(fun.xy)