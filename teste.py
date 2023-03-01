import math
import numpy as np

from num_functions import Edo


def f(x, y):
    y = y[0]
    return [-1.2*y + 7*math.exp(-0.3*x)]
    # return x+y**2

def tem_atributo(classe, atributo):
    if hasattr(classe, atributo):
        print('Tem atributo')
    else:
        print('Não tem atributo')

fun = Edo(f)
fun.IC = [3]
fun.N = 500
fun.intervalo = (0, 2.5)
fun.titulo = 'Equação Diferencial Parcial'
fun.eixo_x = 'Tempo'
fun.eixo_y = 'Concentração'
fun.rungekutta('rk4')
fun.formato_grafico = ['pdf', 'png']
fun.gnuplot('teste')