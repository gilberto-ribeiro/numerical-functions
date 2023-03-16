import math
import numpy as np

from num_functions import Edo


def f(x, y):
    t = x
    Ca = y[0]
    Cb = y[1]
    Cc = y[2]
    Cd = y[3]
    Cain = 0.7
    Cbin = 0.4
    Ccin = 0
    Cdin = 0
    F = 5
    V = 40
    k = 0.855 * 60
    dCa = (F / V) * (Cain - Ca) - k * Ca * Cb
    dCb = (F / V) * (Cbin - Cb) - k * Ca * Cb
    dCc = (F / V) * (Ccin - Cc) + k * Ca * Cb
    dCd = (F / V) * (Cdin - Cd) + k * Ca * Cb
    return [dCa, dCb, dCc, dCd]



func = Edo(f)
func.intervalo = (0, 60)
func.IC = [0, 0, 0, 0.8]
func.N = 1000
func.titulo = 'Concentração por tempo'
func.eixo_x = 'Tempo [min]'
func.eixo_y = 'Concentração [mol/L]'
func.limites_x = (0, 120)
# func.limites_y = (0, 250)
func.legendas = ['C_A', 'C_B', 'C_C', 'C_D']
func.rungekutta()
func.gera_dados()