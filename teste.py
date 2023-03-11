import math
import numpy as np

from num_functions import Edo


def f(x, y):
    t = x
    Ca = y[0]
    Cb = y[1]
    Cc = y[2]
    T = y[3]
    Tj = y[4]
    rho = 880
    cp = 1750
    A = 5
    V = 40
    k0 = 8.2e5
    Ea = 48500
    R = 8.314
    DHr = -72800
    U = 680
    Vj = 0.032
    cpj = 4180
    rhoj = 1000
    Q = 3
    Qj = 0.01
    Cain = 200
    Cbin = 200
    Ccin = 0
    Tin = 300 if t < 20 else 290
    Tjin = 280
    k = k0 * math.exp(-Ea/(R*T))
    r = k * Ca * Cb
    dCa = (1/V) * (Q * (Cain - Ca) - r * V)
    dCb = (1/V) * (Q * (Cbin - Cb) - r * V)
    dCc = (1/V) * (Q * (Ccin - Cc) + r * V)
    dT = (1/(V*rho*cp)) * (Q*rho*cp*(Tin - T) - U*A*(T - Tj) + r*V*(-DHr))
    dTj = (1/(Vj*rhoj*cpj)) * (Qj*rhoj*cpj*(Tjin - Tj) - U*A*(Tj - T))
    return [dCa, dCb, dCc, dT, dTj]



func = Edo(f)
func.intervalo = (0, 100)
func.IC = [49.5, 49.5, 150.5, 307, 282]
func.N = 1000
func.titulo = 'Concentração e temperatura por tempo'
func.eixo_x = 'Tempo [min]'
func.eixo_y = 'Concentração [mol/m^3]'
func.limites_x = (0, 25)
# func.limites_y = (0, 40)
func.legendas = ['C_A', 'C_B', 'C_C', 'T', 'T_j']
func.rungekutta()
func.gera_dados()