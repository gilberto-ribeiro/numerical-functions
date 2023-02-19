import math


#FUNÇÕES COMPLEMENTARES


def subint(a, b, N):
    #Função que retorna uma lista com N subintervalos interitados
    #entre a e b.
    #
    #subint(a, b, N)
    #
    #a e b são os interites do intervalo.
    #N é o número de subintervalos.
    #
    #Desenvolvida por: Gilberto Ribeiro Pinto Júnior.
    #Última atualização: 30/01/2023.
    N = int(N)
    h = (b-a)/N
    x = [float(a)]
    for i in range(N):
        x.append(x[i]+h)
    intervals = [(x[i],x[i+1]) for i in range(N)]
    return {'points': x, 'intervals': intervals, 'h': h}


def subint_tab(x):
    intervals = [(x[i],x[i+1]) for i in range(len(x)-1)]
    h = [abs(x[1]-x[0]) for x in intervals]
    return {'intervals': intervals, 'h': h}


def proxmultiplo(N, d):
    if N%d == 0:
        return N
    else:
        return N+d-N%d


def transpor(matrix):
    # Fonte: https://docs.python.org/3/tutorial/datastructures.html
    n = len(matrix[0])
    transposed = [[row[i] for row in matrix] for i in range(n)]
    return transposed


def dados_edo(data, file='ode'):
    x = data['x']
    y = data['y']
    n = len(x)
    with open('data_'+file, 'w') as fhand:
        for i in range(n):
            row = [x[i]]
            row.extend(y[i])
            row = [str(j) for j in row]
            row = ' '.join(row)+'\n'
            fhand.write(row)
        fhand.close()


def gnuplot_edo(data, file='ode', p={'title': 'Resolução das EDOs', 'axes': ['x', 'y']}):
    dados_edo(data, file)
    x = data['x']
    y = data['y']
    n = len(y[0])
    end = (n-1)*[', ']+['\n']
    with open('gr_'+file+'.gp', 'w', encoding="utf-8") as fhand:
        if 'title' in p:
            fhand.write('set title "{}"\n'.format(p['title']))
        if 'axes' in p:
            fhand.write('set xlabel "{}"\n'.format(p['axes'][0]))
            fhand.write('set ylabel "{}"\n'.format(p['axes'][1]))
        if 'xrange' in p:
            fhand.write('set xrange [{}:{}]\n'.format(str(p['xrange'][0]), str(p['xrange'][1])))
        else:
            fhand.write('set xrange [{}:{}]\n'.format(str(x[0]), str(x[-1])))
        if 'yrange' in p:
            fhand.write('set yrange [{}:{}]\n'.format(str(p['yrange'][0]), str(p['yrange'][1])))
        fhand.write('set terminal pdfcairo\n')
        fhand.write('set output "fig_{}.pdf"\n'.format(file))
        for i in range(n):
            if 'labels' in p and i <= len(p['labels']) - 1:
                label = p['labels'][i] # REVISAR
            else:
                label = 'y[{}]'.format(i)
            if i == 0:
                fhand.write('plot "data_{}" u 1:{} t "{}" w l{}'.format(file, str(i+2), label, end[i]))
            else:
                fhand.write('"" u 1:{} t "{}" w l{}'.format(str(i+2), label, end[i]))
        fhand.write('set output\n')
        fhand.write('set terminal wxt\n')
        fhand.write('replot')
        fhand.close()
        

#FUNÇÕES PARA RESOLUÇÃO DE EQUAÇÕES NÃO LINEARES


def bissec(f, a, b, tolmax=1e-8, imax=100):
    #Função para a resolução de equações não lineares pelo Método da Bissecção.
    #
    #bissec(f, a, b, tolmax, imax)
    #
    #f é a função a ser resolvida: f(x) = 0.
    #a e b são os interites do intervalo.
    #tolmax é a tolerância máxima determinada por abs(b[i]-a[i])/2.
    #imax é o número máximo de iterações até a convergência.
    #
    #Desenvolvida por: Gilberto Ribeiro Pinto Júnior.
    #Última atualização: 30/01/2023.
    a = [a]
    b = [b]
    if f(a[0])*f(b[0]) > 0:
        retval = math.nan
    else:
        x = list()
        tol = list()
        for i in range(imax):
            x.append((a[i]+b[i])/2)
            if f(x[i])*f(a[i]) < 0:
                a.append(a[i])
                b.append(x[i])
            else:
                a.append(x[i])
                b.append(b[i])
            tol.append(abs(b[i]-a[i])/2)
            if tol[i] <= tolmax:
                break
        retval = x[-1]
    return retval


def regfalsi(f, a, b, err=1e-8, imax=100):
    #Função para a resolução de equações não lineares pelo Método da Bissecção.
    #
    #regfalsi(f, a, b, tolmax, imax)
    #
    #f é a função a ser resolvida: f(x) = 0.
    #a e b são os interites do intervalo.
    #err é erro relativo estimado determinado por abs((x[i]-x[i-1])/x[i-1]).
    #imax é o número máximo de iterações até a convergência.
    #
    #Desenvolvida por: Gilberto Ribeiro Pinto Júnior.
    #Última atualização: 30/01/2023.
    a = [a]
    b = [b]
    if f(a[0])*f(b[0]) > 0:
        retval = math.nan
    else:
        x = list()
        errel = [math.nan]
        for i in range(imax):
            x.append((a[i]*f(b[i])-b[i]*f(a[i]))/(f(b[i])-f(a[i])))
            if f(x[i])*f(a[i]) < 0:
                a.append(a[i])
                b.append(x[i])
            else:
                a.append(x[i])
                b.append(b[i])
            if i > 0:
                errel.append(abs((x[i]-x[i-1])/x[i-1]))
            if errel[i] <= err:
                break
        retval = x[-1]
    return retval


def nr(f, df, xest, err=1e-8, imax=100):
    #Função para a resolução de equações não lineares pelo Método de Newton-Rhapson.
    #
    #nr(f, df, xest, err, imax)
    #
    #f é a função a ser resolvida: f(x) = 0.
    #df é a derivada de f.
    #xest é a estimativa inicial para raiz.
    #err é erro relativo estimado determinado por abs((x[i]-x[i-1])/x[i-1]).
    #imax é o número máximo de iterações até a convergência.
    #
    #Desenvolvida por: Gilberto Ribeiro Pinto Júnior.
    #Última atualização: 30/01/2023.
    x = [xest]
    errel = [math.nan]
    for i in range(imax):
        x.append(x[i]-f(x[i])/df(x[i]))
        if i > 0:
            errel.append(abs((x[i]-x[i-1])/x[i-1]))
        if errel[i] <= err:
            break
    retval = x[-1]
    return retval


def sec(f, x1, x2, err=1e-8, imax=100):
    #Função para a resolução de equações não lineares pelo Método da Secante.
    #
    #sec(f, x1, x2, err, imax)
    #
    #f é a função a ser resolvida: f(x) = 0.
    #x1 e x2 são os pontos iniciais na vizinhança da raiz.
    #err é erro relativo estimado determinado por abs((x[i]-x[i-1])/x[i-1]).
    #imax é o número máximo de iterações até a convergência.
    #
    #Desenvolvida por: Gilberto Ribeiro Pinto Júnior.
    #Última atualização: 30/01/2023.
    x = [x1, x2]
    errel = 3*[math.nan]
    for i in range(imax+1):
        if i > 0:
            df = (f(x[i])-f(x[i-1]))/(x[i]-x[i-1])
            x.append(x[i]-f(x[i])/df)
        if i > 2:
            errel.append(abs((x[i]-x[i-1])/x[i-1]))
        if errel[i] <= err:
            break
    retval = x[-1]
    return retval


def pfixo(g, xest, err=1e-8, imax=100):
    #Função para a resolução de equações não lineares pelo Método do Ponto Fixo.
    #
    #pfixo(g, xest, err, imax)
    #
    #g é a função de iteração, em que: g(x) = x.
    #xest é a estimativa inicial para raiz.
    #err é erro relativo estimado determinado por abs((x[i]-x[i-1])/x[i-1]).
    #imax é o número máximo de iterações até a convergência.
    #
    #Desenvolvida por: Gilberto Ribeiro Pinto Júnior.
    #Última atualização: 30/01/2023.
    x = [xest]
    errel = [math.nan]
    for i in range(imax):
        x.append(g(x[i]))
        if i > 0:
            errel.append(abs((x[i]-x[i-1])/x[i-1]))
        if errel[i] <= err:
            break
    retval = x[-1]
    return retval


def raizes(func, a, b, N=100, method=bissec, err=1e-8, imax=100):
    #Função que retorna uma lista com as raízes reais de f(x)=0
    #encontradas entre a e b.
    #
    #Desenvolvida por: Gilberto Ribeiro Pinto Júnior.
    #Última atualização: 30/01/2023.
    inter = subint(a, b, N)['intervals']
    x = list()
    for i in inter:
        root = method(func, i[0], i[1], err, imax)
        if not math.isnan(root):
            x.append(root)
    return x


#FUNÇÕES PARA AJUSTE DE CURVAS E INTERPOLAÇÃO


def reglin(x, y):
    n = len(x)
    Sx = sum(x)
    Sy = sum(y)
    Sxy = sum([x[i]*y[i] for i in range(n)])
    Sxx = sum([x[i]*x[i] for i in range(n)])
    a1 = (n*Sxy-Sx*Sy)/(n*Sxx-Sx*Sx)
    a0 = (Sxx*Sy-Sxy*Sx)/(n*Sxx-Sx*Sx)
    yc = [a1*x+a0 for x in x]
    return {'coeff': (a1, a0), 'yc': yc}
    

#FUNÇÕES PARA INTEGRAÇÃO NUMÉRICA


def ret(f, a, b, N=100):
    subs = subint(a, b, N)
    intervals = subs['intervals']
    h = subs['h']
    I = list()
    for x in intervals:
        I.append(h*f(x[0]))
    Itot = sum(I)
    return Itot


def pcentral(f, a, b, N=100):
    subs = subint(a, b, N)
    intervals = subs['intervals']
    h = subs['h']
    I = list()
    for x in intervals:
        I.append(h*f((x[0]+x[1])/2))
    Itot = sum(I)
    return Itot


def trap(f, a, b, N=100):
    subs = subint(a, b, N)
    intervals = subs['intervals']
    h = subs['h']
    # I = list()
    # for x in intervals[1:]:
    #     I.append(h*f((x[0])))
    # Itot = (h/2)*(f(a)+f(b))+sum(I)
    I = [(h/2)*(f(x[0])+f(x[1])) for x in intervals]
    Itot = sum(I)
    return Itot


def simp13(f, a, b, N=100):
    d = 2
    Nd = proxmultiplo(N, d)/d
    subs = subint(a, b, Nd)
    intervals1 = subs['intervals']
    h = subs['h']/d
    intervals2 = [tuple(subint(i[0], i[1], d)['points']) for i in intervals1]
    I = [(h/3)*(f(x[0])+4*f(x[1])+f(x[2])) for x in intervals2]
    Itot = sum(I)
    return Itot


def simp38(f, a, b, N=100):
    d = 3
    Nd = proxmultiplo(N, d)/d
    subs = subint(a, b, Nd)
    intervals1 = subs['intervals']
    h = subs['h']/d
    intervals2 = [tuple(subint(i[0], i[1], d)['points']) for i in intervals1]
    I = [(3*h/8)*(f(x[0])+3*f(x[1])+3*f(x[2])+f(x[3])) for x in intervals2]
    Itot = sum(I)
    return Itot


def trap_tab(x, y):
    subsx = subint_tab(x)
    subsy = subint_tab(y)
    h = subsx['h']
    y = subsy['intervals']
    I = [(h[i]/2)*(y[i][0]+y[i][1]) for i in range(len(h))]
    Itot = sum(I)
    return Itot


#FUNÇÃO PARA RESOLUÇÃO DE EQUAÇÕES DIFERENCIAIS ORDINÁRIAS


def euler(f, int, IC, N=1000):
    subs = subint(int[0], int[1], N)
    (x, h, y, n) = (subs['points'], subs['h'], IC, len(IC))
    yout = list()
    for i in range(len(x)):
        yout.append(y)
        y = [y[j] + f(x[i], y)[j]*h for j in range(n)]
    return {'x': x, 'y': yout}


def eulermod(f, int, IC, N=1000):
    subs = subint(int[0], int[1], N)
    (x, h, y, n) = (subs['points'], subs['h'], IC, len(IC))
    x.append(x[-1] + h)
    yout = list()
    for i in range(len(x)-1):
        yout.append(y)
        yEu = [y[j] + f(x[i], y)[j]*h for j in range(n)]
        y = [y[j] + ((f(x[i], y)[j] + f(x[i+1], yEu)[j])/2)*h for j in range(n)]
    del x[-1]
    return {'x': x, 'y': yout}


def eulerpcentral(f, int, IC, N=1000):
    subs = subint(int[0], int[1], N)
    (x, h, y, n) = (subs['points'], subs['h'], IC, len(IC))
    yout = list()
    for i in range(len(x)):
        yout.append(y)
        xm = x[i] + h/2
        ym = [y[j] + f(x[i], y)[j]*(h/2) for j in range(n)]
        y = [y[j] + f(xm, ym)[j]*h for j in range(n)]
    return {'x': x, 'y': yout}


def rk2(f, xi, y, h, n, method='heun'):
    if method == 'standard':
        method = 'heun'
    const = {'eulermod': (1/2, 1/2, 1, 1),
             'eulerpcentral': (0, 1, 1/2, 1/2),
             'heun': (1/4, 3/4, 2/3, 2/3)}
    (c1, c2, a2, b21) = const[method]
    k1 = f(xi, y)
    k2 = f(xi + a2*h, [y[i] + b21*k1[i]*h for i in range(n)])
    y = [y[i] + (c1*k1[i] + c2*k2[i])*h for i in range(n)]
    return y


def rk3(f, xi, y, h, n, method='classico'):
    if method == 'standard':
        method = 'classico'
    const = {'classico': (1/6, 4/6, 1/6, 1/2, 1/2, 1, -1, 2),
             'nystrom': (2/8, 3/8, 3/8, 2/3, 2/3, 2/3, 0, 2/3),
             'quaseotimo': (2/9, 3/9, 4/9, 1/2, 1/2, 3/4, 0, 3/4),
             'heun': (1/4, 0, 3/4, 1/3, 1/3, 2/3, 0, 2/3)}
    (c1, c2, c3, a2, b21, a3, b31, b32) = const[method]
    k1 = f(xi, y)
    k2 = f(xi + a2*h, [y[i] + b21*k1[i]*h for i in range(n)])
    k3 = f(xi + a3*h, [y[i] + b31*k1[i]*h + b32*k2[i]*h for i in range(n)])
    y = [y[i] + (c1*k1[i] + c2*k2[i] + c3*k3[i])*h for i in range(n)]
    return y


def rk4(f, xi, y, h, n, method='classico'):
    if method == 'standard':
        method = 'classico'
    const = {'classico': (1/6, 2/6, 2/6, 1/6, 1/2, 1/2, 1/2, 0, 1/2, 1, 0, 0, 1)}
    (c1, c2, c3, c4, a2, b21, a3, b31, b32, a4, b41, b42, b43) = const[method]
    k1 = f(xi, y)
    k2 = f(xi + a2*h, [y[i] + b21*k1[i]*h for i in range(n)])
    k3 = f(xi + a3*h, [y[i] + b31*k1[i]*h + b32*k2[i]*h for i in range(n)])
    k4 = f(xi + a4*h, [y[i] + b41*k1[i]*h + b42*k2[i]*h + b43*k3[i]*h for i in range(n)])
    y = [y[i] + (c1*k1[i] + c2*k2[i] + c3*k3[i] + c4*k4[i])*h for i in range(n)]
    return y


def rungekutta(f, int, IC, N=1000, rk='rk4', method='standard'):
    rkoptions = {'rk2': rk2, 'rk3': rk3, 'rk4': rk4}
    rkfunc = rkoptions[rk]
    subs = subint(int[0], int[1], N)
    (x, h, y, n) = (subs['points'], subs['h'], IC, len(IC))
    yout = list()
    for i in range(len(x)):
        yout.append(y)
        y = rkfunc(f, x[i], y, h, n, method)
        y = corretor(f, x, y, yout, h, n, i)
    return {'x': x, 'y': yout}

def corretor(f, x, y, yout, h, n, i):
    if len(yout) >= 3 and i < len(x)-1:
        xi = x[:i+1]
        for count in range(10):
            yi = y
            fp1 = f(x[i+1], y)
            fi = f(xi[-1], yout[-1])
            fi_1 = f(xi[-2], yout[-2])
            fi_2 = f(xi[-3], yout[-3])
            y = [yout[-1][k] + (h/24)*(9*fp1[k] + 19*fi[k] - 5*fi_1[k] + fi_2[k]) for k in range(n)]
            erro = sum([abs((y[k] - yi[k])/(yi[k])) for k in range(n)])/n
            if erro <= 1e-15:
                print(count)
                break
    return y


def ab2(f, x, y, h, n):
    fi = f(x[-1], y[-1])
    fi_1 = f(x[-2], y[-2])
    y = [y[-1][i] + (h/2)*(3*fi[i] - fi_1[i])*h for i in range(n)]
    return y


def ab3(f, x, y, h, n):
    fi = f(x[-1], y[-1])
    fi_1 = f(x[-2], y[-2])
    fi_2 = f(x[-3], y[-3])
    y = [y[-1][i] + (h/12)*(23*fi[i] - 16*fi_1[i] + 5*fi_2[i])*h for i in range(n)]
    return y


def ab4(f, x, y, h, n):
    fi = f(x[-1], y[-1])
    fi_1 = f(x[-2], y[-2])
    fi_2 = f(x[-3], y[-3])
    fi_3 = f(x[-4], y[-4])
    y = [y[-1][i] + (h/24)*(55*fi[i] - 59*fi_1[i] + 37*fi_2[i] - 9*fi_3[i])*h for i in range(n)]
    return y
    

def adamsbashforth(f, int, IC, N=1000, ab='ab4', rk='rk4', method='standard'):
    aboptions = {'ab2': (ab2, 2), 'ab3': (ab3, 3), 'ab4': (ab4, 4)}
    (abfunc, nrk) = (aboptions[ab][0], aboptions[ab][1])
    rkoptions = {'rk2': rk2, 'rk3': rk3, 'rk4': rk4}
    rkfunc = rkoptions[rk]
    subs = subint(int[0], int[1], N)
    (x, h, y, n) = (subs['points'], subs['h'], IC, len(IC))
    yout = list()
    for i in range(len(x)):
        yout.append(y)
        if i < nrk-1:
            y = rkfunc(f, x[i], y, h, n, method)
        else:
            y = abfunc(f, x[:i+1], yout, h, n)
    return {'x': x, 'y': yout}


def rk4_octave(f, int, IC, N=1000):
    subs = subint(int[0], int[1], N)
    x = subs['points']
    h = subs['h']
    n = len(IC)
    y = IC
    yout = list()
    for i in range(len(x)):
        yout.append(y)
        k1 = [h*f(x[i], y)[j] for j in range(n)]
        k2 = [h*f(x[i] + h/2, [y[k] + k1[k]/2 for k in range(n)])[j] for j in range(n)]
        k3 = [h*f(x[i] + h/2, [y[k] + k2[k]/2 for k in range(n)])[j] for j in range(n)]
        k4 = [h*f(x[i] + h, [y[k] + k3[k] for k in range(n)])[j] for j in range(n)]
        y = [y[j] + (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])/6 for j in range(n)]
    return {'x': x, 'y': yout}