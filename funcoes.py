################## I M P O R T ##################
import importlib
import matplotlib.ticker as tck
import matplotlib.pyplot as plt
import scipy as sci
import scipy.fftpack
import scipy.integrate as integrate
import cmath as cm
import sympy as sym
import sympy.utilities as symu
import numpy as np
from numpy import sqrt, sin, cos, pi
from functools import partial
from math import exp, atan, asin
import sys

################## D E F I N E S ##################
eps = sys.float_info.epsilon

#init_printing()

################ G R A P H I C S ################
def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    def _multiple_formatter(x, pos):
        den = denominator
        num = np.int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            if num==1:
                return r'$%s$'%latex
            elif num==-1:
                return r'$-%s$'%latex
            else:
                return r'$%s%s$'%(num,latex)
        else:
            if num==1:
                return r'$\frac{%s}{%s}$'%(latex,den)
            elif num==-1:
                return r'$\frac{-%s}{%s}$'%(latex,den)
            else:
                return r'$\frac{%s%s}{%s}$'%(num,latex,den)
    return _multiple_formatter

class Multiple:
    def __init__(self, denominator=2, number=np.pi, latex='\pi'):
        self.denominator = denominator
        self.number = number
        self.latex = latex

    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)

    def formatter(self):
        return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))

def plot_fx(f, label, xlim, ax):
    #fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    #plot_fx(func, xlim=(0, 3*T), ax=ax)

    f_ = sym.sympify(f)
    f = sym.lambdify('x', f_, 'numpy')
    x = np.linspace(xlim[0], xlim[1], 500)
    
    g = np.concatenate([f(x), f(x), f(x)])
    x = np.linspace(xlim[0], 3*xlim[1], 1500)
    
    ax.plot(x, g, label=label) #color='red'
    #ax.set_title('')
    ax.set_xlabel('Período ($\omega$t)')
    #ax.set_ylabel(label)
    
    #ax = plt.gca()
    ax.set_aspect('auto')
    ax.axhline(0, color='black', lw=2) #Config eixo
    ax.axvline(0, color='black', lw=2)
    
    ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))

    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    
    ax.set_facecolor('xkcd:white')
    ax.grid(which='major', linestyle='-', linewidth='1.0', color='black')
    ax.grid(which='minor', linestyle='-', linewidth='0.5', color='silver')
    
    ax.grid(True)
    
def plot_fft(f, ax, T=2*np.pi, n_harmonicas = 15, n_pontos = 1000):
    
    harmonicas = calculo_harmonicas(f, T,n_harmonicas, n_pontos)

    x = range(0,harmonicas.size, 1)

    ax.stem(x, harmonicas, '-k')

    plt.title('Resposta em frequência')
    plt.xlabel('Harmônicas')

    ax.axis([0, n_harmonicas, 0, 10 + int(10 * round(float(np.amax(harmonicas))/10))])

    ax.axhline(0, color='black', lw=2)
    ax.axvline(0, color='black', lw=2)

    ax.set_facecolor('xkcd:white')
    ax.yaxis.grid(which='major', linestyle='--', linewidth='1.0', color='silver')

    for a, b in zip( x, harmonicas): 
        if (a >= 0) and (a <= n_harmonicas):
            plt.text(a, b, str(round(b, 3)))
        elif a > n_harmonicas:
            break

def calc_mult(v1, v2, tax):
    m = 1
    teste = 1
    while teste:
        if v1*m < tax*v2:
            m = m*10
        else:
            teste = 0
    return m
    
def plot_info(fv, fi, xlim):
    x = np.linspace(xlim[0], xlim[1], 500)
    
    f_ = sym.sympify(fi)
    f = sym.lambdify('x', f_, 'numpy')
    Ipk = max(f(x))
    
    f_ = sym.sympify(fv)
    f = sym.lambdify('x', f_, 'numpy')
    Vpk = max(f(x))
    
    fp = '(' + fi + ')*(' + fv + ')'
    
    textoV = 'Tensão (V)'
    textoI = 'Corrente (A)'
    
    
    m = calc_mult(Ipk, Vpk, 0.1)
    if m > 1:
        fi = str(m) + '*(' + fi + ')'
        textoI = 'Corrente x'+ str(m) +'(A)'
        
    n = calc_mult(Vpk, Ipk, 0.1)
    if n > 1:
        fv = str(n) + '*(' + fv + ')'
        textoV = 'Tensão x'+ str(n) +'(V)'
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 6), sharey=False, sharex=True)
    
    plot_fx(fv, textoV, xlim, ax=ax1)
    plot_fx(fi, textoI, xlim, ax=ax1)
    ax1.legend()

    #fig, ax2 = plt.subplots(1, 1, figsize=(10, 6), sharey=False, sharex=True)
    plot_fx(fp, "Potencia (W)", xlim, ax=ax2)
    ax2.legend()
    
############### F U N C T I O N S ###############
def calculo_medio(f, T=2*pi, alpha=0, beta=2*pi):
    func = symu.lambdify('x', sym.sympify(f), 'numpy')
    I = integrate.quad(func, alpha, beta)
    return (1/T)*I[0]

def calculo_rms(f, T=2*pi, alpha=0, beta=2*pi):
    f = '(' + f + ')^2'
    func = symu.lambdify('x', sym.sympify(f), 'numpy')
    I = integrate.quad(func, alpha, beta)
    return sqrt((1/T)*I[0])

def calculo_pk(func, ti=0, tf=2*pi):
    f = sym.lambdify('x', sym.sympify(func), 'numpy')
    x = np.linspace(ti, tf, 500)
    return max(f(x))

def rampa_RC(V_pk, teta, wRC, defasagem=0):
    
    return f'{V_pk}*sin({teta})*exp(-(x+{defasagem}-{teta})/({wRC}))'

def calculo_ab(func, valor=0, modo='sympy'):
    
    if modo == 'scipy':
        f = sym.lambdify('x', func, 'numpy')
        ab = sci.optimize.fsolve(f, valor)
        
    elif modo == 'sympy':
        #f = sym.sympify(func)
        eqn = sym.Eq(sym.sympify(func), valor)
        ab = sym.solve(eqn)
    else:
        print("calculo_ab entrada invalida [opções scipy/sympy]")
        
    return ab
  
def calculo_harmonicas(f, T=2*np.pi, n_harmonicas = 15, n_pontos = 1000): 

    f_ = sym.sympify(f)
    f = sym.lambdify('x', f_, 'numpy')

    xf = np.linspace(0, T, n_pontos)
    y = f(xf)  

    v_fft = np.fft.fft(y)

    v_fft = 2.0/n_pontos * np.abs(v_fft)
    v_fft[0] = v_fft[0]/2 
    harmonicas = v_fft[range(0,int((v_fft.size/2)+1), 1)]
    
    return harmonicas[range(0,n_harmonicas+1, 1)]

def calculo_thd(f, T=2*np.pi, n_harmonicas = 500, n_pontos = 1000): 
    
    harmonicas = calculo_harmonicas(f, T, n_harmonicas, n_pontos)
    
    somatorio = harmonicas[0]*harmonicas[0]
    i = 2
    while i < n_harmonicas:
        harmonicas[i] = harmonicas[i]/sqrt(2)
        somatorio += harmonicas[i]*harmonicas[i]
        i += 1
    
    #print(sqrt(somatorio))
    harmonicas[1] = harmonicas[1]/sqrt(2)
    THD = sqrt(somatorio)/harmonicas[1]
    
    return THD

def degrau(a, b, c1 = 1, c2 = -1):
    return f'( {c1}*heaviside(x - {a}, 0) + {c2}*heaviside(x - {b}, 0) + {eps} )'