################## I M P O R T ##################
import numpy as np
import cmath as cm
import scipy as sci
import sympy as sym
import scipy.integrate as integrate
import sympy.utilities as symu

from numpy import sqrt, sin, cos, pi
#from sympy import *

################# F U N C O E S #################
def V_senoidal(wt, Vp):
    return sqrt(2)*Vp*sin(wt)

def V_quadrado(wt, Vp):
    if wt > 2*pi*duty_cycle:
        return 0
    else:
        return Vp

def func_quad(t, f, Vp):
    return f(t, Vp)*f(t, Vp)

################ C A L C U L O S ################
def Vmed(f, T=2*pi, alpha=0, beta=2*pi):
    func = symu.lambdify('x', sym.sympify(f), 'numpy')
    I = integrate.quad(func, alpha, beta)
    return (1/T)*I[0]

def Imed(v, z=1):
    return v/z

def Vrms(f, T=2*pi, alpha=0, beta=2*pi):
    f = '(' + f + ')^2'
    func = symu.lambdify('x', sym.sympify(f), 'numpy')
    I = integrate.quad(func, alpha, beta)
    return sqrt((1/T)*I[0])
