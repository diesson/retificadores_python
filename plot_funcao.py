import numpy as np
from matplotlib import pyplot as plt
from functools import partial

from IPython.display import HTML
from matplotlib import animation
from sympy import *
init_printing()

def plot_fx(f, xlim, ax):
    #fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    #plot_fx(func, xlim=(0, 3*T), ax=ax)
    
    f_ = sympify(f)
    #print('f(x) = {}'.format(f_))
    f = lambdify('x', f_, 'numpy')
    x = np.linspace(xlim[0], xlim[1], 100)
    
    ax.plot(x, f(x), color='red', label=f_)

def plot_tan_fx(f, xs, ax, around=.2):
    #fig, bx = plt.subplots(1, 1, figsize=(6, 6))
    #plot_tan_fx(func, np.arange(0, 3*T, .5), ax=bx, around=.8)
    
    fx = sympify(f).subs('x', 'x0')
    fd = diff(sympify(f)).subs('x', 'x0')
    print('f\'(x0) = {}'.format(fd))
    
    tan_ = sympify('fx + fd * (x - x0)').subs({'fx':fx, 'fd': fd})
    print('tanf(x, x0) = {}'.format(tan_))
    tan = lambdify(('x', 'x0'), tan_, 'numpy')
    
    for x in xs:
        
        x_ = np.linspace(x - around, x + around, 10)
        ax.plot(x_, tan(x_, x))
        
def plot_fx_diodo(func, T, ax, a, b):
    a = a/T
    b = b/T
    
    if a > 0 :
        plot_fx('1.0e-6*x', xlim=(0, a*T), ax=ax)
        
    plot_fx(func, xlim=(a*T, b*T), ax=ax)
    plot_fx('1.0e-6*x', xlim=(b*T, T), ax=ax)
    
    if a > 0 :
        plot_fx('1.0e-6*x', xlim=(T, T + a*T), ax=ax)
        
    plot_fx(func, xlim=(T + a*T, T + b*T), ax=ax)
    plot_fx('1.0e-6*x', xlim=(T + b*T, 2*T), ax=ax)
    
    if a > 0 :
        plot_fx('1.0e-6*x', xlim=(2*T, 2*T + a*T), ax=ax)
        
    plot_fx(func, xlim=(2*T + a*T, 2*T + b*T), ax=ax)
    plot_fx('1.0e-6*x', xlim=(2*T + b*T, 3*T), ax=ax)
    