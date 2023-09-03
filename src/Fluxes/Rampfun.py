import numpy as np
from pylab import * 

       
def rho1(x):
    return x       
       
def rho2(x):
    return x**2 / (1-2*x+2*x**2)

def rho3(x):
    return x**3 / (1 - 3*x + 3*x**2)

def rho4(x):
    return x**4/(2*x**4 - 4*x**3 + 6*x**2 - 4*x + 1)
    
def rhocos(x):
    return 0.5-0.5*cos(pi*x)   # same as in Dritschle's paper 

def rhoexp(x):
    return 1/(1 + exp(-1/(-x + 1.01))*exp(1/(x + 1.0e-2)))   # without overflow
    #return 1./(1. + exp(-1./(-x + 1.0000000001))*exp(1./(x + 1.0e-10)))   #with overflow

def rho_val(sim, x): 
    if sim.rampfun == 'linear':
        return rho1(x)
    elif sim.rampfun == 'quadratic':
        return rho2(x)
    elif sim.rampfun == 'cubic':
        return rho3(x)   
    elif sim.rampfun == 'quartic':
        return rho4(x) 
    elif sim.rampfun == 'cos':
        return rhocos(x)
    elif sim.rampfun == 'exp':
        return rhoexp(x)
