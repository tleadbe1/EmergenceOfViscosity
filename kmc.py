from random import random               #Generates uniform random numbers on (0,1)
from math import exp,log                #exponetial and natural logarithm functions
from scipy.integrate import quad        #General Purpose Numerical Integration
from scipy.optimize import root_scalar  #Method to find roots of a scalar valued function



def kmcForceInd(tFinal,kOff,kOn):
    """
    Generates the bonding and unbonding trajectory of a single bond
    assuming that the rate of unbonding is force independent.

    tFinal: stop simulation once the total time t > tFinal
    kOff:   the rate of unbonding
    kOn:    the rate of bonding
    """
    ts= []
    t = 0.0
    while t < tFinal:
        tOff = -log(random())/kOn 
        tOn = -log(random())/kOff
        ts.extend([t + tOff,t + tOff + tOn]) 
        t += tOff + tOn
    return ts


def kmc(tFinal,v,kOff,kOn,force):
    """
    Generates the bonding and unbonding trajectory of a single bond
    assuming that the rate of unbonding is force dependent.

    tFinal: stop simulation once the total time t > tFinal
    v:      The velocity of the top plate in the stochastic bond model
    kOff:   The zero force rate of unbonding
    kOn:    The rate of bonding
    force:  Magnitude of the force on the bond as a function of bond length.
                Assumes the force is normalize, ie force = a\beta F   
    """
    def getTOn():
        """
        Generates a random time spent bonded from the force dependent unbonding 
        distribution assuming the Bell model.
        """
        rate = lambda s: kOff*exp(force(v*s))
        u = random()
        maxT = log(1/u)/kOff
        rootFunc=lambda s: quad(rate,0.0,s)[0] - log(1/u)
        dt = root_scalar(rootFunc, method = "newton",
                        bracket = [0,maxT],
                        x0 = maxT/2,
                        options = {"xtol":float(1e-6),"rtol":float(1e-6)},
                        fprime = rate).root
        return dt
    ts = []
    t = 0.0 
    while t < tFinal:
        tOff = -log(random())/kOn
        tOn = getTOn()
        ts.extend([t + tOff,t + tOff + tOn])
        t += tOn + tOff
    return ts


