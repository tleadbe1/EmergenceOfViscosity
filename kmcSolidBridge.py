from random import random
from math import exp,log
from scipy.integrate import quad
from scipy.optimize import root_scalar




   
def kmcSolidBridge(tFinal,v,kOff,kOn,fc,nBonds,kEff,kBond):
    """
    Generates a kmc simulation of the stochastic bond model for the solid bridge. 
    tFinal: Simulation runs from time t = 0 until the last event causes  t > tFinal
    v:      Velocity of the AFM cantilever base
    kOff:   The zero force unbonding rate of the solid bridge bonds
    kOn:    The bonding rate of the solid bridge bonds
    fc:     The critical force in the Bell model of bond breaking
    nBonds: The total number of bonds in the solid bridge
    kEff:   The effective spring constant of the AFM arm, and the particle-substrate
            interaction
    kBond:  The spring constant of a bond in the solid bridge

    outputs: ts: The times of each event
             ns: ns[i] is the number of active bonds between time ts[i] and ts[i+1] 
    """

    def forceBond(_t,_n):
        """
        Computes the force on a single bond as a function of time _t, 
        assuming _n bonds are active.
        """
        if _n == 0:
            return 0.
        else:
            return kEff*kBond*v*_t/(kEff + kBond*_n)


    def getNextEvent(_t,_n):
        """
        Randomly select the next event assuming _n bonds are currently active,
        and the current time is _t.
        Returns the time step dt to the next event, and whether the event was 
        a bonding, dn = +1, or an unbonding, dn = -1.
        """
        if _n == 0:
            dt = log(1/random())/nBonds/kOn
            dn = 1
        else:     
            u = random()
            maxT = log(1/u)/((nBonds - _n)*kOn + _n*kOff)
            #use the total rate to chose a time step
            rootFunc = lambda s: (nBonds - _n)*kOn*s + _n*kOff*quad(lambda z: exp(forceBond(_t + z,_n)/fc),0,s)[0] - log(1/u)
            drootFunc = lambda s: (nBonds - _n)*kOn + _n*kOff*exp(forceBond(_t+s,_n)/fc)
            dt = root_scalar(rootFunc, method = "newton",
                        bracket = [0,maxT],
                        x0 = maxT/2,
                        options = {"xtol":float(1e-6),"rtol":float(1e-6)},
                        fprime = drootFunc).root
            # decide whether the next event is a bonding or unbonding
            onRate = (nBonds - _n)*kOn*dt
            offRate = _n*kOff*quad(exp(lambda s: forceBond(_t + s,_n)/fc),0,dt)[0]
            if random() < onRate/(onRate + offRate):
                dn = + 1
            else:
                dn = - 1 
        return dt,dn
 
    t = 0.0
    n = nBonds 
    ts = [t,]
    ns = [n,]
    
    while t < tFinal:
        dt,dn = getNextEvent(t,n)
        t += dt
        n += dn
        ts.append(t)
        ns.append(n)
    return ts,ns

