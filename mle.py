from scipy.optimize import minimize
from scipy.integrate import quad

#Number of bonds in each linear section. $N_i$ in the text
NBONDS   = [25, 12, 18, 19, 14, 9, 4]

#Starting time of each linear section. $T_i$ in the text
CAP_TS   = [0.0, 0.0399, 0.0624, 0.0933, 0.1538, 0.2107, 0.2297]

# Length of time with a constant number of bonds. $t_i$ in the text
LOW_TS   = [0.0243, 0.0197, 0.0250, 0.0527, 0.0346, 0.0111, 0.0312]

# Length of time the bonding or unbonding takes. $\Delta t_i$ in the text
DELTA_TS = [0.0156, 0.0028, 0.0059, 0.0078, 0.0223, 0.0078]


def score(kOn,kOff,fc,v,kBond,kEff):
    """
    Computes the score function (log of the likelihood function) of the experimentally
    observed number of active bonds. 

    kOn:            Bonding rate of a bond in the solid bridge.
    kOff:           Zero force unbonding rate of a bond in the solid bridge.
    fc:             Critical force in the Bell model of force dependent unbonding.
    totNumBonds:    The total number of bonds in the solid bridge.
    v:              Velocity of the AFM cantilever arm.
    kBonds:         Spring constant of the solid bridge bond.
    kEff:           Effective spring constant of the combined interaction of the AFM
                    cantilever and the substrate/microsphere interaction.
    """
    # We assume all the bonds are formed at the begining
    totNumBonds = NS[0]
    def forceBond(_t,_n):
        """
        Computes the force on a single bond as a function of time _t
        assuming _n of the bonds are formed active
        """
        if _n == 0:
            return 0.
        else:
            return kEff*kBond*v*_t/(kEff + kBond*_n)

    def addScore(_t,_t0,_n,_dn):
        """
        Computes the contribution to the score of a single bonding or unbonding event.

        _t:   The length of time between the previous event and the next event
        _t0:  The time the previous event occured
        _n:   The number of active bonds starting at the previous event
        _dn:  Whether the next event is a bonding _dn > 0, or unbonding _dn < 0
        """
        if _dn == 0:
            return 0.0
        onRate = (totNumBonds - _n)*kOn
        offRate = _n*kOff*exp(forceBond(_t + _t0,_n)/fc)
        totRate = (totNumBonds - _n)*kOn*_t + _n*kOff*quad(lambda z: exp(forceBond(_s + _t0,_n)/fc),0.0,_t)[0]
        if _dn > 0:
            return log(onRate) - totRate
        else:
            return log(offRate) - totRate

    S = 0.0
    for i in range(len(CAP_TS) - 1):
        # The last linear section is not used because the epoxy ruptures before
        # more bonds could break.
        Ni = NBONDS[i]
        Ti = CAP_TS[i]
        ti = LOW_TS[i]
        dti= DELTA_TS[i]
        dN = NBONDS[i+1] - Ni
        if dN > 0:
            if dN == 1:
                S += addScore(ti + dti,Ti,Ni,dN)
            else:
                S += addScore(ti,Ti,Ni,dN)
                dN -= 1
                # We assume the bond forming or bond breaking happens at equal intervals 
                dti = dti/dN
                for j in range(dN):
                    S += addScore(dti,Ti + ti + j*dti,Ni + 1+ j,1)
        else:
            if dN == -1:
                S += addScore(ti + dti,Ti,Ni,dN)
            else:
                S += addScore(ti,Ti,Ni,dN)
                dN = abs(dN) - 1
                dti = dti/dN
                for j in range(dN):
                    S += addScore(dti,Ti + ti + j*dti,Ni - 1 - j,-1)
    return S

                
    
def MLE(v,kBond,kEff):
    """
    Uses the score function to compute the MLE kOn,kOff,fc using the
    experimentally observed number of bonds as a function of time.
    
    v:              Velocity of the AFM cantilever arm.
    kBonds:         Spring constant of the solid bridge bond.
    kEff:           Effective spring constant of the combined interaction of the AFM
                    cantilever and the substrate/microsphere interaction.
    """
 
    def optFunc(x):
        val = -score(x[0],x[1],x[2],v,kBond,kEff)
        return val
    # Find the MLE using Nelder Mead optimization method.
    res = minimize(optFunc,x0 = [0.0,0.0,0.0],method = "Nelder-Mead",options = {"disp":True,"maxiter":500},tol = float(1e-6)) 

    return res 

