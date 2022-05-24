import numpy as np
import pykep as pk
import spiceypy as spice

import spice_data as sd
import planetary_data as pd
import numerical_tools as nt

spice.furnsh( sd.leapseconds_kernel )
spice.furnsh( sd.de432s )
spice.furnsh( sd.jup365 )

###### CONFIGURATION #########

T_INIT = '2030-01-01 00:00:00' # ‘YYYY-MM-DD HH:MM:SS’ - Initital time for the search
et_init = spice.str2et( T_INIT )
# print(et_init)

M = 1 # number of spacecraft revolutions on the transfer orbits

##############################


def tot(smaA, smaB, mu):
    """
    The exact time of transfer from A to B is approximated with a multiple of the Hohmann transfer time.
    The assumption is that the exact time of the transfer will be in the interval [ΔT(1−ε ) ΔT(1+ε )] with ε small.
    A: sma of planet A
    smaB: sma of planet B
    mu: gravitational parameter of the central body
    """

    return (2*M+1)*np.pi*np.sqrt((smaA + smaB)**3/(8*mu))



def t_init(bA, bB, k, cb):
    """
    bA: body A (departure)
    bB: body B (arrival)
    k: phasing constraint / number of spacecraft revolutions on the transfer orbits (TBC)
    cb: central body
    """
    state0A = spice.spkgeo(bA['SPICE_ID'], et_init, 'J2000', cb['SPICE_ID'])[0]
    state0B = spice.spkgeo(bB['SPICE_ID'], et_init, 'J2000', cb['SPICE_ID'])[0]

    dTheta0 = nt.vecs2angle(state0A[:3], state0B[:3], deg=False)
    print(dTheta0*nt.r2d)
    deltaT = tot(bA['sma'], bB['sma'], cb['mu']) # seconds
    T_init = ((2*k+1)*np.pi - bB['n']*deltaT - dTheta0)/(bB['n'] - bA['n']) * nt.sec2day # days
    return T_init

if __name__=='__main__':
    print(t_init(pd.earth, pd.mars, 3, pd.sun))

# TODO
# Check if the phasing problem can return negative initial time if the target body has a higher sma than the origin



