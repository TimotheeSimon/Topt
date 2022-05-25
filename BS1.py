import numpy as np
import pykep as pk
import spiceypy as spice

import spice_data as sd
import planetary_data as pd
import numerical_tools as nt
import trajectory_sequences as ts

spice.furnsh( sd.leapseconds_kernel )
spice.furnsh( sd.de432s )
spice.furnsh( sd.jup365 )

###### CONFIGURATION #########

T_INIT = '2030-01-01 00:00:00' # ‘YYYY-MM-DD HH:MM:SS’ - Initital time for the search
et_init = spice.str2et( T_INIT )
# print(et_init)

M = 0 # number of spacecraft revolutions on the transfer orbits

##############################

class TransferLeg:
    def __init__(self, bA, bB, cb=pd.sun):
        self.bA = bA
        self.bB = bB
        self.cb = cb
        self.state0A = spice.spkgeo(self.bA['SPICE_ID'], et_init, 'J2000', self.cb['SPICE_ID'])[0]
        self.state0B = spice.spkgeo(bB['SPICE_ID'], et_init, 'J2000', self.cb['SPICE_ID'])[0]
        self.tot = (2*M+1)*np.pi*np.sqrt((self.bA['sma'] + self.bB['sma'])**3/(8*self.cb['mu']))

    def t_init(self,k=0):
        """
        k: phasing constraint / number of spacecraft revolutions on the transfer orbits (TBC)
        """
        dTheta0 = nt.vecs2angle(self.state0A[:3], self.state0B[:3], deg=False)
        # print(dTheta0*nt.r2d)
        deltaT = self.tot # seconds
        T_init = ((2*k+1)*np.pi - self.bB['n']*deltaT - dTheta0)/(self.bB['n'] - self.bA['n']) * nt.sec2day # days
        return T_init

class Trajectory:
    def __init__(self, sequence):
        self.sequence = sequence
        self.Legs_tot = []
        self.Legs_tinit = []
        self.Legs_tlaunch = []

    def calc_tlaunch(self):
        for i in range(len(self.sequence)-1):
            leg = TransferLeg(self.sequence[i], self.sequence[i+1])
            t_init = leg.t_init()
            self.Legs_tinit.append(t_init)
            self.Legs_tlaunch.append(t_init-sum(self.Legs_tot))
            self.Legs_tot.append(leg.tot)


    def F_merit(self, k_list):
        return

if __name__=='__main__':
    #### Test phasing problem single transfer leg ####
    # transfer = TransferLeg(pd.earth, pd.mars, pd.sun)
    # print(transfer.t_init(0))

    #### Test multiple transfer leg ####

    traj = Trajectory(ts.EVMJ)


# TODO
# Check if the phasing problem can return negative initial time if the target body has a higher sma than the origin



