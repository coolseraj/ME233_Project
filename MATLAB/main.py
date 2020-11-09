import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import random
import math

def BC(XL, XR, X, PV):
    """
    :param XL: X value for the left boundary of domain
    :param XR: X value for the right boundary of domain
    :param X: Particle's current location
    :param PV: Particle's velocity - we might not need this, but I include it just in case
    :return: Location and velocity of the particle after hitting the BCs
    """

    while X < XL or X > XR:
        if X < XL:
            X = 2*XL - X
            PV = -PV
        if X > XR:
            X = 2*XR - X
            PV = -PV
    return X, PV

def check_year(L, Q, und):
    """

    :param L: particle L
    :param Q: partilce Q
    :param und: array of number of undergrads for each year
    :return: year of particle L and Q
    """
    if L <= und[1]:
        year_L = 1
    elif und[1] < L <= und[2]:
        year_L = 2
    elif und[2] < L <= und[3]:
        year_L = 3
    else:
        year_L = 4

    if Q <= und[1]:
        year_Q = 1
    elif und[1] < Q <= und[2]:
        year_Q = 2
    elif und[2] < L <= und[3]:
        year_Q = 3
    else:
        year_Q = 4
    return year_L, year_Q

def main():
    """
        MNM: max # of molecules
        MNC: max # of cells
        MNSC: max # of sub-cells/collision cells
    """
    und = [100, 150, 250, 94]
    MNM = sum(und)
    MNC = 1
    MNSC = 3
    """
    NM: number of molecules
    PP(MNM): x coordinate of molecule M
    PV(MNM): velocity of molecule M which is set to a constant for all particles 
    IP(MNM): sub-cell number of molecule M
             Input: molecule #
             Output: sub-cell it is located at
    IR(MNM): cross-reference array, molecule number in order of sub-cells
             ex/ the 7th molecule in the order of sub-cells is molecule #98
    """
    NM = 0
    PP = np.zeros(MNM)
    PV = 1.7*np.ones(MNM)
    IP = np.zeros(MNM)
    IR = np.zeros(MNM)

    """
    IC(2,MNC)
        1: start address-1 of molecule # in array IR
        2: # of molecules in the cell
    ISC(MNSC): the cell in which the sub-cell lies
    ISCG(2,MNSC): index info on subcell
    1: start address-1 of molecule
    2: # of molecules in the sub-cell
    """
    IC = np.zeros((2, MNC))
    ISC = np.zeros(MNSC)
    ISCG = np.zeros((2, MNSC))

    DTM = 1e-1

    """
    CW: cell width/whole domain size
    NSC: # of sub-cells per cell
    XL,XR: min and max x-coordinate
    """
    CW = 1
    XL = 0
    XR = CW * MNC + XL
    NSC = MNSC

    infected = np.zeros(MNM)
    # initialize the infected people randomly
    for i in range(0,MNM):
        rand = random.uniform(0, 1)
        if rand > 0.5:
            infected[i] = 1
    """
    infection matrix - we will use a probability instead of this
    what it does is that is says what is the probability that a ith row year 
    infect jth row year - I know, I suck at explaining 
    """
    I_matrix = [[0.5, 0.78, 0.15, 0.36],[0.78, 0.3, 0.98, 0.56],[0.15, 0.98, 0.5, 0.6],[0.36, 0.56, 0.6, 0.1]]

    # set subcells
    for M in range(0, NSC):
        L = M
        ISC[L] = 0

    A = MNM #maximum number of particles in each cell
    N = 1
    for M in range(0, MNM):
        rand = random.uniform(0, 1)
        PP[M] = rand * CW
        IP[M] = int(math.ceil(rand * NSC)) - 1

    c = 0
    tot = 50
    infected_tot = np.zeros(tot)
    while c < tot:
        """
        First section: move particles
        We won't use what I wrote here, instead of velocities, we will use particle's
        PDF to move it from one sub-cell to another
        """

        for N in range(0, MNM):
            MC = 0 #initial cell number
            XI = PP[N]
            DX = PV[N] * DTM
            X = XI + DX
            if X<XL or X>XR:
                X, PV[N] = BC(XL, XR, X, PV[N])
            IP[N] = int(math.ceil(X / CW * NSC)) - 1#new sub-cell #
            PP[N] = X

        """
        Second section: indexing - keeping track of which particle is where
        """
        #molecule numbers are arranged in the order of cells and within the cells, in order of sub-cells
        for N in range(0, MNM):
            MSC = int(IP[N])
            ISCG[1][MSC] += 1
            MC = int(ISC[MSC])
            IC[1][MC] += 1
        M = 0
        for N in range(0, MNC): #start address - 1 has been set for the cells
            IC[0][N] = M
            M += IC[1][N]
        M = 0
        for N in range(0, MNSC): #start address - 1 has been set for sub-cells
            ISCG[0][N] = M
            M += ISCG[1][N]
            ISCG[1][N] = 0
        for N in range(0, MNM):
            MSC = int(IP[N])
            ISCG[1][MSC] += 1
            k = int(ISCG[0][MSC] + ISCG[1][MSC]) - 1
            IR[k] = N + 1

        """
        Third section: Collision or contact between people
        """
        for N in range(0, MNSC): #goes through each sub-cell
            AVN = int(ISCG[1][N])  # number of molecules in a sub-cell
            if AVN > 1:
                ASEL = AVN*(AVN - 1)/2 # number of pairs to be selected in a sub-cell for collision
                for ISEL in range(0, int(ASEL)):
                    kk = int(math.ceil(random.uniform(0, 1) * ISCG[1][N]) + ISCG[0][N]) - 1
                    L = int(IR[kk]) - 1 #molecule L from sub-cell N has been randomely chosen
                    choose = 0
                    while choose == 0:
                        kk = int(math.ceil(random.uniform(0, 1) * ISCG[1][N]) + ISCG[0][N]) - 1
                        Q = int(IR[kk]) - 1#picks molecule Q from the same sub-cell
                        if Q == L:
                            choose = 0 #if molecule Q and L are the same, pick another molecule
                        else:
                            choose = 1
                    if abs(PP[L] - PP[Q]) < 0.001: #collision is accepted - infection happens
                        year_L, year_Q = check_year(L, Q, und)
                        if infected[L] == 1 and infected[Q] == 0: #if L was infected and M was not
                            infected[Q] = round(I_matrix[year_L - 1][year_Q - 1])
                        if infected[L] == 0 and infected[Q] == 1:
                            infected[L] = round(I_matrix[year_L - 1][year_Q - 1])


        infected_tot[c] = sum(infected)
        c += 1
        print(c)
        IC = np.zeros((2, MNC))
        ISCG = np.zeros((2, MNSC))

        plt.plot(c - 1, infected_tot[c - 1], 'bo:')
    print(infected_tot[c - 1])
    plt.show()

if __name__ == "__main__":
    main()

