from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
import os,sys

dmfit_info={"OLD":("fortran_save2/gammamc_dif.dat",np.array([10, 25, 50, 80.3, 91.2, 100, 150, 176, 200, 250, 350, 500, 750, 1000, 1500, 2000, 3000, 5000]),250/10/np.log(10)),
           "NEW":("../../data/gammamc_dif.dat",np.array([2,4,6,8,10, 25, 50, 80.3, 91.2, 100, 150, 176, 200, 250, 350, 500, 750, 1000, 1500, 2000, 3000, 5000, 7000, 10000]),1.)}

channels={"1":"e+e-","2":"mu+mu-","3":"tau+tau-","4":"b bbar","5":"t tbar","6":"gg","7":"W+W-","8":"ZZ","9":"c cbar","10":"s sbar","11":"u ubar","12":"d dbar"}

def energies(MX):
    z=(np.arange(0,250)+0.5)/250.
    return MX*np.power(10.,10.*(z-1.))

def logE_MX(MX):
    z=(np.arange(0,250)+0.5)/250.
    return 10.*(z-1.)

def readTable(MX, CH,dmfit="NEW"):
    filename,masses,cst=dmfit_info[dmfit]
    gammadiff=np.loadtxt(filename)
    if CH==1:
        chref=9
    elif CH==2:
        chref=7
    elif CH==3:
        chref=4
    elif CH==4:
        chref=2
    elif CH==5:
        chref=3
    elif CH==6:
        chref=8
    elif CH==7:
        chref=5
    elif CH==8:
        chref=6
    elif CH==9:
        chref=1
    elif CH in [10,11,12]:
        chref=CH
    idx=np.where(masses<=MX)[0][-1]
    offset = (chref-1)*len(masses)+idx
    print idx, offset, len(masses)
    return offset,gammadiff[offset]*cst


def myplot(MX,CH):
    #direct read from table
    off,ednde=readTable(MX,CH)
    plt.plot(logE_MX(MX),ednde,"o")
    
    # read through the fortran code
    if os.path.exists("./test.exe"):
        forfile="test_MX%d_CH%d.txt"%(MX,CH)
        if not os.path.exists(forfile):
            os.system("./test.exe %d %d > %s"%(MX,CH,forfile))
        d=np.loadtxt(forfile)
        plt.plot(d[:,2],d[:,3]*energies(MX))

if __name__ == "__main__":
    filename,masses,cst=dmfit_info[sys.argv[1]]
    gammadiff=np.loadtxt(filename)

#MX=12
    CH=int(sys.argv[2])
    for MX in [2,4,6]:
        off,ednde=readTable(MX,CH)
        energies=MX*10**logE_MX(MX)
        if np.any(ednde>0):
            plt.loglog(energies[energies>0.1],ednde[energies>0.1]/energies[energies>0.1],label="MX=%d"%MX)
        
    plt.title("channel : %s"%channels[str(CH)])
    plt.xlabel("E [GeV]")
    plt.ylabel("dN/dE [per annihilation]")
    plt.legend()
    plt.show()
