from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
import os,sys

dmfit_info={"OLD":("fortran_save2/gammamc_dif.dat",np.array([10, 25, 50, 80.3, 91.2, 100, 150, 176, 200, 250, 350, 500, 750, 1000, 1500, 2000, 3000, 5000]),250/10/np.log(10)),
           "NEW":("gammamc_dif.dat",np.array([2,4,6,8,10, 25, 50, 80.3, 91.2, 100, 150, 176, 200, 250, 350, 500, 750, 1000, 1500, 2000, 3000, 5000, 7000, 10000]),1.)}

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
    if CH==4:
        chref=2
    elif CH in [10,11,12]:
        chref=CH
    idx=np.where(masses<=MX)[0][-1]
    print idx, len(masses)
    offset = (chref-1)*len(masses)+idx
    return offset,gammadiff[offset]*cst


def plot(MX,CH):
    #direct read from table
    off,ednde=readTable(MX,CH)
    plt.plot(logE_MX(MX),ednde,"o")
    
    # read through the fortran code
    forfile="test_MX%d_CH%d.txt"%(MX,CH)
    if not os.path.exists(forfile):
        os.system("./test.exe %d %d > %s"%(MX,CH,forfile))
    d=np.loadtxt(forfile)
    plt.plot(d[:,2],d[:,3]*energies(MX))

if __name__ == "__main__":
    filename,masses,cst=dmfit_info[sys.argv[1]]
    gammadiff=np.loadtxt(filename)

#MX=12
    CH=4
    for MX in masses:
        plot(MX,CH)
        
    plt.show()
