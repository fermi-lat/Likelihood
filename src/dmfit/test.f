c-----DMFIT, v 0.1
c     
c     this test program tests the dmfit subroutine and function
c     it also shows an application to KK dark matter

      program test

      implicit none

      character*32 arg
c-----variables needed to run the routines: mass
      real*8 MX,z
c     annihilation channel, smoothing parameter, switch diff/integral (see below)
      integer CH,HSM,SWITCH,ii
c     remember that the following channels, corresponding to WIMP annihilation modes
c     are implemented:
c
c     1: e^+e^-
c     2: mu^+ mu^-
c     3: tau^+ tau^-
c     4: b \bar b
c     5: t \bar t
c     6: gluon gluon
c     7: W+ W-
c     8: Z Z
c     9: c \bar c

c     gamma ray energy, in GeV
      real*8 EE

c-----outputs: 
      real*8 DNDE,DNDM

c-----the actual dmfit functions
      real*8 dmfit_de,dmfit_dm


c-----test case: a 100 GeV WIMP in b\bar b, same smoothing as DarkSUSY's default
      CALL getarg(1,arg)
      read(arg,*) MX
      CALL getarg(2,arg)
      read(arg,*) CH
c      MX=6.d0
c      CH=12
      HSM=0

c-----let's look at an energy of 10 GeV
      EE=0.2d0

      call dmfit_load ("new_gammamc_dif.dat")

c     now let's get the integrated flux

c$$$      DNDE=dmfit_de(MX,CH,EE)
c$$$      DNDM=dmfit_dm(MX,CH,EE)

      write(*,*) '# mass',MX
      write(*,*) '# channel',CH
c      write(*,*) 'energy',EE

c$$$      write(*,*) 
c$$$      write(*,*) 'dN/dE:',DNDE
c$$$      write(*,*) 'd(dN/dE)/dM:',DNDM
      write(*,*) '# idx energy Log(E/MX) EdN/dE d(dN/dE)/dM'
      do ii=0,249
         z=(ii+0.5)/250.
         EE=MX*10.**(10*(z-1.))
         DNDE=dmfit_de(MX,CH,EE)
c          write(*,*) DNDE
c         DNDM=dmfit_dm(MX,CH,EE)
         write(*,*) ii, EE,log10(EE)-log10(MX),DNDE!, DNDM
      enddo

      stop
      end

c-----these are the routines needed to be included:

      include 'dmfit_load.f'
      include 'dmfit_func.f'
      include 'dmfit_comm.f'

