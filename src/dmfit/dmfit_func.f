c-----DMFIT dN/dE routine 

      REAL*8 FUNCTION dmfit_de(MX,CH,EE)

      IMPLICIT NONE

c-----ARGUMENTS:

c-----WIMP mass
c     since the MC simulations run between 10 and 5000 GeV, we impose 
c     the corresponding cuts
      REAL*8 MX

c-----DM annihilation mode
      INTEGER CH

c-----the following channels, corresponding to WIMP annihilation modes
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
c     10: cosmo b \bar b
c     11: cosmo gam gam Line
c
c     Notice that if a channel is not kinematically available
c     (e.g., ch 7 with MX < M_W), then automatically we fit
c     with default channel 4. This applies to CH=5,7,8

c-----Smoothing procedure (hard-wired in this routine!)
      INTEGER HSM
c     HSM = 0: simply takes the closest point in the table
c     HSM = 1: interpolates between 3 points
c     HSM = 2: interpolates between 5 points

c-----User-supplied photon energy, in GeV
      REAL*8 EE

c-----Photon flux (differential)
      REAL*8 FLUX

c-----whether or not to smooth the interpolation
      integer hasmooth
c-----masses of top, W, Z
      real*8 mt,mw,mz
c-----lowest mass index
c      real*8 milow(10)

c-----this is the reference channel
      integer chref

c-----variables to initialize tables
      integer i
c,ntype
c      integer j,k,l
      integer zi,m1i,m2i,zn
      real*8 z,ndec,zpl,tmp
      real*8 mi(24),mp1,mp2
      real*8 phi1,phi2
      real*8 yieldget
      real*8 lge

c-----number of decades tabulated
      parameter(ndec=10.0d0)  
      parameter(lge=0.434294481903d0)

      real*8 zindex(-1:250,2)

c-----backup masses and energies for low energy extrapolation
      real*8 EEold,MXold
c-----function that computes the differential \gamma flux from e+e-
      real*8 llg

c-----data tables     
      real phidif
      common/hasim/phidif(-1:250,24,12),hasmooth


********************************************************************
c     backup the input energy and mass values
         EEold=EE
         MXold=MX

c-----mass cut: lower limit, all channels but not e+e-
c$$$      if(MX.lt.10.d0.and.CH.ne.1) then
c$$$c         write(*,*) 'WARNING: MASS IS TOO LOW!'
c$$$c         write(*,*) 'USING EXTRAPOLATION OF MC DATA!'
c$$$c-----here we do the barbaric extrapolation to lower masses
c$$$         EE=EE/MX*10.d0
c$$$         MX=10.d0
c$$$c         stop
c$$$      endif

      HSM=0

c-----imposes the smoothing to use in the interpolation
      if (hsm.eq.0) then
      hasmooth=0
      elseif(hsm.eq.1) then
      hasmooth=1
      elseif(hsm.eq.2) then
      hasmooth=2
      else
      hasmooth=0
      endif

c-----set the masses for top, W, Z    
      mt=174.3d0
      mw=80.33d0
      mz=91.187d0

c-----switches to CH 4 if mass limit inconsistent
      if(CH.eq.5) then
         if(MX.lt.mt) then
            CH=4
            write(*,*) 'ERROR: CHANNEL NOT KINEMATICALLY ALLOWED!'
            write(*,*) '       SWITHCING TO DEFAULT, B\bar B!'
         endif
      endif

      if(CH.eq.7) then
         if(MX.lt.mw) then
            CH=4
            write(*,*) 'ERROR: CHANNEL NOT KINEMATICALLY ALLOWED!'
            write(*,*) '       SWITHCING TO DEFAULT, B\bar B!'
         endif
      endif

      if(CH.eq.8) then
         if(MX.lt.mz) then
            CH=4
            write(*,*) 'ERROR: CHANNEL NOT KINEMATICALLY ALLOWED!'
            write(*,*) '       SWITHCING TO DEFAULT, B\bar B!'
         endif
      endif

c-----translate the channels

c-----this is the default channel, b \bar b
      chref=4

c$$$      if(CH.eq.1) then
c$$$c-----for the e+e- channel, go ahead and compute it!
c$$$         dmfit_de=llg(EE/MX,MX,0.511d-3)
c$$$         return
c$$$
c$$$c-----for muon channel at low mass, the MonteCarlo statistics 
c$$$c-----resulting in the DMFIT table is too low, so we switch to
c$$$c-----the same equation as the electron channel instead.
c$$$      elseif(CH.eq.2.and.MX.lt.10.d0) then
c$$$         dmfit_de=llg(EE/MX,MX,0.1057d0)
c$$$         return
c$$$      elseif(CH.eq.2) then
c$$$         chref=7
c$$$      elseif(CH.eq.3) then
c$$$         chref=4
c$$$      elseif(CH.eq.4) then
c$$$         chref=2
c$$$      elseif(CH.eq.5) then
c$$$         chref=3
c$$$      elseif(CH.eq.6) then
c$$$         chref=8
c$$$      elseif(CH.eq.7) then
c$$$         chref=5
c$$$      elseif(CH.eq.8) then
c$$$         chref=6
c$$$      elseif(CH.eq.9) then
c$$$         chref=1
c$$$      elseif(CH.eq.10) then
c$$$         chref=9
c$$$      elseif(CH.eq.11) then
c$$$         chref=10
c$$$      else
c$$$         chref=4
c$$$      endif
      if(CH.eq.1) then
c-----e+e-
         chref=9
      elseif(CH.eq.2) then
c-----mu+mu-
         chref=7
      elseif(CH.eq.3) then
c-----tau+tau-
         chref=4
      elseif(CH.eq.4) then
c-----b bbar
         chref=2
      elseif(CH.eq.5) then
c------t tbar
         chref=3
      elseif(CH.eq.6) then
c-----g g
         chref=8
      elseif(CH.eq.7) then
c-----W+W-
         chref=5
      elseif(CH.eq.8) then
c-----ZZ
         chref=6
      elseif(CH.eq.9) then
c-----mu+mu-
         chref=1
      elseif(CH.eq.10) then
c-----s sbar
         chref=10
      elseif(CH.eq.11) then
c-----u ubar
         chref=11
      elseif(CH.eq.12) then
c-----d dbar
         chref=12
      else
         chref=4
      endif

c********************************************************************
c-----this initializes and loads the MC simulation data for gammas

c-----masses for simulation corresponding to mass index i
c$$$        mi(1)=10.0
c$$$        mi(2)=25.0
c$$$        mi(3)=50.0
c$$$        mi(4)=80.3
c$$$        mi(5)=91.2
c$$$        mi(6)=100.0
c$$$        mi(7)=150.0
c$$$        mi(8)=176.0
c$$$        mi(9)=200.0
c$$$        mi(10)=250.0
c$$$        mi(11)=350.0
c$$$        mi(12)=500.0
c$$$        mi(13)=750.0
c$$$        mi(14)=1000.0
c$$$        mi(15)=1500.0
c$$$        mi(16)=2000.0
c$$$        mi(17)=3000.0
c$$$        mi(18)=5000.0
        mi(1)=2.d0
        mi(2)=4.d0
        mi(3)=6.d0
        mi(4)=8.d0
        mi(5)=10.d0
        mi(6)=25.d0
        mi(7)=50.d0
        mi(8)=80.3d0
        mi(9)=91.2d0
        mi(10)=100.d0
        mi(11)=150.d0
        mi(12)=176.d0
        mi(13)=200.d0
        mi(14)=250.d0
        mi(15)=350.d0
        mi(16)=500.d0
        mi(17)=750.d0
        mi(18)=1000.d0
        mi(19)=1500.d0
        mi(20)=2000.d0
        mi(21)=3000.d0
        mi(22)=5000.d0
        mi(23)=7000.d0
        mi(24)=10000.d0

c-----initialize eindex array where the energies for the bins are stored
c-----integrated yields (lower end of each bin)
        zn=250
        do i=-1,zn
          zindex(i,1)=dble(i)/dble(zn)
        enddo
        do i=-1,zn
          zindex(i,2)=dble(i)/dble(zn)+0.5d0/dble(zn)
        enddo

c********************************************************************

c-----interpolation: lower energy

      if(EE.ge.MX) then
         dmfit_de=0.d0
         return
      endif

          z=(log10(EE/MX)+ndec)/ndec

        call ifind(z,zindex(-1,2),zpl,zi,-1,zn-1)

        if (zi.eq.-5.or.zi.ge.zn) then
         dmfit_de=0.d0
          return
        endif

        call ifind(MX,mi(1),tmp,m1i,1,23)
        mp1=mi(m1i)
        m2i=m1i+1
        mp2=mi(m2i)

        if (MX.ge.mi(24)) then
          m1i=24
          m2i=24
          mp1=mi(24)
          mp2=mp1

          FLUX =
     &      (1.0-zpl)*yieldget(zi,m1i,chref)+
     &      zpl*yieldget(zi+1,m1i,chref)

c$$$          FLUX=FLUX*lge/(ndec*EE)
          FLUX=FLUX/EE
        else
c           write(*,*) CH,zi,m1i,EE,yieldget(zi,m1i,chref)
          phi1 =
     &      (1.0-zpl)*yieldget(zi,m1i,chref)+
     &      zpl*yieldget(zi+1,m1i,chref)
          phi2 =
     &      (1.0-zpl)*yieldget(zi,m2i,chref)+
     &      zpl*yieldget(zi+1,m2i,chref)
          FLUX = phi1 + (phi2-phi1)*(MX-mp1)*(MX+mp1)/
     &      ((mp2-mp1)*(mp2+mp1))

c$$$          FLUX=FLUX*lge/(ndec*EE)
          FLUX=FLUX/EE
        endif

 105    continue

c$$$        if(MXOLD.lt.10.d0.and.CH.ne.1) then
c$$$           FLUX=FLUX/(MXOLD/10.d0)
c$$$        endif


        dmfit_de=FLUX

c-----here we restore the input values if they were altered
         EE=EEold
         MX=MXold


        return

        end



c-----DMFIT d(dN/dE)/dM routine 

      REAL*8 FUNCTION dmfit_dm(MX,CH,EE)

      IMPLICIT NONE

c-----ARGUMENTS:

      REAL*8 MX,EE
      INTEGER CH
      
c-----the other function...

      REAL*8 dmfit_de

c-----auxiliary variables

      REAL*8 delta,m1,m2,d1,d2,inc1,inc2

********************************************************************
c     this sets the "delta" for the computation of the derivative

      delta=0.0001d0

      m1=MX/(1.d0+delta)
      d1=MX-m1

      m2=MX*(1.d0+delta)
      d2=m2-MX

      if(m1.gt.EE) then
      inc1=dmfit_de(MX,CH,EE)-dmfit_de(m1,CH,EE)
      inc1=inc1/d1

      inc2=dmfit_de(m2,CH,EE)-dmfit_de(MX,CH,EE)
      inc2=inc2/d2

      dmfit_dm=(inc1+inc2)/2.d0
      else
      inc2=dmfit_de(m2,CH,EE)-dmfit_de(MX,CH,EE)
      inc2=inc2/d2

      dmfit_dm=inc2
      endif

      return
      end

c********************************************************************
c-----DMFIT \int_E^MX (dN/dE) dE

      REAL*8 FUNCTION dmfit_deint(MX,CH,EE)

      IMPLICIT NONE

c-----ARGUMENTS:

      REAL*8 MX,EE
      INTEGER CH
      
c-----the other function, which we want to integrate...

      REAL*8 dmfit_de

c-----auxiliary variables

      REAL*8 delta,integral,ene1,ene2,dnde1,dnde2
      INTEGER nstep,ii


c     sets to 0 the integral
      integral=0.d0

c     defines the number of steps - say 100
      nstep=100

c     case where energy is above the mass - flux is 0
      if (EE.ge.MX) then
         dmfit_deint=0.d0
         return
      endif

c     integration
      do ii=0,nstep

         ene1=EE*dexp((dlog(MX)-dlog(EE))*dble(ii)/dble(nstep))
         ene2=EE*dexp((dlog(MX)-dlog(EE))*dble(ii+1)/dble(nstep))

         delta=ene2-ene1

         dnde1=dmfit_de(MX,CH,ene1)
         dnde2=dmfit_de(MX,CH,ene2)

         integral=integral+(dnde1+dnde2)*delta/2.d0

       enddo

       dmfit_deint=integral

       return
       end


c********************************************************************
c-----DMFIT \int_E^MX (dN/dM) dE

      REAL*8 FUNCTION dmfit_dmint(MX,CH,EE)

      IMPLICIT NONE

c-----ARGUMENTS:

      REAL*8 MX,EE
      INTEGER CH
      
c-----the other function, which we want to integrate...

      REAL*8 dmfit_dm

c-----auxiliary variables

      REAL*8 delta,integral,ene1,ene2,dnde1,dnde2
      INTEGER nstep,ii


c     sets to 0 the integral
      integral=0.d0

c     defines the number of steps - say 100
      nstep=100

c     case where energy is above the mass - flux is 0
      if (EE.ge.MX) then
         dmfit_dmint=0.d0
         return
      endif

c     integration
      do ii=0,nstep

         ene1=EE*dexp((dlog(MX)-dlog(EE))*dble(ii)/dble(nstep))
         ene2=EE*dexp((dlog(MX)-dlog(EE))*dble(ii+1)/dble(nstep))

         delta=ene2-ene1

         dnde1=dmfit_dm(MX,CH,ene1)
         dnde2=dmfit_dm(MX,CH,ene2)

         integral=integral+(dnde1+dnde2)*delta/2.d0

       enddo

       dmfit_dmint=integral

       return
       end
