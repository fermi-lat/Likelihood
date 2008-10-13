c********************************************************************
c********************************************************************
*** the smoothing is controlled by the parameter hasmooth in the
*** following manner.
*** hasmooth = 0 - no smoothing
***            1 - smoothing of zi-1,zi and zi+1 bins if z>.3
***            2 - smoothing of zi-2,zi-1,zi,zi+1 and zi+2 if z>0.3
*****************************************************************************

      real*8 function yieldget(zi,mxi,ch)

      implicit none
      integer zi,mxi,ch,fi,ftype,istat,zsmstart,zimin
      integer zn
      integer hasmooth
      real phidif
      common/hasim/phidif(-1:250,18,10),hasmooth

      zn=250

c-----only differential flux

      if(hasmooth.eq.1) then
          if (zi.ge.1.and.zi.le.(zn-1)) then
            yieldget=0.25d0*dble(phidif(zi-1,mxi,ch))
     &        +0.5d0*dble(phidif(zi,mxi,ch))
     &        +0.25d0*dble(phidif(zi+1,mxi,ch))
          elseif (zi.eq.0) then
            yieldget=0.75d0*dble(phidif(zi,mxi,ch))+
     &        0.25d0*dble(phidif(zi+1,mxi,ch))
          elseif (zi.eq.zn) then
            yieldget=0.75d0*dble(phidif(zi,mxi,ch))+
     &        0.25d0*dble(phidif(zi-1,mxi,ch))
          endif

       elseif(hasmooth.eq.2) then
          if (zi.le.(zn-2)) then
            yieldget=0.10d0*dble(phidif(zi-2,mxi,ch))
     &        +0.225d0*dble(phidif(zi-1,mxi,ch))
     &        +0.350d0*dble(phidif(zi,mxi,ch))
     &        +0.225d0*dble(phidif(zi+1,mxi,ch))
     &        +0.100d0*dble(phidif(zi+2,mxi,ch))
          elseif (zi.eq.zn-1) then
            yieldget=0.10d0*dble(phidif(zi-2,mxi,ch))
     &        +0.225d0*dble(phidif(zi-1,mxi,ch))
     &        +0.450d0*dble(phidif(zi,mxi,ch))
     &        +0.225d0*dble(phidif(zi+1,mxi,ch))
          else
            yieldget=0.10d0*dble(phidif(zi-2,mxi,ch))
     &        +0.225d0*dble(phidif(zi-1,mxi,ch))
     &        +0.675d0*dble(phidif(zi,mxi,ch))
          endif
          else
          yieldget=dble(phidif(zi,mxi,ch))
          endif

          return
          end

c********************************************************************
c********************************************************************

      subroutine ifind(value,array,ipl,ii,imin,imax)

      implicit none

      integer imin,imax,i,inew,imint,imaxt,iold,ii
      real*8 ipl,value,array(imin+10:imax+11)    ! +10 to avoid sign problems

      imint=imin+10
      imaxt=imax+10

      if (value.lt.array(imint).or.value.ge.array(imaxt+1)) then
        ii=-5
        return
      endif

      iold=0
      inew=0
      i=(imaxt+imint)/2

 10   if (value.ge.array(i).and.value.lt.array(i+1)) then
        ii=i-10
        ipl=(value-array(i))/(array(i+1)-array(i))
        return
      endif

      if (value.gt.array(i)) then
        inew=(i+1+imaxt)/2
        imint=i
      else
        inew=(imint+i)/2
        imaxt=i
      endif

      i=inew
      if (iold.eq.inew) then
        stop
      endif
      iold=inew

      goto 10

      end

c********************************************************************

      function llg(x,mx,ml)

      real*8 llg,x,mx,ml
      real*8 alpha,pi

      pi=3.141592653589793238d0
      alpha=1.d0/128.d0

      if(x.lt.1.d0-1.d-16) then
      llg=alpha/pi*(x*x-2.d0*x+2.d0)/x*dlog((1.d0-x)*(mx/ml)**2)
      else
      llg=0.d0
      endif

      llg=llg/mx

      return
      end
c********************************************************************
