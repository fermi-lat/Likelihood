      SUBROUTINE dmfit_load(filename)

      IMPLICIT NONE

      INTEGER j,k,l
      integer zn,hasmooth

c-----data tables     
      real phidif
      common/hasim/phidif(-1:250,24,12),hasmooth
c-----lowest mass index
      real*8 milow(12)

      character*(*) filename

        zn=250
c...lowest mass index for channel j
        milow(1)=1   ! c c-bar
        milow(2)=1   ! b b-bar
        milow(3)=1   ! t t-bar
        milow(4)=1   ! tau+ tau-
        milow(5)=1   ! w+ w-
        milow(6)=1   ! z z
        milow(7)=1   ! mu+ mu-
        milow(8)=1   ! gluons
        milow(9)=1   ! gluons
        milow(10)=1   ! gluons
        milow(11)=1   ! gluons
        milow(12)=1   ! gluons

c-----clear the tables
        do j=1,12
          do k=1,24
            do l=0,250
              phidif(l,k,j)=0.0d0
            enddo
          enddo
        enddo

c-----load the table for differential flux
c        open(unit=13,file='gammamc_dif.dat',status='old',
        open(unit=13,file=filename,status='old',
     & form='formatted')

        do j=1,12
          do k=1,24
            if (k.ge.milow(j)) then
                read(13,*) (phidif(l,k,j),l=0,zn-1)
            endif
          enddo
        enddo
      close(13)

        do j=1,12
          do k=1,24
              do l=0,zn-1 ! correct units of dyield / dz
                phidif(l,k,j)=phidif(l,k,j)!/
c     &            (1.0d0/dble(zn))
              enddo

             phidif(-1,k,j)=phidif(0,k,j)
             phidif(zn,k,j)=phidif(zn-1,k,j)

             enddo

         enddo

c         dmfit_load=1

c 2000 format(1000(1x,e12.6))

      return 
      end
