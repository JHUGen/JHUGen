      subroutine run_0(j,DetGr,Xtwiddle0,Gtwiddle,Shat1,N0)
      implicit none
C---Fixes D0 according to 5.41 Denner-Dittmaier
C---with corrections of order Delta Di
      include 'pvDnames.f'
      include 'pvDv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      integer ep,N0,j
      double precision DetGr,Xtwiddle0(3),Gtwiddle(3,3)
      double complex Shat1(3,-2:0)
       
      do ep=-2,0
      Dv(dd0+N0,ep)=-(
     . +Gtwiddle(j,1)*Shat1(1,ep)
     . +Gtwiddle(j,2)*Shat1(2,ep)
     . +Gtwiddle(j,3)*Shat1(3,ep)
     . -DetGr*Dv(di(j)+N0,ep))/Xtwiddle0(j)
      enddo

      return
      end
