      subroutine run_ii(j,i1,i2,DetGr,Xtwiddle0,Gtwiddle,Shat3,N0)
      implicit none
C---  Fixes Dii using Eq. 5.45
C     knowing D00i with correction of order Delta Diii
      include 'pvDnames.f'
      include 'pvDv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      integer ep,N0,j,i1,i2,n
      double precision DetGr,Xtwiddle0(3),Gtwiddle(3,3)
      double complex Shat3(3,z2max,-2:0),Shat3s(3,z2max,-2:0)
       
      do ep=-2,0
      do n=1,3
      Shat3s(n,z2(i1,i2),ep)=Shat3(n,z2(i1,i2),ep)
     .   -2d0*(delta(n,i1)*Dv(N0+dzzi(i2),ep)
     .        +delta(n,i2)*Dv(N0+dzzi(i1),ep))
      enddo 
         
      Dv(dii(z2(i1,i2))+N0,ep)=-(
     . +Gtwiddle(j,1)*Shat3s(1,z2(i1,i2),ep)
     . +Gtwiddle(j,2)*Shat3s(2,z2(i1,i2),ep)
     . +Gtwiddle(j,3)*Shat3s(3,z2(i1,i2),ep)
     . -DetGr*Dv(diii(z3(j,i1,i2))+N0,ep))/Xtwiddle0(j)

      enddo

      return
      end
