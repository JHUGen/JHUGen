      subroutine run_iiii(j,i1,i2,i3,i4,DetGr,Xtwiddle0,Gtwiddle,
     . Shat5,N0)
C---Fixes Diiii according to extension of Denner-Dittmaier
C---knowing D00iii and dropping terms of order Delta Diiiii
      implicit none
      include 'pvDnames.f'
      include 'pvDv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      integer ep,N0,j,i1,i2,i3,i4,n
      double precision DetGr,Xtwiddle0(3),Gtwiddle(3,3)
      double complex Shat5(3,z4max,-2:0),Shat5s(3,z4max,-2:0)
       
      do ep=-2,0
      do n=1,3
      Shat5s(n,z4(i1,i2,i3,i4),ep)=Shat5(n,z4(i1,i2,i3,i4),ep)
     .   -2d0*(delta(n,i1)*Dv(N0+dzziii(z3(i2,i3,i4)),ep)
     .        +delta(n,i2)*Dv(N0+dzziii(z3(i1,i3,i4)),ep)
     .        +delta(n,i3)*Dv(N0+dzziii(z3(i1,i2,i4)),ep)
     .        +delta(n,i4)*Dv(N0+dzziii(z3(i1,i2,i3)),ep))
      enddo 

      Dv(diiii(z4(i1,i2,i3,i4))+N0,ep)=-(
     . +Gtwiddle(j,1)*Shat5s(1,z4(i1,i2,i3,i4),ep)
     . +Gtwiddle(j,2)*Shat5s(2,z4(i1,i2,i3,i4),ep)
     . +Gtwiddle(j,3)*Shat5s(3,z4(i1,i2,i3,i4),ep)
     . -DetGr*Dv(diiiii(z5(j,i1,i2,i3,i4))+N0,ep))/Xtwiddle0(j)
      enddo

      return
      end
