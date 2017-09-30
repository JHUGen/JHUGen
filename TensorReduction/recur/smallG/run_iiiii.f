      subroutine run_iiiii(j,i1,i2,i3,i4,i5,DetGr,Xtwiddle0,Gtwiddle,
     . Shat6,N0)
      implicit none
C---Fixes Diiiii according to extension of Denner-Dittmaier
c---knowing D00iiii with a correction of order Delta*Diiiiii
      include 'pvDnames.f'
      include 'pvDv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      integer ep,N0,j,i1,i2,i3,i4,i5,n
      double precision DetGr,Xtwiddle0(3),Gtwiddle(3,3)
      double complex Shat6s(3,z5max,-2:0),Shat6(3,z5max,-2:0)
       

      do ep=-2,0
      do n=1,3
      Shat6s(n,z5(i1,i2,i3,i4,i5),ep)=Shat6(n,z5(i1,i2,i3,i4,i5),ep)
     .   -2d0*(delta(n,i1)*Dv(N0+dzziiii(z4(i2,i3,i4,i5)),ep)
     .        +delta(n,i2)*Dv(N0+dzziiii(z4(i1,i3,i4,i5)),ep)
     .        +delta(n,i3)*Dv(N0+dzziiii(z4(i1,i2,i4,i5)),ep)
     .        +delta(n,i4)*Dv(N0+dzziiii(z4(i1,i2,i3,i5)),ep)
     .        +delta(n,i5)*Dv(N0+dzziiii(z4(i1,i2,i3,i4)),ep))
      enddo 

      Dv(diiiii(z5(i1,i2,i3,i4,i5))+N0,ep)=-(
     . +Gtwiddle(j,1)*Shat6s(1,z5(i1,i2,i3,i4,i5),ep)
     . +Gtwiddle(j,2)*Shat6s(2,z5(i1,i2,i3,i4,i5),ep)
     . +Gtwiddle(j,3)*Shat6s(3,z5(i1,i2,i3,i4,i5),ep)
     . -DetGr*Dv(diiiiii(z6(j,i1,i2,i3,i4,i5))+N0,ep))
     . /Xtwiddle0(j)

      enddo

      return
      end
