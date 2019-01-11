      subroutine run_iii(j,i1,i2,i3,DetGr,Xtwiddle0,Gtwiddle,
     . Shat4,N0)
      implicit none
C---  Fixes Diii using 5.48
c---  knowing D00ii with a correction of order Delta*Diiii

      include 'pvDnames.f'
      include 'pvDv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      integer ep,N0,j,i1,i2,i3,n
      double precision DetGr,Xtwiddle0(3),Gtwiddle(3,3)
      double complex Shat4(3,z3max,-2:0),Shat4s(3,z3max,-2:0)
       
      do ep=-2,0
      do n=1,3
      Shat4s(n,z3(i1,i2,i3),ep)=Shat4(n,z3(i1,i2,i3),ep)
     .   -2d0*(delta(n,i1)*Dv(N0+dzzii(z2(i2,i3)),ep)
     .        +delta(n,i2)*Dv(N0+dzzii(z2(i1,i3)),ep)
     .        +delta(n,i3)*Dv(N0+dzzii(z2(i1,i2)),ep))
      enddo 

      Dv(diii(z3(i1,i2,i3))+N0,ep)=-(
     . +Gtwiddle(j,1)*Shat4s(1,z3(i1,i2,i3),ep)
     . +Gtwiddle(j,2)*Shat4s(2,z3(i1,i2,i3),ep)
     . +Gtwiddle(j,3)*Shat4s(3,z3(i1,i2,i3),ep)
     . -DetGr*Dv(diiii(z4(j,i1,i2,i3))+N0,ep))
     . /Xtwiddle0(j)
      enddo

      return
      end
