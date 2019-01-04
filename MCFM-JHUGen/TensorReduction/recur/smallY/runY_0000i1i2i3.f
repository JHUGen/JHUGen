      subroutine runY_0000i1i2i3(k,l,i1,i2,i3,Xtwiddle,Gtwiddle,
     . Shat7zz,N0)
      implicit none
C---  Expression for Eq. 5.60d
C---  Calculates D0000i1i2i3, requires D0000li1i2,D0000li2i3,D0000li3i1
C---  Small terms of order Xtwiddle(0,k)*D00iiii,Xtwiddle(0,0)*D00iiiii
C---  Denominator Gtwiddle(k,l)
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,i1,i2,i3,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat7zz(np,z4max,-2:0)

      if (  (i1 .eq. l) .or. (i2 .eq. l).or. (i3 .eq. l)
     . .or. (i1 .eq. 0) .or. (i2 .eq. 0).or. (i3 .eq. 0)) then
      return
      endif

      do ep=-2,0
      Dv(dzzzziii(z3(i1,i2,i3))+N0,ep)=
     . +(-2d0*Gtwiddle(k,i1)*Dv(dzzzziii(z3(l,i2,i3))+N0,ep)
     .   -2d0*Gtwiddle(k,i2)*Dv(dzzzziii(z3(l,i1,i3))+N0,ep)
     .   -2d0*Gtwiddle(k,i3)*Dv(dzzzziii(z3(l,i1,i2))+N0,ep)
     . +Gtwiddle(k,1)*Shat7zz(1,z4(l,i1,i2,i3),ep)
     . +Gtwiddle(k,2)*Shat7zz(2,z4(l,i1,i2,i3),ep)
     . +Gtwiddle(k,3)*Shat7zz(3,z4(l,i1,i2,i3),ep)
     . +Xtwiddle(k,0)*Dv(dzziiii(z4(l,i1,i2,i3))+N0,ep)
     . -Xtwiddle(0,0)*Dv(dzziiiii(z5(k,l,i1,i2,i3))+N0,ep))
     . /(2d0*Gtwiddle(k,l))
      enddo


      return
      end
  



