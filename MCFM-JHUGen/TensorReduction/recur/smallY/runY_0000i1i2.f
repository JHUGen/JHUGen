      subroutine runY_0000i1i2(k,l,i1,i2,Xtwiddle,Gtwiddle,Shat6zz,N0)
      implicit none
C---  Expression for Eq. 5.58c
C---  Calculates D0000i1i2, requires D0000li1,D0000li2
C---  Small terms of order Xtwiddle(0,k)*D00iii,Xtwiddle(0,0)*D00iiii
C---  Denominator Gtwiddle(k,l)
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,i1,i2,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat6zz(np,z3max,-2:0)

      if (  (i1 .eq. l) .or. (i2 .eq. l)
     . .or. (i1 .eq. 0) .or. (i2 .eq. 0)) then
      return
      endif

      do ep=-2,0

      Dv(dzzzzii(z2(i1,i2))+N0,ep)=
     .(-2d0*Gtwiddle(k,i1)*Dv(dzzzzii(z2(l,i2))+N0,ep)
     . -2d0*Gtwiddle(k,i2)*Dv(dzzzzii(z2(l,i1))+N0,ep)
     . +Gtwiddle(k,1)*Shat6zz(1,z3(l,i1,i2),ep)
     . +Gtwiddle(k,2)*Shat6zz(2,z3(l,i1,i2),ep)
     . +Gtwiddle(k,3)*Shat6zz(3,z3(l,i1,i2),ep)
     . +Xtwiddle(0,k)*Dv(dzziii(z3(l,i1,i2))+N0,ep)
     . -Xtwiddle(0,0)*Dv(dzziiii(z4(k,l,i1,i2))+N0,ep)
     . )/(2d0*Gtwiddle(k,l))
      enddo


      return
      end
  



