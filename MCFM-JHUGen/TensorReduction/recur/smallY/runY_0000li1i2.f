      subroutine runY_0000li1i2(k,l,i1,i2,Xtwiddle,Gtwiddle,Shat7zz,N0)
      implicit none
C---  Expression for Eq. 5.60c
C---  Calculates D0000li1i2, requires D0000lli1
C---  Small terms of order Xtwiddle(0,k)*D00iiii,Xtwiddle(0,0)*D00iiiii
C---  Denominator Gtwiddle(k,l)
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,i1,i2,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat7zz(np,z4max,-2:0)

      if (  (i1 .eq. l) .or. (i2 .eq. l)
     . .or. (i1 .eq. 0) .or. (i2 .eq. 0)) then
      return
      endif

      do ep=-2,0
      Dv(dzzzziii(z3(l,i1,i2))+N0,ep)=
     .(-2d0*Gtwiddle(k,i1)*Dv(dzzzziii(z3(l,l,i2))+N0,ep)
     . -2d0*Gtwiddle(k,i2)*Dv(dzzzziii(z3(l,l,i1))+N0,ep)
     . +Gtwiddle(k,1)*Shat7zz(1,z4(l,l,i1,i2),ep)
     . +Gtwiddle(k,2)*Shat7zz(2,z4(l,l,i1,i2),ep)
     . +Gtwiddle(k,3)*Shat7zz(3,z4(l,l,i1,i2),ep)
     . +Xtwiddle(0,k)*Dv(dzziiii(z4(l,l,i1,i2))+N0,ep)
     . -Xtwiddle(0,0)*Dv(dzziiiii(z5(k,l,l,i1,i2))+N0,ep))
     . /(4d0*Gtwiddle(k,l))
 
      enddo


      return
      end
  



