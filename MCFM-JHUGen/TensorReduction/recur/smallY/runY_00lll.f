      subroutine runY_00lll(k,l,Xtwiddle,Gtwiddle,Shat5,N0)
      implicit none
C---  Expression for Eq. 5.60a
C---  Calculates D00lll
C---  Small terms of order Xtwiddle(0,k)*Diiii,Xtwiddle(0,0)*Diiiii
C---  Denominator Gtwiddle(k,l)
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat5(np,z4max,-2:0)

      do ep=-2,0

      Dv(dzziii(z3(l,l,l))+N0,ep)=
     . (Gtwiddle(k,1)*Shat5(1,z4(l,l,l,l),ep)
     . +Gtwiddle(k,2)*Shat5(2,z4(l,l,l,l),ep)
     . +Gtwiddle(k,3)*Shat5(3,z4(l,l,l,l),ep)
     . +Xtwiddle(0,k)*Dv(diiii(z4(l,l,l,l))+N0,ep)
     . -Xtwiddle(0,0)*Dv(diiiii(z5(k,l,l,l,l))+N0,ep))/(8*Gtwiddle(k,l))
      enddo


      return
      end
  



