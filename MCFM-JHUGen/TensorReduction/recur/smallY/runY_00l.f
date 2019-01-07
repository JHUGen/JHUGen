      subroutine runY_00l(k,l,Xtwiddle,Gtwiddle,Shat3,N0)
      implicit none
C---  Expression for Eq. 5.56a
C---  Calculates D00l
C---  Small terms of order Xtwiddle(0,k)*Dii,Xtwiddle(0,0)*Diii
C---  Denominator Gtwiddle(k,l)
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat3(np,z2max,-2:0)

      do ep=-2,0
      Dv(dzzi(l)+N0,ep)=
     . (Gtwiddle(k,1)*Shat3(1,z2(l,l),ep)
     .  +Gtwiddle(k,2)*Shat3(2,z2(l,l),ep)
     .  +Gtwiddle(k,3)*Shat3(3,z2(l,l),ep)
     .  +Xtwiddle(0,k)*Dv(dii(z2(l,l))+N0,ep)
     .  -Xtwiddle(0,0)*Dv(diii(z3(k,l,l))+N0,ep))/(4*Gtwiddle(k,l))
      enddo


      return
      end
  



