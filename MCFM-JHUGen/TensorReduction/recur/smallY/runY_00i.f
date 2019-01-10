      subroutine runY_00i(k,l,i1,Xtwiddle,Gtwiddle,Shat3,N0)
      implicit none
C---  Expression for Eq. 5.56b
C---  Calculates D00i1
C---  Requires D00l
C---  Small terms of order Xtwiddle(0,k)*Dli1,Xtwiddle(0,0)*Dkli1
C---  Denominator Gtwiddle(k,l)
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat3(np,z2max,-2:0)

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      Dv(dzzi(i1)+N0,ep)=
     . (-2*Gtwiddle(k,i1)*Dv(dzzi(l)+N0,ep)
     . +Gtwiddle(k,1)*Shat3(1,z2(l,i1),ep)
     . +Gtwiddle(k,2)*Shat3(2,z2(l,i1),ep)
     . +Gtwiddle(k,3)*Shat3(3,z2(l,i1),ep)
     . +Xtwiddle(0,k)*Dv(dii(z2(l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Dv(diii(z3(k,l,i1))+N0,ep))/(2*Gtwiddle(k,l))

      enddo

      return
      end
  



