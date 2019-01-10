      subroutine runY_00llli(k,l,i1,Xtwiddle,Gtwiddle,Shat6,N0)
      implicit none
C---  Expression for extension of Eq. 5.60b
C---  Calculates D00llli, requires D00llll
C---  Small terms of order Xtwiddle(0,k)*Diiiii,Xtwiddle(0,0)*Diiiiii
C---  Denominator Gtwiddle(k,l)
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat6(np,z5max,-2:0)

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      Dv(dzziiii(z4(l,l,l,i1))+N0,ep)=
     . (-2d0*Gtwiddle(k,i1)*Dv(dzziiii(z4(l,l,l,l))+N0,ep)
     . +Gtwiddle(k,1)*Shat6(1,z5(l,l,l,l,i1),ep)
     . +Gtwiddle(k,2)*Shat6(2,z5(l,l,l,l,i1),ep)
     . +Gtwiddle(k,3)*Shat6(3,z5(l,l,l,l,i1),ep)
     . +Xtwiddle(k,0)*Dv(diiiii(z5(l,l,l,l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Dv(diiiiii(z6(k,l,l,l,l,i1))+N0,ep))
     . /(8d0*Gtwiddle(k,l))
 
      enddo


      return
      end
  



