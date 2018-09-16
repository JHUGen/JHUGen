      subroutine runY_00lllli(k,l,i1,Xtwiddle,Gtwiddle,Shat7,N0)
      implicit none
C---  Expression for extension of Eq. 5.60b
C---  Calculates D00llli, requires D00lllll
C---  Small terms of order Xtwiddle(0,k)*Diiiiii,Xtwiddle(0,0)*Diiiiiii
C---  Denominator Gtwiddle(k,l)
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat7(np,z6max,-2:0)

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      Dv(dzziiiii(z5(l,l,l,l,i1))+N0,ep)=
     . (-2d0*Gtwiddle(k,i1)*Dv(dzziiiii(z5(l,l,l,l,l))+N0,ep)
     . +Gtwiddle(k,1)*Shat7(1,z6(l,l,l,l,l,i1),ep)
     . +Gtwiddle(k,2)*Shat7(2,z6(l,l,l,l,l,i1),ep)
     . +Gtwiddle(k,3)*Shat7(3,z6(l,l,l,l,l,i1),ep)
     . +Xtwiddle(k,0)*Dv(diiiiii(z6(l,l,l,l,l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Dv(diiiiiii(z7(k,l,l,l,l,l,i1))+N0,ep))
     . /(10d0*Gtwiddle(k,l))
 
      enddo


      return
      end
  



