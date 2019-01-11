      subroutine runY_0000lli(k,l,i1,Xtwiddle,Gtwiddle,Shat7zz,N0)
      implicit none
C---  Expression for Eq. 5.60b
C---  Calculates D0000lli, requires D0000lll
C---  Small terms of order Xtwiddle(0,k)*D00iiii,Xtwiddle(0,0)*D00iiiii
C---  Denominator Gtwiddle(k,l)
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat7zz(np,z4max,-2:0)

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      Dv(dzzzziii(z3(l,l,i1))+N0,ep)=
     . (-2d0*Gtwiddle(k,i1)*Dv(dzzzziii(z3(l,l,l))+N0,ep)
     . +Gtwiddle(k,1)*Shat7zz(1,z4(l,l,l,i1),ep)
     . +Gtwiddle(k,2)*Shat7zz(2,z4(l,l,l,i1),ep)
     . +Gtwiddle(k,3)*Shat7zz(3,z4(l,l,l,i1),ep)
     . +Xtwiddle(k,0)*Dv(dzziiii(z4(l,l,l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Dv(dzziiiii(z5(k,l,l,l,i1))+N0,ep))
     . /(6d0*Gtwiddle(k,l))
 
      enddo


      return
      end
  



