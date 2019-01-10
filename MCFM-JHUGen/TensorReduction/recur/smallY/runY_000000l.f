      subroutine runY_000000l(k,l,Xtwiddle,Gtwiddle,Shat7zzzz,N0)
      implicit none
C---  Expression for extension of Eq. 5.60a
C---  Calculates D00llll
C---  Small terms of order Xtwiddle(0,k)*Diiiii,Xtwiddle(0,0)*Diiiiii
C---  Denominator Gtwiddle(k,l)
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat7zzzz(np,z2max,-2:0)

      do ep=-2,0
      Dv(dzzzzzzi(l)+N0,ep)=
     . (Gtwiddle(k,1)*Shat7zzzz(1,z2(l,l),ep)
     . +Gtwiddle(k,2)*Shat7zzzz(2,z2(l,l),ep)
     . +Gtwiddle(k,3)*Shat7zzzz(3,z2(l,l),ep)
     . +Xtwiddle(k,0)*Dv(dzzzzii(z2(l,l))+N0,ep)
     . -Xtwiddle(0,0)*Dv(dzzzziii(z3(k,l,l))+N0,ep))
     . /(4d0*Gtwiddle(k,l))
      enddo

      return
      end
  



