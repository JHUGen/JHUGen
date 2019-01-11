      subroutine runY_0000ll(k,l,Xtwiddle,Gtwiddle,Shat6zz,N0)
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
      double complex Shat6zz(np,z3max,-2:0)

      do ep=-2,0
      Dv(dzzzzii(z2(l,l))+N0,ep)=
     . (Gtwiddle(k,1)*Shat6zz(1,z3(l,l,l),ep)
     . +Gtwiddle(k,2)*Shat6zz(2,z3(l,l,l),ep)
     . +Gtwiddle(k,3)*Shat6zz(3,z3(l,l,l),ep)
     . +Xtwiddle(k,0)*Dv(dzziii(z3(l,l,l))+N0,ep)
     . -Xtwiddle(0,0)*Dv(dzziiii(z4(k,l,l,l))+N0,ep))
     . /(6d0*Gtwiddle(k,l))
      enddo

      return
      end
  



