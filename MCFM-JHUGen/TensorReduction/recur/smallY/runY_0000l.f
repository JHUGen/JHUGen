      subroutine runY_0000l(k,l,Xtwiddle,Gtwiddle,Shat5zz,N0)
      implicit none
C---  Expression for D0000 obtained from 5.50, following the comment after 
C     5.50 on how to add adding additional "00" pairs
C---  (similar to Eq. 5.56a but with "00" added) 
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat5zz(np,z2max,-2:0)

      do ep=-2,0

      Dv(dzzzzi(l)+N0,ep) = -(
     . -(Gtwiddle(k,1)*Shat5zz(1,z2(l,l),ep)
     .  +Gtwiddle(k,2)*Shat5zz(2,z2(l,l),ep)
     .  +Gtwiddle(k,3)*Shat5zz(3,z2(l,l),ep)
     .  +Xtwiddle(0,k)*Dv(dzzii(z2(l,l))+N0,ep)
     .  -Xtwiddle(0,0)*Dv(dzziii(z3(k,l,l))+N0,ep)))/(4d0*Gtwiddle(k,l))

      enddo

      return
      end
  

