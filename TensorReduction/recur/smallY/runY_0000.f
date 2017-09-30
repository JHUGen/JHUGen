      subroutine runY_0000(k,l,Xtwiddle,Gtwiddle,Shat4zz,N0)
      implicit none
C---  Expression for D0000 obtained from 5.50, following the comment after 
C     5.50 on how to add adding additional "00" pairs
C---  Calculates D0000
C---  Small terms of order Xtwiddle(0,k)*D00i,Xtwiddle(0,0)*D00ii
C---  Denominator Gtwiddle(k,l)
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat4zz(np,z1max,-2:0)

      do ep=-2,0
         Dv(dd0000+N0,ep) =  
     .    -(-(Gtwiddle(k,1)*Shat4zz(1,l,ep)
     .    +Gtwiddle(k,2)*Shat4zz(2,l,ep)
     .    +Gtwiddle(k,3)*Shat4zz(3,l,ep)
     .    +Xtwiddle(k,0)*Dv(dzzi(l)+N0,ep)
     .    -Xtwiddle(0,0)*Dv(dzzii(z2(k,l))+N0,ep)
     .        ))/(2d0*Gtwiddle(k,l))

      enddo

      return
      end
  
