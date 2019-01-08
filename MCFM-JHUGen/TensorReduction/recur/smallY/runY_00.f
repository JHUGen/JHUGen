      subroutine runY_00(k,l,Xtwiddle,Gtwiddle,Shat2,N0)
      implicit none
C---  Expression for Eq. 5.54
C---  Calculates D00
C---  Small terms of order Xtwiddle(0,k)*Di,Xtwiddle(0,0)*Dii
C---  Denominator Gtwiddle(k,l)
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat2(np,z1max,-2:0)

      do ep=-2,0

      Dv(dd00+N0,ep)=
     .  (Gtwiddle(k,1)*Shat2(1,l,ep)
     .  +Gtwiddle(k,2)*Shat2(2,l,ep)
     .  +Gtwiddle(k,3)*Shat2(3,l,ep)
     .  +Xtwiddle(0,k)*Dv(di(l)+N0,ep)
     .  -Xtwiddle(0,0)*Dv(dii(z2(k,l))+N0,ep))/(2*Gtwiddle(k,l))
 
c      write(6,*) 'runY_00:Dv(dd00+N0,ep)',Dv(dd00+N0,ep)

      enddo

      return
      end
  



