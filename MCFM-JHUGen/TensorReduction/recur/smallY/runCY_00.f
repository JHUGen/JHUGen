      subroutine runCY_00(k,l,Xtwiddle,Gtwiddle,Shat2,N0)
      implicit none
C---  Expression for Eq. 5.54
C---  Calculates C00
C---  Small terms of order Xtwiddle(0,k)*Ci,Xtwiddle(0,0)*Cii
C---  Denominator Gtwiddle(k,l)
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,k,l,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat2(np,z1max,-2:0)

      do ep=-2,0

      Cv(cc00+N0,ep)=
     .  (Gtwiddle(k,1)*Shat2(1,l,ep)
     .  +Gtwiddle(k,2)*Shat2(2,l,ep)
     .  +Xtwiddle(0,k)*Cv(ci(l)+N0,ep)
     .  -Xtwiddle(0,0)*Cv(cii(z2(k,l))+N0,ep))/(2*Gtwiddle(k,l))
 
c      write(6,*) 'runCY_00:Cv(cc00+N0,ep)',Cv(cc00+N0,ep)

      enddo

      return
      end
  



