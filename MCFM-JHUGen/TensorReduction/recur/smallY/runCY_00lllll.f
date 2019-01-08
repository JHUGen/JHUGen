      subroutine runCY_00lllll(k,l,Xtwiddle,Gtwiddle,Shat7,N0)
      implicit none
C---  Expression for extension of Eq. 5.60a
C---  Calculates C00llll
C---  Small terms of order Xtwiddle(0,k)*Ciiiiii,Xtwiddle(0,0)*Ciiiiiii
C---  Denominator Gtwiddle(k,l)
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,k,l,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat7(np,z6max,-2:0)

      do ep=-2,0
      Cv(czziiiii(z5(l,l,l,l,l))+N0,ep)=
     . (Gtwiddle(k,1)*Shat7(1,z6(l,l,l,l,l,l),ep)
     . +Gtwiddle(k,2)*Shat7(2,z6(l,l,l,l,l,l),ep)
     . +Xtwiddle(k,0)*Cv(ciiiiii(z6(l,l,l,l,l,l))+N0,ep)
     . -Xtwiddle(0,0)*Cv(ciiiiiii(z7(k,l,l,l,l,l,l))+N0,ep))
     . /(12d0*Gtwiddle(k,l))
      enddo

      return
      end
  



