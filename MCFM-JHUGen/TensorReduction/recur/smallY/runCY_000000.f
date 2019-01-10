      subroutine runCY_000000(k,l,Xtwiddle,Gtwiddle,Shat6zzzz,N0)
      implicit none
C---  Expression for C0000 obtained from 5.50, following the comment after 
C     5.50 on how to add adcing adcitional "00" pairs
C---  Calculates C0000
C---  Small terms of order Xtwiddle(0,k)*C00i,Xtwiddle(0,0)*C00ii
C---  Denominator Gtwiddle(k,l)
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,k,l,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat6zzzz(np,z1max,-2:0)

      do ep=-2,0
         Cv(cc000000+N0,ep) =  
     . -(-(Gtwiddle(k,1)*Shat6zzzz(1,l,ep)
     .    +Gtwiddle(k,2)*Shat6zzzz(2,l,ep)
     .    +Xtwiddle(k,0)*Cv(czzzzi(l)+N0,ep)
     .    -Xtwiddle(0,0)*Cv(czzzzii(z2(k,l))+N0,ep)
     .        ))/(2d0*Gtwiddle(k,l))

      enddo

      return
      end
  
