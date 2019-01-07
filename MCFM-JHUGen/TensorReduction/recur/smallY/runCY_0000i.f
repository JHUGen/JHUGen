      subroutine runCY_0000i(k,l,i1,Xtwiddle,Gtwiddle,Shat5zz,N0)
      implicit none
C---  Expression for C0000i obtained from 5.50, following the comment after 
C     5.50 on how to add adcing adcitional "00" pairs
C---  (similar to Eq. 5.56b but with "00" added) 
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat5zz(np,z2max,-2:0)

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      Cv(czzzzi(i1)+N0,ep)=
     . (-2*Gtwiddle(k,i1)*Cv(czzzzi(l)+N0,ep)
     . +Gtwiddle(k,1)*Shat5zz(1,z2(l,i1),ep)
     . +Gtwiddle(k,2)*Shat5zz(2,z2(l,i1),ep)
     . +Xtwiddle(0,k)*Cv(czzii(z2(l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Cv(czziii(z3(k,l,i1))+N0,ep))/(2*Gtwiddle(k,l))

      enddo

      return
      end
  



