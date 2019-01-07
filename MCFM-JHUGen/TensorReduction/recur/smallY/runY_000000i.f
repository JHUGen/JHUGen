      subroutine runY_000000i(k,l,i1,Xtwiddle,Gtwiddle,Shat7zzzz,N0)
      implicit none
C---  Expression for D000000i obtained from 5.50, following the comment after 
C     5.50 on how to add adding additional "00" pairs
C---  (similar to Eq. 5.56b but with "0000" added) 
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat7zzzz(np,z2max,-2:0)

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      Dv(dzzzzzzi(i1)+N0,ep)=
     . (-2d0*Gtwiddle(k,i1)*Dv(dzzzzzzi(l)+N0,ep)
     . +Gtwiddle(k,1)*Shat7zzzz(1,z2(l,i1),ep)
     . +Gtwiddle(k,2)*Shat7zzzz(2,z2(l,i1),ep)
     . +Gtwiddle(k,3)*Shat7zzzz(3,z2(l,i1),ep)
     . +Xtwiddle(0,k)*Dv(dzzzzii(z2(l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Dv(dzzzziii(z3(k,l,i1))+N0,ep)
     . )/(2d0*Gtwiddle(k,l))

      enddo

      return
      end
  



