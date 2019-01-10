      subroutine runCY_00lli(k,l,i1,Xtwiddle,Gtwiddle,Shat5,N0)
      implicit none
C---  Expression for Eq. 5.60b
C---  Calculates C00lli, requires C00lll
C---  Small terms of order Xtwiddle(0,k)*Ciiii,Xtwiddle(0,0)*Ciiiii
C---  Denominator Gtwiddle(k,l)
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat5(np,z4max,-2:0)

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      Cv(czziii(z3(l,l,i1))+N0,ep)=
     . (-2*Gtwiddle(k,i1)*Cv(czziii(z3(l,l,l))+N0,ep)
     . +Gtwiddle(k,1)*Shat5(1,z4(l,l,l,i1),ep)
     . +Gtwiddle(k,2)*Shat5(2,z4(l,l,l,i1),ep)
     . +Xtwiddle(0,k)*Cv(ciiii(z4(l,l,l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Cv(ciiiii(z5(k,l,l,l,i1))+N0,ep))
     . /(6*Gtwiddle(k,l))
 
      enddo


      return
      end
  



