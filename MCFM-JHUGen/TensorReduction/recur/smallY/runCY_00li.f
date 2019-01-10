      subroutine runCY_00li(k,l,i1,Xtwiddle,Gtwiddle,Shat4,N0)
      implicit none
C---  Expression for Eq. 5.58b
C---  Calculates C00li, requires C00ll
C---  Small terms of order Xtwiddle(0,k)*Ciii,Xtwiddle(0,0)*Ciiii
C---  Denominator Gtwiddle(k,l)
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat4(np,z3max,-2:0)

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      Cv(czzii(z2(l,i1))+N0,ep)=
     . (-2*Gtwiddle(k,i1)*Cv(czzii(z2(l,l))+N0,ep)
     . +Gtwiddle(k,1)*Shat4(1,z3(l,l,i1),ep)
     . +Gtwiddle(k,2)*Shat4(2,z3(l,l,i1),ep)
     . +Xtwiddle(0,k)*Cv(ciii(z3(l,l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Cv(ciiii(z4(k,l,l,i1))+N0,ep))/(4*Gtwiddle(k,l))
 
      enddo


      return
      end
  



