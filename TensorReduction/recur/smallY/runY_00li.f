      subroutine runY_00li(k,l,i1,Xtwiddle,Gtwiddle,Shat4,N0)
      implicit none
C---  Expression for Eq. 5.58b
C---  Calculates D00li, requires D00ll
C---  Small terms of order Xtwiddle(0,k)*Diii,Xtwiddle(0,0)*Diiii
C---  Denominator Gtwiddle(k,l)
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat4(np,z3max,-2:0)

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      Dv(dzzii(z2(l,i1))+N0,ep)=
     . (-2*Gtwiddle(k,i1)*Dv(dzzii(z2(l,l))+N0,ep)
     . +Gtwiddle(k,1)*Shat4(1,z3(l,l,i1),ep)
     . +Gtwiddle(k,2)*Shat4(2,z3(l,l,i1),ep)
     . +Gtwiddle(k,3)*Shat4(3,z3(l,l,i1),ep)
     . +Xtwiddle(0,k)*Dv(diii(z3(l,l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Dv(diiii(z4(k,l,l,i1))+N0,ep))/(4*Gtwiddle(k,l))
 
      enddo


      return
      end
  



