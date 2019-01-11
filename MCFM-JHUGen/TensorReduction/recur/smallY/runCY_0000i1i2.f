      subroutine runCY_0000i1i2(k,l,i1,i2,Xtwiddle,Gtwiddle,Shat6zz,N0)
      implicit none
C---  Expression for Eq. 5.58c
C---  Calculates C0000i1i2, requires C0000li1,C0000li2
C---  Small terms of order Xtwiddle(0,k)*C00iii,Xtwiddle(0,0)*C00iiii
C---  Denominator Gtwiddle(k,l)
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,k,l,i1,i2,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat6zz(np,z3max,-2:0)

      if (  (i1 .eq. l) .or. (i2 .eq. l)
     . .or. (i1 .eq. 0) .or. (i2 .eq. 0)) then
      return
      endif

      do ep=-2,0

      Cv(czzzzii(z2(i1,i2))+N0,ep)=
     .(-2d0*Gtwiddle(k,i1)*Cv(czzzzii(z2(l,i2))+N0,ep)
     . -2d0*Gtwiddle(k,i2)*Cv(czzzzii(z2(l,i1))+N0,ep)
     . +Gtwiddle(k,1)*Shat6zz(1,z3(l,i1,i2),ep)
     . +Gtwiddle(k,2)*Shat6zz(2,z3(l,i1,i2),ep)
     . +Xtwiddle(0,k)*Cv(czziii(z3(l,i1,i2))+N0,ep)
     . -Xtwiddle(0,0)*Cv(czziiii(z4(k,l,i1,i2))+N0,ep)
     . )/(2d0*Gtwiddle(k,l))
      enddo


      return
      end
  



