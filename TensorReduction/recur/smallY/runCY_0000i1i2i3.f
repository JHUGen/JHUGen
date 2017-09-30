      subroutine runCY_0000i1i2i3(k,l,i1,i2,i3,Xtwiddle,Gtwiddle,
     . Shat7zz,N0)
      implicit none
C---  Expression for Eq. 5.60d
C---  Calculates C0000i1i2i3, requires C0000li1i2,C0000li2i3,C0000li3i1
C---  Small terms of order Xtwiddle(0,k)*C00iiii,Xtwiddle(0,0)*C00iiiii
C---  Denominator Gtwiddle(k,l)
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,k,l,i1,i2,i3,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat7zz(np,z4max,-2:0)

      if (  (i1 .eq. l) .or. (i2 .eq. l).or. (i3 .eq. l)
     . .or. (i1 .eq. 0) .or. (i2 .eq. 0).or. (i3 .eq. 0)) then
      return
      endif

      do ep=-2,0
      Cv(czzzziii(z3(i1,i2,i3))+N0,ep)=
     . +(-2d0*Gtwiddle(k,i1)*Cv(czzzziii(z3(l,i2,i3))+N0,ep)
     .   -2d0*Gtwiddle(k,i2)*Cv(czzzziii(z3(l,i1,i3))+N0,ep)
     .   -2d0*Gtwiddle(k,i3)*Cv(czzzziii(z3(l,i1,i2))+N0,ep)
     . +Gtwiddle(k,1)*Shat7zz(1,z4(l,i1,i2,i3),ep)
     . +Gtwiddle(k,2)*Shat7zz(2,z4(l,i1,i2,i3),ep)
     . +Xtwiddle(k,0)*Cv(czziiii(z4(l,i1,i2,i3))+N0,ep)
     . -Xtwiddle(0,0)*Cv(czziiiii(z5(k,l,i1,i2,i3))+N0,ep))
     . /(2d0*Gtwiddle(k,l))
      enddo


      return
      end
  



