      subroutine runCY_00i1i2i3i4(k,l,i1,i2,i3,i4,Xtwiddle,Gtwiddle,
     . Shat6,N0)
      implicit none
C---  Expression for extension of Eq. 5.60c
C---  Calculates C00i1i2i3i4, requires C00li1i2i3
C---  Small terms of order Xtwiddle(0,k)*Ciiiii,Xtwiddle(0,0)*Ciiiiii
C---  Denominator Gtwiddle(k,l)
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,k,l,i1,i2,i3,i4,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat6(np,z5max,-2:0)

      if (  (i1.eq.l) .or. (i2.eq.l) .or. (i3.eq.l) .or. (i4.eq.l)
     . .or. (i1.eq.0) .or. (i2.eq.0) .or. (i3.eq.0) .or. (i4.eq.0)) then
      return
      endif

      do ep=-2,0
      Cv(czziiii(z4(i1,i2,i3,i4))+N0,ep)=
     .(-2d0*Gtwiddle(k,i1)*Cv(czziiii(z4(l,i2,i3,i4))+N0,ep)
     . -2d0*Gtwiddle(k,i2)*Cv(czziiii(z4(l,i1,i3,i4))+N0,ep)
     . -2d0*Gtwiddle(k,i3)*Cv(czziiii(z4(l,i1,i2,i4))+N0,ep)
     . -2d0*Gtwiddle(k,i4)*Cv(czziiii(z4(l,i1,i2,i3))+N0,ep)
     . +Gtwiddle(k,1)*Shat6(1,z5(l,i1,i2,i3,i4),ep)
     . +Gtwiddle(k,2)*Shat6(2,z5(l,i1,i2,i3,i4),ep)
     . +Xtwiddle(k,0)*Cv(ciiiii(z5(l,i1,i2,i3,i4))+N0,ep)
     . -Xtwiddle(0,0)*Cv(ciiiiii(z6(k,l,i1,i2,i3,i4))+N0,ep))
     . /(2d0*Gtwiddle(k,l))
 
      enddo


      return
      end
  



