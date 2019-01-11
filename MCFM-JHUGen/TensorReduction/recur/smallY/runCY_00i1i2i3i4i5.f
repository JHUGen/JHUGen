      subroutine runCY_00i1i2i3i4i5(k,l,i1,i2,i3,i4,i5,
     . Xtwiddle,Gtwiddle,Shat7,N0)
      implicit none
C---  Expression for extension of Eq. 5.60c
C---  Calculates C00i1i2i3i4i5, requires C00li1i2i3i4
C---  Small terms of order Xtwiddle(0,k)*Ciiiiii,Xtwiddle(0,0)*Ciiiiiii
C---  Denominator Gtwiddle(k,l)
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,k,l,i1,i2,i3,i4,i5,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat7(np,z6max,-2:0)

      if ( (i1.eq.l).or.(i2.eq.l).or.(i3.eq.l).or.(i4.eq.l).or.(i5.eq.l)
     ..or. (i1.eq.0).or.(i2.eq.0).or.(i3.eq.0).or.(i4.eq.0).or.(i5.eq.0)
     .    ) then
      return
      endif

      do ep=-2,0
      Cv(czziiiii(z5(i1,i2,i3,i4,i5))+N0,ep)=
     .(-2d0*Gtwiddle(k,i1)*Cv(czziiiii(z5(l,i2,i3,i4,i5))+N0,ep)
     . -2d0*Gtwiddle(k,i2)*Cv(czziiiii(z5(l,i1,i3,i4,i5))+N0,ep)
     . -2d0*Gtwiddle(k,i3)*Cv(czziiiii(z5(l,i1,i2,i4,i5))+N0,ep)
     . -2d0*Gtwiddle(k,i4)*Cv(czziiiii(z5(l,i1,i2,i3,i5))+N0,ep)
     . -2d0*Gtwiddle(k,i5)*Cv(czziiiii(z5(l,i1,i2,i3,i4))+N0,ep)
     . +Gtwiddle(k,1)*Shat7(1,z6(l,i1,i2,i3,i4,i5),ep)
     . +Gtwiddle(k,2)*Shat7(2,z6(l,i1,i2,i3,i4,i5),ep)
     . +Xtwiddle(k,0)*Cv(ciiiiii(z6(l,i1,i2,i3,i4,i5))+N0,ep)
     . -Xtwiddle(0,0)*Cv(ciiiiiii(z7(k,l,i1,i2,i3,i4,i5))+N0,ep))
     . /(2d0*Gtwiddle(k,l))
 
      enddo

      return
      end
  



