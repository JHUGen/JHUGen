      subroutine runCY_00lli1i2(k,l,i1,i2,Xtwiddle,Gtwiddle,Shat6,N0)
      implicit none
C---  Expression for extension of Eq. 5.60c
C---  Calculates C00lli1i2, requires C00llli1
C---  Small terms of order Xtwiddle(0,k)*Ciiiii,Xtwiddle(0,0)*Ciiiiii
C---  Denominator Gtwiddle(k,l)
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,k,l,i1,i2,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat6(np,z5max,-2:0)

      if (  (i1 .eq. l) .or. (i2 .eq. l)
     . .or. (i1 .eq. 0) .or. (i2 .eq. 0)) then
      return
      endif

      do ep=-2,0
      Cv(czziiii(z4(l,l,i1,i2))+N0,ep)=
     .(-2d0*Gtwiddle(k,i1)*Cv(czziiii(z4(l,l,l,i2))+N0,ep)
     . -2d0*Gtwiddle(k,i2)*Cv(czziiii(z4(l,l,l,i1))+N0,ep)
     . +Gtwiddle(k,1)*Shat6(1,z5(l,l,l,i1,i2),ep)
     . +Gtwiddle(k,2)*Shat6(2,z5(l,l,l,i1,i2),ep)
     . +Xtwiddle(k,0)*Cv(ciiiii(z5(l,l,l,i1,i2))+N0,ep)
     . -Xtwiddle(0,0)*Cv(ciiiiii(z6(k,l,l,l,i1,i2))+N0,ep))
     . /(6d0*Gtwiddle(k,l))
 
      enddo


      return
      end
  



