      subroutine runY_i1i2i3i4i5(i,j,i1,i2,i3,i4,i5,f,Xtwiddle,
     . Gtt,Gtwiddle,Shat6,Czero5,N0)
      implicit none
C---  Expression for Eq. 5.61
      include 'TRconstants.f' 
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,i,j,i1,i2,i3,i4,i5,n,m,np
      parameter(np=3)
      double precision 
     . Xtwiddle(0:np,0:np),Gtwiddle(np,np),f(np),Gtt(np,np,np,np)
      double complex Shat6(np,z5max,-2:0),Czero5(z5max,-2:0),
     . bit,pole


      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gtt(i,n,j,m)*f(n)
     . *(Shat6(m,z5(i1,i2,i3,i4,i5),ep)
     . -2d0*delta(m,i1)*Dv(dzziiii(z4(i2,i3,i4,i5))+N0,ep)
     . -2d0*delta(m,i2)*Dv(dzziiii(z4(i1,i3,i4,i5))+N0,ep)
     . -2d0*delta(m,i3)*Dv(dzziiii(z4(i1,i2,i4,i5))+N0,ep)
     . -2d0*delta(m,i4)*Dv(dzziiii(z4(i1,i2,i3,i5))+N0,ep)
     . -2d0*delta(m,i5)*Dv(dzziiii(z4(i1,i2,i3,i4))+N0,ep))
      enddo
      enddo

      pole=czip
      if (ep .gt. -2) pole=-4d0*Dv(dzziiiii(z5(i1,i2,i3,i4,i5))+N0,ep-1)
      Dv(diiiii(z5(i1,i2,i3,i4,i5))+N0,ep)=
     .  (Gtwiddle(i,j)
     . *(12d0*Dv(dzziiiii(z5(i1,i2,i3,i4,i5))+N0,ep)
     .  +pole-Czero5(z5(i1,i2,i3,i4,i5),ep))
     .  +bit+Xtwiddle(0,j)*Dv(diiiiii(z6(i,i1,i2,i3,i4,i5))+N0,ep)
     .   )/Xtwiddle(i,j)

      enddo
      return
      end
  



