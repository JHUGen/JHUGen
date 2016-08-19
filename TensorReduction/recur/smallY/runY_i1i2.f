      subroutine runY_i1i2(i,j,i1,i2,f,Xtwiddle,Gtt,Gtwiddle,Shat3,
     . Czero2,N0)
      implicit none
C---  Expression for Eq. 5.59
      include 'TRconstants.f' 
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,i,j,i1,i2,n,m,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np),f(np),
     . Gtt(np,np,np,np)
      double complex Shat3(np,z2max,-2:0),Czero2(z2max,-2:0),
     . bit,pole


      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gtt(i,n,j,m)*f(n)
     . *(Shat3(m,z2(i1,i2),ep)
     . -2d0*delta(m,i1)*Dv(dzzi(i2)+N0,ep)
     . -2d0*delta(m,i2)*Dv(dzzi(i1)+N0,ep))
      enddo
      enddo

      pole=czip
      if (ep .gt. -2) pole=-4*Dv(dzzii(z2(i1,i2))+N0,ep-1)

      Dv(dii(z2(i1,i2))+N0,ep)=
     .  (Gtwiddle(i,j)
     . *(6d0*Dv(dzzii(z2(i1,i2))+N0,ep)+pole-Czero2(z2(i1,i2),ep))
     . +bit+Xtwiddle(0,j)*Dv(diii(z3(i,i1,i2))+N0,ep))/Xtwiddle(i,j)

      enddo
      return
      end
  



