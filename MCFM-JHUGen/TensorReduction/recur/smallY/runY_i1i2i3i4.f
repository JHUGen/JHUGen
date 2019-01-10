      subroutine runY_i1i2i3i4(i,j,i1,i2,i3,i4,f,Xtwiddle,Gtt,Gtwiddle,
     . Shat5,Czero4,N0)
      implicit none
C---  Expression for Eq. 5.61
      include 'TRconstants.f' 
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,i,j,i1,i2,i3,i4,n,m,np
      parameter(np=3)
      double precision 
     . Xtwiddle(0:np,0:np),Gtwiddle(np,np),f(np),Gtt(np,np,np,np)
      double complex Shat5(np,z4max,-2:0),Czero4(z4max,-2:0),
     . bit,pole


      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gtt(i,n,j,m)*f(n)
     . *(Shat5(m,z4(i1,i2,i3,i4),ep)
     . -2d0*delta(m,i1)*Dv(dzziii(z3(i2,i3,i4))+N0,ep)
     . -2d0*delta(m,i2)*Dv(dzziii(z3(i1,i3,i4))+N0,ep)
     . -2d0*delta(m,i3)*Dv(dzziii(z3(i1,i2,i4))+N0,ep)
     . -2d0*delta(m,i4)*Dv(dzziii(z3(i1,i2,i3))+N0,ep))
      enddo
      enddo

      pole=czip
      if (ep .gt. -2) pole=-4d0*Dv(dzziiii(z4(i1,i2,i3,i4))+N0,ep-1)
      Dv(diiii(z4(i1,i2,i3,i4))+N0,ep)=
     .  (Gtwiddle(i,j)
     . *(10d0*Dv(dzziiii(z4(i1,i2,i3,i4))+N0,ep)
     .  +pole-Czero4(z4(i1,i2,i3,i4),ep))
     .  +bit+Xtwiddle(0,j)*Dv(diiiii(z5(i,i1,i2,i3,i4))+N0,ep)
     .   )/Xtwiddle(i,j)

      enddo
      return
      end
  



