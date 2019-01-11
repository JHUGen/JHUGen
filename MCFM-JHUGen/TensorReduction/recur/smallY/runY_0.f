      subroutine runY_0(i,j,f,Xtwiddle,Gtwiddle,Gtt,Shat1,Czero0,N0)
      implicit none
C---  Expression for Eq. 5.55
C---  Calculates D0, requires D00
C---  Small terms of order Xtwiddle(0,j)*Di
C---  Denominator Xtwiddle(i,j)
      include 'TRconstants.f' 
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,i,j,n,m,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np),f(np),
     . Gtt(np,np,np,np)
      double complex Shat1(np,-2:0),Czero0(-2:0),
     . bit,pole

      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gtt(i,n,j,m)*f(n)*Shat1(m,ep)
      enddo
      enddo

      pole=czip
      if (ep .gt. -2) pole=-4*Dv(dd00+N0,ep-1)
      Dv(dd0+N0,ep)=
     .  (Gtwiddle(i,j)*(2d0*Dv(dd00+N0,ep)+pole-Czero0(ep))
     . +bit+Xtwiddle(0,j)*Dv(di(i)+N0,ep))/Xtwiddle(i,j)
      enddo

c      write(6,*) 'Dv(dd0)',Dv(dd0+N0,-2),Dv(dd0+N0,-1),Dv(dd0+N0,0)
c      write(6,*) 'bit',Dv(dd0+N0,-2),Dv(dd0+N0,-1),Dv(dd0+N0,0)

      return
      end
  



