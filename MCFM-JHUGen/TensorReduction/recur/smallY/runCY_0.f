      subroutine runCY_0(i,j,f,Xtwiddle,Gtwiddle,Gtt,Shat1,Bzero0,N0)
      implicit none
C---  Expression for Eq. 5.55
C---  Calculates C0, requires C00
C---  Small terms of order Xtwiddle(0,j)*Ci
C---  Denominator Xtwiddle(i,j)
      include 'TRconstants.f' 
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,i,j,n,m,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np),f(np),
     . Gtt(np,np,np,np)
      double complex Shat1(np,-2:0),Bzero0(-2:0),
     . bit,pole

      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gtt(i,n,j,m)*f(n)*Shat1(m,ep)
      enddo
      enddo

      pole=czip
      if (ep .gt. -2) pole=-4d0*Cv(cc00+N0,ep-1)
      Cv(cc0+N0,ep)=
     .  (Gtwiddle(i,j)*(4d0*Cv(cc00+N0,ep)+pole-Bzero0(ep))
     . +bit+Xtwiddle(0,j)*Cv(ci(i)+N0,ep))/Xtwiddle(i,j)
c      write(6,*) Gtwiddle(i,j)*(4d0*Cv(cc00+N0,ep))/Xtwiddle(i,j)
c      write(6,*) Gtwiddle(i,j)*(+pole)/Xtwiddle(i,j)
c      write(6,*) Gtwiddle(i,j)*(-Bzero0(ep))/Xtwiddle(i,j)
c      write(6,*) (bit)/Xtwiddle(i,j)
c      write(6,*) (Xtwiddle(0,j)*Cv(ci(i)+N0,ep))/Xtwiddle(i,j)
c      write(6,*)
      enddo

c      write(6,*) 'Cv(cc0)',Cv(cc0+N0,-2),Cv(cc0+N0,-1),Cv(cc0+N0,0)
c      write(6,*) 'bit',Cv(cc0+N0,-2),Cv(cc0+N0,-1),Cv(cc0+N0,0)

      return
      end
  



