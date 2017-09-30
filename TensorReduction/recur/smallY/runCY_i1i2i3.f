      subroutine runCY_i1i2i3(i,j,i1,i2,i3,
     . f,Xtwiddle,Gtt,Gtwiddle,Shat4,Bzero3,N0)
      implicit none
C---  Expression for Eq. 5.61
      include 'TRconstants.f' 
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,i,j,i1,i2,i3,n,m,np
      parameter(np=2)
      double precision 
     . Xtwiddle(0:np,0:np),Gtwiddle(np,np),f(np),Gtt(np,np,np,np)
      double complex Shat4(np,z3max,-2:0),Bzero3(z3max,-2:0),
     . bit,pole


      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gtt(i,n,j,m)*f(n)
     . *(Shat4(m,z3(i1,i2,i3),ep)
     . -2d0*delta(m,i1)*Cv(czzii(z2(i2,i3))+N0,ep)
     . -2d0*delta(m,i2)*Cv(czzii(z2(i3,i1))+N0,ep)
     . -2d0*delta(m,i3)*Cv(czzii(z2(i1,i2))+N0,ep))
      enddo
      enddo

      pole=czip
      if (ep .gt. -2) pole=-4*Cv(czziii(z3(i1,i2,i3))+N0,ep-1)
      Cv(ciii(z3(i1,i2,i3))+N0,ep)=
     .  (Gtwiddle(i,j)
     . *(10*Cv(czziii(z3(i1,i2,i3))+N0,ep)+pole-Bzero3(z3(i1,i2,i3),ep))
     . +bit+Xtwiddle(0,j)*Cv(ciiii(z4(i,i1,i2,i3))+N0,ep))/Xtwiddle(i,j)

      enddo
      return
      end
  



