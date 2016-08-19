      subroutine runCY_i(i,j,i1,f,Xtwiddle,Gtt,Gtwiddle,Shat2,Bzero1,N0)
      implicit none
C---  Expression for Eq. 5.57
C---  Calculates Ci, requires C00i
C---  Small terms of order Xtwiddle(0,j)*Cii
C---  Denominator Xtwiddle(i,j)

      include 'TRconstants.f' 
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,i,j,i1,n,m,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),
     . Gtwiddle(np,np),f(np),Gtt(np,np,np,np)
      double complex Shat2(np,z1max,-2:0),Bzero1(z1max,-2:0),
     . bit,pole


      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gtt(i,n,j,m)*f(n)
     . *(Shat2(m,i1,ep)-2d0*delta(m,i1)*Cv(cc00+N0,ep))
      enddo
      enddo

      pole=czip
      if (ep .gt. -2) pole=-4d0*Cv(czzi(i1)+N0,ep-1)
      Cv(ci(i1)+N0,ep)=  
     .  (Gtwiddle(i,j)*(6d0*Cv(czzi(i1)+N0,ep)+pole-Bzero1(i1,ep))
     . +bit+Xtwiddle(0,j)*Cv(cii(z2(i,i1))+N0,ep))/Xtwiddle(i,j)
      enddo
      return
      end
  



