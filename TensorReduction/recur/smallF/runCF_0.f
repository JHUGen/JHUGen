      subroutine runCF_0(m0sq,Gr,Bzero0,N0)
C---  Expression for rearrangement of Eq. 5.65
C---  Calculates C0, requires C00
C---  Small terms of order Gr(i,j)*Cij
C---  Denominator m0sq
      implicit none
      include 'TRconstants.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'pvweenumber.f' 
      integer ep,N0,m,n,np
      parameter(np=2)
      double precision m0sq,Gr(np,np)
      double complex Bzero0(-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Cv(cii(z2(n,m))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Cv(cc00+N0,ep-1)
      
      Cv(cc0+N0,ep)=
     . (8d0*Cv(cc00+N0,ep)
     . -pole
     . -2d0*Bzero0(ep)
     . +bit)/(2d0*m0sq)

      enddo
      
      return
      end
