      subroutine runCP_00(m0sq,Gr,Bzero0,N0)
      implicit none
      include 'TRconstants.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
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
      
      Cv(cc00+N0,ep)=
     . (pole
     . +2d0*Bzero0(ep)
     . +2d0*m0sq*Cv(cc0+N0,ep)
     . -bit)/8d0
      enddo
      
      return
      end
