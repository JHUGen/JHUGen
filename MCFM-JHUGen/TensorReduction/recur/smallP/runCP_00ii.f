      subroutine runCP_00ii(i1,i2,m0sq,Gr,Bzero2,N0)
      implicit none
      include 'TRconstants.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      integer ep,N0,i1,i2,m,n,np
      parameter(np=2)
      double precision m0sq,Gr(np,np)
      double complex Bzero2(z2max,-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Cv(ciiii(z4(n,m,i1,i2))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Cv(czzii(z2(i1,i2))+N0,ep-1)
      
      Cv(czzii(z2(i1,i2))+N0,ep)=
     . (pole
     . +2d0*Bzero2(z2(i1,i2),ep)
     . +2d0*m0sq*Cv(cii(z2(i1,i2))+N0,ep)
     . -bit)/16d0
      enddo
      
      return
      end
