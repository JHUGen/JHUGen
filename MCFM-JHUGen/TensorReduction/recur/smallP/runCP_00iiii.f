      subroutine runCP_00iiii(i1,i2,i3,i4,m0sq,Gr,Bzero4,N0)
      implicit none
      include 'TRconstants.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      integer ep,N0,i1,i2,i3,i4,m,n,np
      parameter(np=2)
      double precision m0sq,Gr(np,np)
      double complex Bzero4(z4max,-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Cv(ciiiiii(z6(n,m,i1,i2,i3,i4))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Cv(czziiii(z4(i1,i2,i3,i4))+N0,ep-1)
      
      Cv(czziiii(z4(i1,i2,i3,i4))+N0,ep)=
     . (pole
     . +2d0*Bzero4(z4(i1,i2,i3,i4),ep)
     . +2d0*m0sq*Cv(ciiii(z4(i1,i2,i3,i4))+N0,ep)
     . -bit)/24d0
      enddo
      
      return
      end
