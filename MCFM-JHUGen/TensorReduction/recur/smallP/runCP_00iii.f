      subroutine runCP_00iii(i1,i2,i3,m0sq,Gr,Bzero3,N0)
      implicit none
      include 'TRconstants.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      integer ep,N0,i1,i2,i3,m,n,np
      parameter(np=2)
      double precision m0sq,Gr(np,np)
      double complex Bzero3(z3max,-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Cv(ciiiii(z5(n,m,i1,i2,i3))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Cv(czziii(z3(i1,i2,i3))+N0,ep-1)
      
      Cv(czziii(z3(i1,i2,i3))+N0,ep)=
     . (pole
     . +2d0*Bzero3(z3(i1,i2,i3),ep)
     . +2d0*m0sq*Cv(ciii(z3(i1,i2,i3))+N0,ep)
     . -bit)/20d0
      enddo
      
      return
      end
