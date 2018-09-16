      subroutine runCP_00i(i1,m0sq,Gr,Bzero1,N0)
      implicit none
      include 'TRconstants.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      integer ep,N0,i1,m,n,np
      parameter(np=2)
      double precision m0sq,Gr(np,np)
      double complex Bzero1(z1max,-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Cv(ciii(z3(n,m,i1))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Cv(czzi(i1)+N0,ep-1)
      
      Cv(czzi(i1)+N0,ep)=
     . (pole
     . +2d0*Bzero1(i1,ep)
     . +2d0*m0sq*Cv(ci(i1)+N0,ep)
     . -bit)/12d0
      enddo
      
      return
      end
