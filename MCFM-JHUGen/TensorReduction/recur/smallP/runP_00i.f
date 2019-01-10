      subroutine runP_00i(i1,m0sq,Gr,Czero1,N0)
      implicit none
      include 'TRconstants.f'
      include 'pvDnames.f'
      include 'pvDv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      integer ep,N0,i1,m,n,np
      parameter(np=3)
      double precision m0sq,Gr(np,np)
      double complex Czero1(z1max,-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Dv(diii(z3(n,m,i1))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Dv(dzzi(i1)+N0,ep-1)
      
      Dv(dzzi(i1)+N0,ep)=
     . (pole
     . +2d0*Czero1(i1,ep)
     . +2d0*m0sq*Dv(di(i1)+N0,ep)
     . -bit)/12d0
      enddo
      
      return
      end
