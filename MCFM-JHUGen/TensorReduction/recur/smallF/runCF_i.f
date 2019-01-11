      subroutine runCF_i(i1,m0sq,Gr,Bzero1,N0)
C---  Expression for rearrangement of Eq. 5.67
C---  Calculates Ci, requires C00i
C---  Small terms of order Gr(i,j)*Cijk
C---  Denominator m0sq
      implicit none
      include 'TRconstants.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'pvweenumber.f' 
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
      
      Cv(ci(i1)+N0,ep)=    
     . (12d0*Cv(czzi(i1)+N0,ep)
     . -pole
     . -2d0*Bzero1(i1,ep)
     . +bit)/(2d0*m0sq)

      enddo
      
      return
      end
