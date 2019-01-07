      subroutine runCP_0000(Gr,S0000,N0)
      implicit none
      include 'TRconstants.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      integer ep,N0,m,n,np
      parameter(np=2)
      double precision Gr(np,np)
      double complex S0000(-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Cv(czzii(z2(n,m))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Cv(cc0000+N0,ep-1)

c--- note: we have simplified the recursion relation by using DD (5.9)
c--- so that S0000 = 2B00(0)+2*m0sq*C00      
      Cv(cc0000+N0,ep)=
     . (pole
     . +S0000(ep)
     . -bit)/12d0
      enddo
      
      return
      end
