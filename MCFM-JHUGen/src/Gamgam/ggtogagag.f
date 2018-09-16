      function ggtogagag()
      implicit none
      include 'types.f'
      real(dp):: ggtogagag
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'zprods_com.f'
      real(dp):: fac,Qsumsq
      complex(dp):: Asum(2,2,2,2,2),
     & A1(2,2,2,2,2),A2(2,2,2,2,2),A3(2,2,2,2,2),
     & A4(2,2,2,2,2),A5(2,2,2,2,2),A6(2,2,2,2,2),
     & A7(2,2,2,2,2),A8(2,2,2,2,2),A9(2,2,2,2,2),
     & A10(2,2,2,2,2),A11(2,2,2,2,2),A12(2,2,2,2,2)
      integer:: h1,h2,h3,h4,h5
      Qsumsq=Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2
C-----factor includes statistical factor and is averaged over initial colors and spins
C-----cf hep-ph/0206194v2, Eq.(8)
      fac=esq**2*gsq**3/(fourpi)**4*xn/V*Qsumsq**2     
c      call spinoru(5,p,za,zb)
      call Aboxfill(1,2,3,4,5,za,zb,A1)
      call Aboxfill(1,2,4,3,5,za,zb,A2)
      call Aboxfill(1,3,2,4,5,za,zb,A3)
      call Aboxfill(3,1,2,4,5,za,zb,A4)

      call Aboxfill(1,4,2,3,5,za,zb,A5)
      call Aboxfill(3,1,4,2,5,za,zb,A6)
      call Aboxfill(1,3,4,2,5,za,zb,A7)
      call Aboxfill(1,4,3,2,5,za,zb,A8)

      call Aboxfill(4,1,2,3,5,za,zb,A9)
      call Aboxfill(4,1,3,2,5,za,zb,A10)
      call Aboxfill(3,4,1,2,5,za,zb,A11)
      call Aboxfill(4,3,1,2,5,za,zb,A12)


      ggtogagag=0._dp
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
      Asum(h1,h2,h3,h4,h5)=
     &  A1(h1,h2,h3,h4,h5)+A2(h1,h2,h3,h4,h5)+A3(h1,h2,h3,h4,h5)
     & +A4(h1,h2,h3,h4,h5)+A5(h1,h2,h3,h4,h5)+A6(h1,h2,h3,h4,h5)
     & +A7(h1,h2,h3,h4,h5)+A8(h1,h2,h3,h4,h5)+A9(h1,h2,h3,h4,h5)
     & +A10(h1,h2,h3,h4,h5)+A11(h1,h2,h3,h4,h5)+A12(h1,h2,h3,h4,h5)
      ggtogagag=ggtogagag+abs(Asum(h1,h2,h3,h4,h5))**2

      enddo
      enddo
      enddo
      enddo
      enddo

c--- apply overall factor
      ggtogagag=fac*ggtogagag
 
      return
      end
