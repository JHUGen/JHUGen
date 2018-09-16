      subroutine ttggHsingdriver(q,ampsq)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      real(dp):: q(mxpart,4),ampsq
      integer:: h1,h2,h3,h4
      complex(dp):: cab(-2:-1),cba(-2:-1),c00ab(-2:-1),c00ba(-2:-1),
     & ampAB(2,2,2,2),ampBA(2,2,2,2),
     & ampABs(2,2,2,2),ampBAs(2,2,2,2),ampABd(2,2,2,2),ampBAd(2,2,2,2),
     & tmp

      call ttggHamp(q,3,4,1,2,1,1,ampAB)
      call ttggHamp(q,3,4,2,1,1,1,ampBA)
      call singcoeffg(q,3,4,1,2,cab,cba,c00ab,c00ba)

      do h1=1,2
      do h2=1,2
      h3=1
      h4=2
      tmp=ampBA(h1,h2,h3,h4)
      ampBA(h1,h2,h3,h4)=ampBA(h1,h2,h4,h3)
      ampBA(h1,h2,h4,h3)=tmp
      enddo
      enddo


      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      ampABs(h1,h2,h3,h4)=
     & +cab(-2)*epinv**2*ampAB(h1,h2,h3,h4)
     & +cab(-1)*epinv*ampAB(h1,h2,h3,h4)
      ampBAs(h1,h2,h3,h4)=
     & +cba(-2)*epinv**2*ampBA(h1,h2,h3,h4)
     & +cba(-1)*epinv*ampBA(h1,h2,h3,h4)
      ampABd(h1,h2,h3,h4)=
     & +c00ab(-2)*epinv**2*ampAB(h1,h2,h3,h4)
     & +c00ab(-1)*epinv*ampAB(h1,h2,h3,h4)
      ampBAd(h1,h2,h3,h4)=
     & +c00ba(-2)*epinv**2*ampBA(h1,h2,h3,h4)
     & +c00ba(-1)*epinv*ampBA(h1,h2,h3,h4)
      enddo
      enddo
      enddo
      enddo

      ampsq=0d0
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
c  V*xn/4*(cab*mab*mabC+cba*mba*mbaC)
c  +V * ( 1/2*c00ba*mba*(mbaC+mabC)
c        +1/2*c00ab*mab*(mbaC+mabC))
c  -V/xn * (1/4*cba*mba*(mbaC+mabC)
c          +1/4*cab*mab*(mbaC+mabC));
C-----Factor of V*xn/4 removed
      ampsq=ampsq
     & +real(ampABs(h1,h2,h3,h4)*conjg(ampAB(h1,h2,h3,h4)))
     & +real(ampBAs(h1,h2,h3,h4)*conjg(ampBA(h1,h2,h3,h4)))
     & +2d0/xn*real((ampABd(h1,h2,h3,h4)+ampBAd(h1,h2,h3,h4))
     &  *conjg(ampAB(h1,h2,h3,h4)+ampBA(h1,h2,h3,h4)))
     & -1d0/xn**2*real((ampABs(h1,h2,h3,h4)+ampBAs(h1,h2,h3,h4))
     &  *conjg(ampAB(h1,h2,h3,h4)+ampBA(h1,h2,h3,h4)))
      enddo
      enddo
      enddo
      enddo
      return
      end
