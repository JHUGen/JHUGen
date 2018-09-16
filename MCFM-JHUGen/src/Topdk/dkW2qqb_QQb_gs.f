      subroutine dkW2qqb_QQb_gs(p,msqc)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     November, 2011.                                                  *
*     calculate the subtraction term for radiation in the              *
*     anti-top quark decay for the process                             *
*                                                                      *
*     q(-p1)+qbar(-p2) = nu(p3)+e+(p4)+b(p5)                           *
*                        +bbar(p6)+e-(p7)+nubar(p8)+g(p9)              *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained. B-quark is taken to be either massless or massive  *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ptilde.f'
      include 'qcdcouple.f'
      include 'alfacut.f'
      include 'qqgg.f'
      include 'incldip.f'
      real(dp):: msqc(maxd,-nf:nf,-nf:nf),
     & p(mxpart,4)
      real(dp)::
     & msq79_8(-nf:nf,-nf:nf),msq89_7(-nf:nf,-nf:nf),
     & dummyv(-nf:nf,-nf:nf),sub79_8(4),sub89_7(4),dsubv
      integer:: j,k,nd
      external qqb_QQbdk,donothing_gvec

      ndmax=2

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msqc(nd,j,k)=0._dp
        incldip(nd)=.true.
      enddo
      enddo
      enddo

      call dips(1,p,7,9,8,sub79_8,dsubv,msq79_8,dummyv,
     & qqb_QQbdk,donothing_gvec)
      call dips(2,p,8,9,7,sub89_7,dsubv,msq89_7,dummyv,
     & qqb_QQbdk,donothing_gvec)
      do j=-nf,nf
        k=-j
        msqc(1,j,k)=sub79_8(qq)*msq79_8(j,k)*2._dp*cf
        msqc(2,j,k)=sub89_7(qq)*msq89_7(j,k)*2._dp*cf
      enddo

      return
      end

