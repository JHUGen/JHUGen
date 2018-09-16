      subroutine dkW1qqb_QQb_gs(p,msqc)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     November, 2011.                                                  *
*     calculate the subtraction term for radiation in the              *
*     top quark decay for the process                                  *
*                                                                      *
*     q(-p1)+qbar(-p2) = nu(p3)+e+(p4)+b(p5)                           *
*                        +bbar(p6)+e-(p7)+nubar(p8)+g(p9)              *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained. B-quark is taken to be either massive or massless. *
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
      real(dp):: msqc(maxd,-nf:nf,-nf:nf),p(mxpart,4)
      real(dp)::
     & msq39_4(-nf:nf,-nf:nf),msq49_3(-nf:nf,-nf:nf),
     & dummyv(-nf:nf,-nf:nf),sub39_4(4),sub49_3(4),dsubv
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

      call dips(1,p,3,9,4,sub39_4,dsubv,msq39_4,dummyv,
     & qqb_QQbdk,donothing_gvec)
      call dips(2,p,4,9,3,sub49_3,dsubv,msq49_3,dummyv,
     & qqb_QQbdk,donothing_gvec)
      do j=-nf,nf
        k=-j
        msqc(1,j,k)=sub39_4(qq)*msq39_4(j,k)*2._dp*cf
        msqc(2,j,k)=sub49_3(qq)*msq49_3(j,k)*2._dp*cf
      enddo

      return
      end

