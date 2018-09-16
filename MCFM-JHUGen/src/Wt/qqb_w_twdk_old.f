      subroutine qqb_w_twdk_old(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: Francesco Tramontano                                     *
*     February, 2005.                                                  *
*    Matrix element squared and averaged over initial colours and spins*
*--- W+t production (nwz=-1)                                           *
*     q(-p1) + qbar(-p2) --> W + t(p567)                               *
*                            |   |                                     *
*                            |   --> nu(p5) + e^+(p6) + b(p7)          *
*                            |                                         *
*                            --> e^-(p3) + nubar(p4)                   *
*--- W+tbar production (nwz=+1)                                        *
*     q(-p1) + qbar(-p2) --> W + t(p567)                               *
*                            |   |                                     *
*                            |   --> e^-(p5) + nubar(p6) + b(p7)       *
*                            |                                         *
*                            --> nu(p3) + e^+(p4)                      *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'nwz.f'
      integer:: j,k,i3,i4,i5,i6,iq
      real(dp):: p(mxpart,4),tvec(4),msq(-nf:nf,-nf:nf),
     & msq_gq,msq_qg,fac
      complex(dp):: ampl_gq(2),ampl_qg(2)
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      common/zabprods/zab,zba

c---initialize
      msq_gq=zero
      msq_qg=zero

      msq(:,:)=zero

c--- set up lepton variables depending on whether it's t or tbar
      if     (nwz == -1) then
        i3=3
	i4=4
	i5=5
	i6=6
	iq=1 ! top quark
      elseif (nwz == +1) then
        i3=4
	i4=3
	i5=6
	i6=5
	iq=-1 ! antitop quark
      else
        write(6,*) 'Error in qqb_w_twdk, nwz is not +1 or -1 :   ',nwz
	stop
      endif

c--- overall factor
      fac=aveqg*(V/two)*two*gsq*gwsq**4

      do j=1,4
         tvec(j)=p(5,j)+p(6,j)+p(7,j)
      enddo

c---fill matrices of spinor products

      call spinoru(7,p,za,zb)
      call spinork(7,p,zab,zba,tvec)

      call wamp(mt,twidth,1,2,i3,i4,i5,i6,7,ampl_gq)
      call wamp(mt,twidth,2,1,i3,i4,i5,i6,7,ampl_qg)

c--- sum over gluon helicities
      do j=1,2
      msq_gq=msq_gq+abs(ampl_gq(j))**2
      msq_qg=msq_qg+abs(ampl_qg(j))**2
      enddo

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      msq(j,k)=zero
      if     ((j == +5*iq) .and. (k == 0)) then
          msq(j,k)=fac*msq_qg
      elseif ((j == 0) .and. (k == +5*iq)) then
          msq(j,k)=fac*msq_gq
      endif
      enddo
      enddo

      return
      end


