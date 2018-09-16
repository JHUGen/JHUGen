      subroutine dkqqb_w_twdk_g_old(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: Francesco Tramontano                                     *
*     February, 2005.                                                  *
*     Real correction to W+t, radiation in the decay                   *
*                                                                      *
*    Matrix element squared and averaged over initial colours and spins*
*--- W+t production (nwz=-1)                                           *
*     q(-p1) + qbar(-p2) --> W + t(p5678)                              *
*                            |   |                                     *
*                            |   --> nu(p5) + e^+(p6) + b(p7) + f(p8)  *
*                            |                                         *
*                            --> e^-(p3) + nubar(p4)                   *
*--- W+tbar production (nwz=+1)                                        *
*     q(-p1) + qbar(-p2) --> W + t(p5678)                              *
*                            |   |                                     *
*                            |   --> e^-(p5) + nubar(p6) + b(p7) +f(p8)*      *
*                            |                                         *
*                            --> nu(p3) + e^+(p4)                      *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zprods_com.f'
      include 'nwz.f'
      integer:: i,j,k,nu,jpart,i3,i4,i5,i6,iq
      real(dp):: p(mxpart,4),q(mxpart,4),msq(-nf:nf,-nf:nf),
     & fac,dot,tautg,msq_gq,msq_qg
      complex(dp):: ampp(2,2),ampd(2,2),tot(2,2)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=zero
      enddo
      enddo

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
        write(6,*) 'Error in qqb_w_twdk_gdk, nwz is not +1 or -1 : ',nwz
	stop
      endif

      fac=aveqg*four*gsq**2*gwsq**4

c--- calculate the amplitudes for a g+q initial state
      tautg=two*(dot(p,1,5)+dot(p,1,6)+dot(p,1,7)+dot(p,1,8))
      do nu=1,4
        do jpart=1,8
          q(jpart,nu)=p(jpart,nu)
        enddo
      enddo
      do nu=1,4
        q(9,nu)=q(5,nu)+q(6,nu)+q(7,nu)+q(8,nu)-mt**2/tautg*q(1,nu)
      enddo

      call spinoru(9,q,za,zb)
c--- Note: call to tree now passes top mass as a parameter
      call tree(mt,1,2,i3,i4,9,ampp)
      call gs_wc_dg(q,1,2,i3,i4,i5,i6,7,8,9,ampd)

c--- Construct factored form of amplitudes by adding the helicities of
c--- the heavy quark. Perform sum over the squares of gluon helicities
      msq_gq=zero
      do i=1,2
        do j=1,2
          tot(i,j)=ampp(1,i)*ampd(1,j)+ampp(2,i)*ampd(2,j)
          msq_gq=msq_gq+xn*cf**2*abs(tot(i,j))**2
        enddo
      enddo

c--- calculate the amplitudes for a q+g initial state
      tautg=two*(dot(p,2,5)+dot(p,2,6)+dot(p,2,7)+dot(p,2,8))
      do nu=1,4
        q(9,nu)=q(5,nu)+q(6,nu)+q(7,nu)+q(8,nu)-mt**2/tautg*q(2,nu)
      enddo

c--- Construct factored form of amplitudes by adding the helicities of
c--- the heavy quark. Perform sum over the squares of gluon helicities
      call spinoru(9,q,za,zb)
      call tree(mt,2,1,i3,i4,9,ampp)
      call gs_wc_dg(q,2,1,i3,i4,i5,i6,7,8,9,ampd)

      msq_qg=zero
      do i=1,2
        do j=1,2
          tot(i,j)=ampp(1,i)*ampd(1,j)+ampp(2,i)*ampd(2,j)
          msq_qg=msq_qg+xn*cf**2*abs(tot(i,j))**2
        enddo
      enddo

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      if     ((j == 5*iq) .and. (k == 0)) then
          msq(j,k)=fac*msq_qg
      elseif ((j == 0) .and. (k == 5*iq)) then
          msq(j,k)=fac*msq_gq
      endif
      enddo
      enddo

      return
      end
