      subroutine dkqqb_w_twdk_g(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: Keith Ellis                                              *
*     May, 2012.                                                       *
*     Real correction to W+t, radiation in the decay                   *
*                                                                      *
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
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'nwz.f'
      integer:: j,k,ht,hb,hg,h2,hc,ha
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: fac,msq_qg,msq_gq,msq_qbg,msq_gqb
      complex(dp):: prop
      complex(dp):: mtop(2,2,2),manti(2,2,2),
     & mqg(2,2,2),mgq(2,2,2),mqbg(2,2,2),mgqb(2,2,2),
     & mtotqg(2,2,2,2),mtotgq(2,2,2,2),mtotqbg(2,2,2,2),mtotgqb(2,2,2,2)

c---initialize
      msq(:,:)=zero
      mtotqg(:,:,:,:)=czip
      mtotgq(:,:,:,:)=czip
      mtotqbg(:,:,:,:)=czip
      mtotgqb(:,:,:,:)=czip


      if (nwz == -1) then
         call Wtoponshell(1,2,p,1,mqg)
         call Wtoponshell(2,1,p,1,mgq)
         call tdecayg(p,5,6,7,8,mtop)
         do hb=1,2
         do hg=1,2
         do h2=1,2
         do hc=1,2
         do ht=1,2
         mtotqg(hb,hg,h2,hc)=mtotqg(hb,hg,h2,hc)
     &    +mtop(hb,hg,ht)*mqg(ht,h2,hc)
         mtotgq(hb,hg,h2,hc)=mtotgq(hb,hg,h2,hc)
     &    +mtop(hb,hg,ht)*mgq(ht,h2,hc)
         enddo
         enddo
         enddo
         enddo
         enddo
      else
         call Watoponshell(1,2,p,-1,mqbg)
         call Watoponshell(2,1,p,-1,mgqb)
         call adecayg(p,5,6,7,8,manti)
         do hb=1,2
         do hg=1,2
         do h2=1,2
         do hc=1,2
         do ha=1,2
         mtotqbg(hb,hg,h2,hc)=mtotqbg(hb,hg,h2,hc)
     &    +mqbg(hb,h2,ha)*manti(ha,hg,hc)
         mtotgqb(hb,hg,h2,hc)=mtotgqb(hb,hg,h2,hc)
     &    +mgqb(hb,h2,ha)*manti(ha,hg,hc)
         enddo
         enddo
         enddo
         enddo
         enddo
      endif

      prop=cplx2(zip,mt*twidth)
      fac=aveqg*V*gwsq**4*gsq/abs(prop)**2*gsq*V/xn

      msq_qg=zero
      msq_gq=zero
      msq_qbg=zero
      msq_gqb=zero
      do hb=1,2
      do hg=1,2
      do h2=1,2
      do hc=1,2
      msq_qg=msq_qg+abs(mtotqg(hb,hg,h2,hc))**2
      msq_gq=msq_gq+abs(mtotgq(hb,hg,h2,hc))**2
      msq_qbg=msq_qbg+abs(mtotqbg(hb,hg,h2,hc))**2
      msq_gqb=msq_gqb+abs(mtotgqb(hb,hg,h2,hc))**2
      enddo
      enddo
      enddo
      enddo

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      if     ((j == +5) .and. (k == 0)) then
          msq(j,k)=fac*msq_qg
      elseif ((j == -5) .and. (k == 0)) then
          msq(j,k)=fac*msq_qbg
      elseif ((j == 0) .and. (k == +5)) then
          msq(j,k)=fac*msq_gq
      elseif ((j == 0) .and. (k == -5)) then
          msq(j,k)=fac*msq_gqb
      endif
      enddo
      enddo

      return
      end


