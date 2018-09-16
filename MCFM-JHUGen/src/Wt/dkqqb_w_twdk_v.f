      subroutine dkqqb_w_twdk_v(p,msq)
************************************************************************
*     Author: Keith Ellis                                              *
*     May, 2012.                                                       *
*    Virtual corrections in the decay, averaged over initial           *
*     colours and spins                                                *
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
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'nwz.f'
      integer j,k,ht,hb,h2,hc,ha
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf)
      double precision fac,msq_qg,msq_gq,msq_qbg,msq_gqb
      double complex prop
      double complex mtop(2,2),manti(2,2),mtopv(2,2),mantiv(2,2),
     & mqg(2,2,2),mgq(2,2,2),mqbg(2,2,2),mgqb(2,2,2),
     & mtotqg(2,2,2),mtotgq(2,2,2),mtotqbg(2,2,2),mtotgqb(2,2,2),
     & mtotqgv(2,2,2),mtotgqv(2,2,2),mtotqbgv(2,2,2),mtotgqbv(2,2,2)

c---initialize
      msq(:,:)=0d0
      mtotqg(:,:,:)=czip
      mtotgq(:,:,:)=czip
      mtotqbg(:,:,:)=czip
      mtotgqb(:,:,:)=czip
      mtotqgv(:,:,:)=czip
      mtotgqv(:,:,:)=czip
      mtotqbgv(:,:,:)=czip
      mtotgqbv(:,:,:)=czip

      if (nwz .eq. -1) then
         call Wtoponshell(1,2,p,0,mqg)
         call Wtoponshell(2,1,p,0,mgq)
         call tdecay(p,5,6,7,mtop)
         call tdecay_v(p,5,6,7,mtopv)
         do hb=1,2
         do h2=1,2
         do hc=1,2
         do ht=1,2
         mtotqg(hb,h2,hc)=mtotqg(hb,h2,hc)+mtop(hb,ht)*mqg(ht,h2,hc)
         mtotgq(hb,h2,hc)=mtotgq(hb,h2,hc)+mtop(hb,ht)*mgq(ht,h2,hc)
         mtotqgv(hb,h2,hc)=mtotqgv(hb,h2,hc)+mtopv(hb,ht)*mqg(ht,h2,hc)
         mtotgqv(hb,h2,hc)=mtotgqv(hb,h2,hc)+mtopv(hb,ht)*mgq(ht,h2,hc)
         enddo
         enddo
         enddo
         enddo
      else
         call Watoponshell(1,2,p,0,mqbg)
         call Watoponshell(2,1,p,0,mgqb)
         call adecay(p,5,6,7,manti)
         call adecay_v(p,5,6,7,mantiv)
         do hb=1,2
         do h2=1,2
         do hc=1,2
         do ha=1,2
         mtotqbg(hb,h2,hc)=mtotqbg(hb,h2,hc)+mqbg(hb,h2,ha)*manti(ha,hc)
         mtotgqb(hb,h2,hc)=mtotgqb(hb,h2,hc)+mgqb(hb,h2,ha)*manti(ha,hc)
         mtotqbgv(hb,h2,hc)=mtotqbgv(hb,h2,hc)
     &                     +mqbg(hb,h2,ha)*mantiv(ha,hc)
         mtotgqbv(hb,h2,hc)=mtotgqbv(hb,h2,hc)
     &                     +mgqb(hb,h2,ha)*mantiv(ha,hc)
         enddo
         enddo
         enddo
         enddo
      endif
      
      prop=dcmplx(zip,mt*twidth)
      fac=aveqg*V*gwsq**4*gsq/abs(prop)**2*ason2pi*CF

      msq_qg=0d0
      msq_gq=0d0
      msq_qbg=0d0
      msq_gqb=0d0
      do hb=1,2
      do h2=1,2
      do hc=1,2
      msq_qg=msq_qg+dble(dconjg(mtotqg(hb,h2,hc))*mtotqgv(hb,h2,hc))
      msq_gq=msq_gq+dble(dconjg(mtotgq(hb,h2,hc))*mtotgqv(hb,h2,hc))
      msq_qbg=msq_qbg+dble(dconjg(mtotqbg(hb,h2,hc))*mtotqbgv(hb,h2,hc))
      msq_gqb=msq_gqb+dble(dconjg(mtotgqb(hb,h2,hc))*mtotgqbv(hb,h2,hc))
      enddo
      enddo
      enddo

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      msq(j,k)=0d0
      if     ((j .eq. +5) .and. (k .eq. 0)) then
          msq(j,k)=fac*msq_qg
      elseif ((j .eq. -5) .and. (k .eq. 0)) then
          msq(j,k)=fac*msq_qbg
      elseif ((j .eq. 0) .and. (k .eq. +5)) then
          msq(j,k)=fac*msq_gq
      elseif ((j .eq. 0) .and. (k .eq. -5)) then
          msq(j,k)=fac*msq_gqb
      endif
      enddo
      enddo

      return
      end


