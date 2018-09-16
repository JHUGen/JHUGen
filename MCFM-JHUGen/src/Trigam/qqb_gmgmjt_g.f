      subroutine qqb_gmgmjt_g(p,msq)
      implicit none
      include 'types.f'
c--- matrix element squared for the process
c---    q(p1) + q~(p2) --> gam(p3) + gam(p4) + g(p5) + g(p6)
c--- (and all crossings)
c---
c--- J. M. Campbell, March 2013

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      integer:: j,k
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),
     & ampsq_2gam2g,qqb,qbq,qg,gq,qbg,gqb,gg,fac,cfac,qqb_4q,
     & ampsq_qq(2,2),ampsqid_qq(2,2),
     & ampsq_qqb(2,2),ampsqid_qqb(2,2),ampsq_qqb_anni(2,2),
     & xoppo,xsame
      integer,parameter::ii(6)=(/1,2,1,2,1,2/)


      call spinoru(6,p,za,zb)

c--- overall coupling, color and photon identical particle factors
      fac=8._dp*esq**2*gsq**2*cf*xn**2/2._dp

c      qbq=aveqq*fac*ampsq_2gam2g(5,6,3,4,2,1,za,zb)

c--- extra stat. factor for identical gluons in qqb case
      qqb=aveqq*fac/2._dp*ampsq_2gam2g(6,5,3,4,1,2,za,zb)
      qg=aveqg*fac*ampsq_2gam2g(2,6,3,4,1,5,za,zb)
      gq=aveqg*fac*ampsq_2gam2g(1,6,3,4,2,5,za,zb)
      gg=avegg*fac*ampsq_2gam2g(1,2,3,4,6,5,za,zb)
      qbq=qqb
      qbg=qg
      gqb=gq

      msq(:,:)=0._dp
      do j=1,nf
        cfac=Q(j)**4
        msq(j,-j)=cfac*qqb
        msq(-j,j)=cfac*qbq
        msq(j,0)=cfac*qg
        msq(0,j)=cfac*gq
        msq(-j,0)=cfac*qbg
        msq(0,-j)=cfac*gqb
      enddo
c--- assume 2 up-type flavors of quark and (nf-2) down-type
      msq(0,0)=gg*(2._dp*Q(2)**4+(nf-2)*Q(1)**4)


c--- overall coupling, color and photon identical particle factors
      fac=aveqq*8._dp*esq**2*gsq**2*cf*xn/2._dp
c--- 4-quark matrix elements
      call ampsq_2gam2q(1,5,2,6,3,4,za,zb,ampsq_qq,ampsqid_qq)
      call ampsq_2gam2q(1,5,6,2,3,4,za,zb,ampsq_qqb,ampsqid_qqb)
      call ampsq_2gam2q(1,2,6,5,3,4,za,zb,ampsq_qqb_anni,ampsqid_qqb)


c--- protection for har.e-_dpcoding used below
      if (nf .ne. 5) then
        write(6,*) 'Routine qqb_gmgmjt_g.f har.e-_dpcoded for nf=5!'
        stop
      endif

      do j=1,nf
      do k=1,nf

      if (j == k) then
c--- q-q identical
        msq(j,k)=fac/2._dp*ampsqid_qq(ii(j),ii(k))
c--- q-qb with annihilation allowed
        if (mod(j,2) == 1) then ! i.e. down-type quark
           xoppo=2._dp
           xsame=2._dp
        else                      ! i.e. up-type quark
           xoppo=3._dp
           xsame=1._dp
        endif
        qqb_4q=fac*(
     &    xoppo*ampsq_qqb_anni(ii(j),3-ii(j))  ! e.g. uub ->ddb
     &   +xsame*ampsq_qqb_anni(ii(j),ii(j))    ! e.g. uub ->ccb
     &   +ampsqid_qqb(ii(j),ii(j)))            ! e.g. uub ->uub
      else
c--- q-q non-identical
        msq(j,k)=fac*ampsq_qq(ii(j),ii(k))
c--- q-qb with no annihilation
        qqb_4q=fac*ampsq_qqb(ii(j),ii(k))
      endif
      msq(j,-k)=msq(j,-k)+qqb_4q

c--- qb-q and qb-qb (by symmetry)
      msq(-j,k)=msq(-j,k)+qqb_4q
      msq(-j,-k)=msq(j,k)

      enddo
      enddo

c      msq(2,1)=aveqq*fac*ampsq_qq(2,1)
c      msq(1,2)=aveqq*fac*ampsq_qq(1,2)
c      msq(2,4)=aveqq*fac*ampsq_qq(2,2)
c      msq(1,3)=aveqq*fac*ampsq_qq(1,1)

c      msq(1,1)=aveqq*fac*ampsqid_qq(1,1)
c      msq(2,2)=aveqq*fac*ampsqid_qq(2,2)

      return
      end


