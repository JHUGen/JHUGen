      subroutine qqb_wtbndk(p,msq)
************************************************************************
*    W+t+b production, with massive b                                  *
*    Matrix element squared and averaged over initial colours and spins*
*--- W+t production (nwz=-1)                                           *
*     q(-p1) + qbar(-p2) --> W + t(p5) + b(p6)                         *
*                            |                                         *
*                            |                                         *
*                            |                                         *
*                            --> e^-(p3) + nubar(p4)                   *
************************************************************************
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer ia,ig,ib,it,j,k
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),fac,dot,prop,
     . msq_gg
      double complex ampgg_ga(2,2,2,2),ampgg_ag(2,2,2,2)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      fac=gsq**2*gwsq**2

      prop=(two*dot(p,3,4)-wmass**2)**2+(wmass*wwidth)**2

c--- calculate amplitudes for g+g in initial state
      call BBamps_nores(p,1,2,3,4,5,6,ampgg_ag)
      call BBamps_nores(p,2,1,3,4,5,6,ampgg_ga)
c      call BBamps(p,1,2,3,4,5,6,ampgg_ag)
c      call BBamps(p,2,1,3,4,5,6,ampgg_ga)
      
      msq_gg=0d0
c--- sum over helicities of gluons and massive quarks
      do ib=1,2
      do it=1,2
      do ia=1,2
      do ig=1,2
        msq_gg=msq_gg+xn*cf**2*(
     .        + abs(ampgg_ag(ia,ig,ib,it))**2
     .        + abs(ampgg_ga(ig,ia,ib,it))**2
     .        - one/xn/cf*dble(ampgg_ag(ia,ig,ib,it)
     .                 *dconjg(ampgg_ga(ig,ia,ib,it))))
      enddo
      enddo
      enddo
      enddo

      msq_gg=msq_gg/prop

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      msq(j,k)=0d0
      if ((j .eq. 0) .and. (k .eq. 0)) then
          msq(j,k)=avegg*fac*msq_gg
      endif
      enddo
      enddo

      return
      end

