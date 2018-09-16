      subroutine dkWqqb_QQb_v(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     January, 2012.                                                   *
*     calculate the element squared for the virtual corrections        *
*     to the hadronic decay of the W for the process                   *
*                                                                      *
*     q(-p1) +qbar(-p2)=bbar(p6)+e-(p7)+nubar(p8)+nu(p3)+e+(p4)+b(p5)  *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'plabel.f'
      integer:: j,k,hb,hc,h12,j1,j2,h1,h2,j1max,j2min
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,qqb,gg
      complex(dp)::  prop,
     & mtop(2,2),manti(2,2),mprod(2,2,2),mtot(2,2,2),
     & mabtot(2,2,2,2),mbatot(2,2,2,2),mqed(2,2,2,2),
     & mab(2,2,2,2),mba(2,2,2,2),
     & mtopv(2,2),mantiv(2,2),mtotv(2,2,2),
     & mabtotv(2,2,2,2),mbatotv(2,2,2,2),mqedv(2,2,2,2)
       parameter(j1max=2,j2min=1)  ! NEED TO RE-CHECK THESE

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call toppaironshell(p,0,mprod,mab,mba)
      call tdecay(p,3,4,5,mtop)
      call adecay(p,7,8,6,manti)

      do j1=1,2
      do j2=1,2
      mtopv(j1,j2)=czip
      mantiv(j1,j2)=czip
      enddo
      enddo
      if (plabel(3) == 'pp') then
        call tdecayW_v(p,3,4,5,mtopv)
      endif
      if (plabel(7) == 'pp') then
        call adecayW_v(p,7,8,6,mantiv)
      endif

      do hb=1,2
      do hc=1,2
      do h12=1,2
      mtot(hb,hc,h12)=czip
      mtotv(hb,hc,h12)=czip
      do j1=1,j1max
      do j2=j2min,2
      mtot(hb,hc,h12)=mtot(hb,hc,h12)+
     & mtop(hb,j1)*mprod(j1,h12,j2)*manti(j2,hc)
      mtotv(hb,hc,h12)=mtotv(hb,hc,h12)
     &+mtopv(hb,j1)*mprod(j1,h12,j2)*manti(j2,hc)
     &+mtop(hb,j1)*mprod(j1,h12,j2)*mantiv(j2,hc)
      enddo
      enddo
      enddo
      enddo
      enddo
      prop=cplx2(zip,mt*twidth)**2
      fac=V*gwsq**4*gsq**2/abs(prop)**2*ason2pi*CF
c--- include factor for hadronic decays of W
      if (plabel(3) == 'pp') fac=2._dp*xn*fac
      if (plabel(7) == 'pp') fac=2._dp*xn*fac
      qqb=0._dp
      do hb=1,2
      do hc=1,2
      do h12=1,2
      qqb=qqb+fac*aveqq*real(conjg(mtot(hb,hc,h12))*mtotv(hb,hc,h12))
      enddo
      enddo
      enddo

      do hb=1,2
      do hc=1,2
      do h1=1,2
      do h2=1,2
      mabtot(hb,h1,h2,hc)=czip
      mbatot(hb,h1,h2,hc)=czip
      mabtotv(hb,h1,h2,hc)=czip
      mbatotv(hb,h1,h2,hc)=czip

      do j1=1,j1max
      do j2=j2min,2
      mabtot(hb,h1,h2,hc)=mabtot(hb,h1,h2,hc)+
     & mtop(hb,j1)*mab(j1,h1,h2,j2)*manti(j2,hc)
      mbatot(hb,h1,h2,hc)=mbatot(hb,h1,h2,hc)+
     & mtop(hb,j1)*mba(j1,h1,h2,j2)*manti(j2,hc)
      mqed(hb,h1,h2,hc)=mabtot(hb,h1,h2,hc)+mbatot(hb,h1,h2,hc)

      mabtotv(hb,h1,h2,hc)=mabtotv(hb,h1,h2,hc)
     & +mtopv(hb,j1)*mab(j1,h1,h2,j2)*manti(j2,hc)
     & +mtop(hb,j1)*mab(j1,h1,h2,j2)*mantiv(j2,hc)
      mbatotv(hb,h1,h2,hc)=mbatotv(hb,h1,h2,hc)
     & +mtopv(hb,j1)*mba(j1,h1,h2,j2)*manti(j2,hc)
     & +mtop(hb,j1)*mba(j1,h1,h2,j2)*mantiv(j2,hc)
      mqedv(hb,h1,h2,hc)=mabtotv(hb,h1,h2,hc)+mbatotv(hb,h1,h2,hc)
      enddo
      enddo
      enddo
      enddo

      enddo
      enddo

      gg=0._dp
      do hb=1,2
      do hc=1,2
      do h1=1,2
      do h2=1,2
      gg=gg+fac*avegg*xn*(
     & +real(conjg(mbatot(hb,h1,h2,hc))*mbatotv(hb,h1,h2,hc))
     & +real(conjg(mabtot(hb,h1,h2,hc))*mabtotv(hb,h1,h2,hc))
     & +real(-conjg(mqed(hb,h1,h2,hc))*mqedv(hb,h1,h2,hc)/xnsq))
      enddo
      enddo
      enddo
      enddo

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if ((j < 0) .or. (j > 0)) then
          msq(j,-j)=qqb
      elseif (j == 0) then
          msq(0,0)=gg
      endif
      enddo

      return
      end
