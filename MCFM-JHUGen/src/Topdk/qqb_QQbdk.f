      subroutine qqb_QQbdk(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     January, 2012.                                                   *
*     calculate the element squared and subtraction terms              *
*     for the process                                                  *
*                                                                      *
*     q(-p1) +qbar(-p2)=bbar(p6)+e-(p7)+nubar(p8)+nu(p3)+e+(p4)+b(p5)  *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
*                                                                      *
*     NOTE: this routine is a replacement for qqb_QQbdk_old.f,         *
*           including the effect of the b-quark mass. In the massive   *
*           case it's approximately 1.5 times slower than that routine *
*           [could be sped up by only returning (1,...,2) amplitudes   *
*            in the routine toppaironshell.f, if necessary]            *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'msq_cs.f'
      include 'plabel.f'
      integer:: j,k,cs,hb,hc,h12,j1,j2,h1,h2,j1max,j2min
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,qqb
      complex(dp)::  prop,
     & mtop(2,2),manti(2,2),mprod(2,2,2),mtot(2,2,2),
     & mabtot(2,2,2,2),mbatot(2,2,2,2),mqed(2,2,2,2),
     & mab(2,2,2,2),mba(2,2,2,2)
       parameter(j1max=1,j2min=2)

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      do cs=0,2
      msq_cs(cs,j,k)=0._dp
      enddo
      enddo
      enddo

      call toppaironshell(p,0,mprod,mab,mba)
      call tdecay(p,3,4,5,mtop)
      call adecay(p,7,8,6,manti)

      do hb=1,2
      do hc=1,2
      do h12=1,2
      mtot(hb,hc,h12)=czip
      do j1=1,j1max
      do j2=j2min,2
      mtot(hb,hc,h12)=mtot(hb,hc,h12)+
     & mtop(hb,j1)*mprod(j1,h12,j2)*manti(j2,hc)
      enddo
      enddo
c      write(6,*) hb,hc,h12,abs(mtot(hb,hc,h12))
      enddo
      enddo
      enddo
      prop=cplx2(zip,mt*twidth)**2
      fac=V*gwsq**4*gsq**2/abs(prop)**2
c--- include factor for hadronic decays of W
      if (plabel(3) == 'pp') fac=2._dp*xn*fac
      if (plabel(7) == 'pp') fac=2._dp*xn*fac

      qqb=0._dp
      do hb=1,2
      do hc=1,2
      do h12=1,2
      qqb=qqb+fac*aveqq*abs(mtot(hb,hc,h12))**2
      enddo
      enddo
      enddo

      do hb=1,2
      do hc=1,2
      do h1=1,2
      do h2=1,2
      mabtot(hb,h1,h2,hc)=czip
      mbatot(hb,h1,h2,hc)=czip

      do j1=1,j1max
      do j2=j2min,2
      mabtot(hb,h1,h2,hc)=mabtot(hb,h1,h2,hc)+
     & mtop(hb,j1)*mab(j1,h1,h2,j2)*manti(j2,hc)
      mbatot(hb,h1,h2,hc)=mbatot(hb,h1,h2,hc)+
     & mtop(hb,j1)*mba(j1,h1,h2,j2)*manti(j2,hc)
      mqed(hb,h1,h2,hc)=mabtot(hb,h1,h2,hc)+mbatot(hb,h1,h2,hc)
      enddo
      enddo
      enddo
      enddo

      enddo
      enddo

      msq_cs(1,0,0)=0._dp
      msq_cs(2,0,0)=0._dp
      msq_cs(0,0,0)=0._dp
      do hb=1,2
      do hc=1,2
      do h1=1,2
      do h2=1,2
      msq_cs(1,0,0)=msq_cs(1,0,0)+fac*avegg
     & *xn*(abs(mbatot(hb,h1,h2,hc))**2)
      msq_cs(2,0,0)=msq_cs(2,0,0)+fac*avegg
     & *xn*(abs(mabtot(hb,h1,h2,hc))**2)
      msq_cs(0,0,0)=msq_cs(0,0,0)+fac*avegg
     & *xn*(-abs(mqed(hb,h1,h2,hc))**2/xnsq)
      enddo
      enddo
      enddo
      enddo

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if ((j < 0) .or. (j > 0)) then
          msq(j,-j)=qqb
C Division of quark into color structures is arbitrary
          msq_cs(1,j,-j)=qqb/3._dp
          msq_cs(2,j,-j)=qqb/3._dp
          msq_cs(0,j,-j)=qqb/3._dp
      elseif (j == 0) then
          msq(0,0)=msq_cs(1,0,0)+msq_cs(2,0,0)+msq_cs(0,0,0)
C msq_cs(1,0,0)+msq_cs(2,0,0)+msq_cs(0,0,0)
      endif
      enddo

      return
      end
