      subroutine dkW1qqb_QQb_g(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     January, 2012.                                                   *
*     calculate the element squared and subtraction terms              *
*     for the process                                                  *
*                                                                      *
*     q(-p1) +qbar(-p2)=nu(p3)+e+(p4)+b(p5)+bbar(p6)+e-(p7)+nubar(p8)  *
*     +g(p9) radiated from top in decay                                *
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
      integer:: j,k,hb,hc,h12,j1,j2,h1,h2,hg,j2min
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,qqb,gg
      complex(dp)::  prop
      complex(dp)::  manti(2,2),mprod(2,2,2),
     & mab(2,2,2,2),mba(2,2,2,2)
      complex(dp):: Wtop(2,2,2),Wtot(2,2,2,2),
     & Wabtot(2,2,2,2,2),Wbatot(2,2,2,2,2),Wqed(2,2,2,2,2)
      parameter(j2min=2)

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call toppaironshell(p,1,mprod,mab,mba)
      call tdecayWg(p,3,4,5,9,Wtop)
      call adecay(p,7,8,6,manti)

      prop=cplx2(zip,mt*twidth)**2
      fac=V*gwsq**4*gsq**2/abs(prop)**2*gsq*V/xn
c--- include factor for hadronic decays of W
      if (plabel(3) == 'pp') fac=2._dp*xn*fac
      if (plabel(7) == 'pp') fac=2._dp*xn*fac

      do hb=1,2
      do hg=1,2
      do hc=1,2
      do h12=1,2
      Wtot(hb,hg,hc,h12)=czip
      do j1=1,2
      do j2=j2min,2
      Wtot(hb,hg,hc,h12)=Wtot(hb,hg,hc,h12)+
     & Wtop(hb,hg,j1)*mprod(j1,h12,j2)*manti(j2,hc)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      qqb=0._dp
      do hb=1,2
      do hg=1,2
      do hc=1,2
      do h12=1,2
      qqb=qqb+fac*aveqq*abs(Wtot(hb,hg,hc,h12))**2
      enddo
      enddo
      enddo
      enddo

c--- gg amplitudes
      do hb=1,2
      do hg=1,2
      do hc=1,2
      do h1=1,2
      do h2=1,2
      Wabtot(hb,hg,h1,h2,hc)=czip
      Wbatot(hb,hg,h1,h2,hc)=czip

      do j1=1,2
      do j2=j2min,2
      Wabtot(hb,hg,h1,h2,hc)=Wabtot(hb,hg,h1,h2,hc)+
     & Wtop(hb,hg,j1)*mab(j1,h1,h2,j2)*manti(j2,hc)
      Wbatot(hb,hg,h1,h2,hc)=Wbatot(hb,hg,h1,h2,hc)+
     & Wtop(hb,hg,j1)*mba(j1,h1,h2,j2)*manti(j2,hc)
      Wqed(hb,hg,h1,h2,hc)=Wabtot(hb,hg,h1,h2,hc)+Wbatot(hb,hg,h1,h2,hc)
      enddo
      enddo

      enddo
      enddo
      enddo
      enddo
      enddo

      gg=0._dp
      do hb=1,2
      do hg=1,2
      do hc=1,2
      do h1=1,2
      do h2=1,2
      gg=gg+fac*avegg*xn
     & *(abs(Wabtot(hb,hg,h1,h2,hc))**2+abs(Wbatot(hb,hg,h1,h2,hc))**2
     & -abs(Wqed(hb,hg,h1,h2,hc))**2/xnsq)
      enddo
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
