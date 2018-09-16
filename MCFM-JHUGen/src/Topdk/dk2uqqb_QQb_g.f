      subroutine dk2uqqb_QQb_g(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     January, 2012.                                                   *
*     calculate the element squared and subtraction terms              *
*     for the process                                                  *
*                                                                      *
*     q(-p1) +qbar(-p2)=nu(p3)+e+(p4)+b(p5)+bbar(p6)+e-(p7)+nubar(p8)  *
*     +g(p9) radiated from antitop quark in decay                      *
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
      integer:: j,k,hb,hc,h12,j1,j2,h1,h2,hg
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,qqb,gg,msqtop,msqatop
      complex(dp)::  prop,
     & mtop(2,2),manti(2,2,2),mprod(2,2,2),
     & mab(2,2,2,2),mba(2,2,2,2),mqed(2,2,2,2)

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call toppaironshell(p,-1,mprod,mab,mba)
      call tdecay(p,3,4,5,mtop)
      call adecayg(p,7,8,6,9,manti)

C----calculate spin-averaged decay matrix elements
      msqtop=0._dp
      do hb=1,2
      do j1=1,2
      msqtop=msqtop+abs(mtop(hb,j1))**2
      enddo
      enddo

      msqatop=0._dp
      do hc=1,2
      do hg=1,2
      do j2=1,2
      msqatop=msqatop+abs(manti(j2,hg,hc))**2
      enddo
      enddo
      enddo

      prop=cplx2(zip,mt*twidth)**2
      fac=V*gwsq**4*gsq**2/abs(prop)**2*gsq*V/xn*(msqtop*msqatop/4._dp)
c--- include factor for W hadronic decays
      if (plabel(3) == 'pp') fac=2._dp*xn*fac
      if (plabel(7) == 'pp') fac=2._dp*xn*fac

c--- q-qbar amplitudes
      qqb=0._dp
      do j1=1,2
      do j2=1,2
      do h12=1,2
      qqb=qqb+fac*aveqq*abs(mprod(j1,h12,j2))**2
      enddo
      enddo
      enddo

c--- gg amplitudes
      mqed(:,:,:,:)=mab(:,:,:,:)+mba(:,:,:,:)
      gg=0._dp
      do j1=1,2
      do j2=1,2
      do h1=1,2
      do h2=1,2
      gg=gg+fac*avegg*xn
     &  *(abs(mab(j1,h1,h2,j2))**2+abs(mba(j1,h1,h2,j2))**2
     & -abs(mqed(j1,h1,h2,j2))**2/xnsq)
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
