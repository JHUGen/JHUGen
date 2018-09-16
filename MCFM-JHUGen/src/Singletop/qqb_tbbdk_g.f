      subroutine qqb_tbbdk_g(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     January, 2012.                                                   *
*     calculate the LO element squared and subtraction terms           *
*     for the process                                                  *
*                                                                      *
*     u(-p1) +dbar(-p2)=t(nu(p3)+e+(p4)+b(p5))+bbar(p6)+g(p7)          *
*     or                                                               *
*     d(-p1) +ubar(-p2)=t~(e-(p3)+nu~(p4)+bb(p5))+b(p6)+g(p7)          *
*                                                                      *
*     Top (antitop) is kept strictly on-shell                          *
*     although all spin correlations are retained.                     *
*                                                                      *
*     NOTE: this routine is a replacement for qqb_tbb_g.f, including   *
*           the effect of the b-quark mass. In the massless case it is *
*           approximately 25% slower than that routine                 *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ckm.f'
      include 'nwz.f'
      integer:: j,k,hb,hc,hg,jmax,jmin
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,qa,aq,gq,qg,ag,ga
      complex(dp)::  prop
      complex(dp)::  mdecay(2,2),
     & miprodqa(2,2,2),miprodaq(2,2,2),
     & miprodga(2,2,2),miprodag(2,2,2),
     & miprodqg(2,2,2),miprodgq(2,2,2),
     & mfprodqa(2,2,2),mfprodaq(2,2,2),
     & mfprodga(2,2,2),mfprodag(2,2,2),
     & mfprodqg(2,2,2),mfprodgq(2,2,2),
     & mitotqa(2,2,2),mitotaq(2,2,2),
     & mitotga(2,2,2),mitotag(2,2,2),
     & mitotqg(2,2,2),mitotgq(2,2,2),
     & mftotqa(2,2,2),mftotaq(2,2,2),
     & mftotga(2,2,2),mftotag(2,2,2),
     & mftotqg(2,2,2),mftotgq(2,2,2)
      parameter(jmax=1,jmin=2)

      prop=cplx2(zip,mt*twidth)
      fac=xn**2*gwsq**4/abs(prop)**2*gsq*V/xn

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      if (nwz == +1) then
        call schantoponshellg(1,2,7,p,0,miprodqa,mfprodqa)
        call schantoponshellg(2,1,7,p,0,miprodaq,mfprodaq)
        call schantoponshellg(2,7,1,p,0,miprodgq,mfprodgq)
        call schantoponshellg(7,2,1,p,0,miprodga,mfprodga)
        call schantoponshellg(7,1,2,p,0,miprodag,mfprodag)
        call schantoponshellg(1,7,2,p,0,miprodqg,mfprodqg)
        call tdecay(p,3,4,5,mdecay)
      elseif (nwz == -1) then
        call schanatoponshellg(1,2,7,p,0,miprodqa,mfprodqa)
        call schanatoponshellg(2,1,7,p,0,miprodaq,mfprodaq)
        call schanatoponshellg(2,7,1,p,0,miprodgq,mfprodgq)
        call schanatoponshellg(7,2,1,p,0,miprodga,mfprodga)
        call schanatoponshellg(7,1,2,p,0,miprodag,mfprodag)
        call schanatoponshellg(1,7,2,p,0,miprodqg,mfprodqg)
        call adecay(p,3,4,5,mdecay)
      endif

c--- Calculate complete amplitudes and square
      qa=0._dp
      aq=0._dp
      ga=0._dp
      gq=0._dp
      ag=0._dp
      qg=0._dp
      do hb=1,2
      do hg=1,2
      do hc=1,2
      mitotqa(hb,hg,hc)=czip
      mitotaq(hb,hg,hc)=czip
      mitotga(hb,hg,hc)=czip
      mitotgq(hb,hg,hc)=czip
      mitotag(hb,hg,hc)=czip
      mitotqg(hb,hg,hc)=czip
      mftotqa(hb,hg,hc)=czip
      mftotaq(hb,hg,hc)=czip
      mftotga(hb,hg,hc)=czip
      mftotgq(hb,hg,hc)=czip
      mftotag(hb,hg,hc)=czip
      mftotqg(hb,hg,hc)=czip
      if (nwz == +1) then
      do j=1,jmax
      mitotqa(hb,hg,hc)=mitotqa(hb,hg,hc)+mdecay(hb,j)*miprodqa(j,hg,hc)
      mitotaq(hb,hg,hc)=mitotaq(hb,hg,hc)+mdecay(hb,j)*miprodaq(j,hg,hc)
      mitotgq(hb,hg,hc)=mitotgq(hb,hg,hc)+mdecay(hb,j)*miprodgq(j,hg,hc)
      mitotga(hb,hg,hc)=mitotga(hb,hg,hc)+mdecay(hb,j)*miprodga(j,hg,hc)
      mitotag(hb,hg,hc)=mitotag(hb,hg,hc)+mdecay(hb,j)*miprodag(j,hg,hc)
      mitotqg(hb,hg,hc)=mitotqg(hb,hg,hc)+mdecay(hb,j)*miprodqg(j,hg,hc)
      mftotqa(hb,hg,hc)=mftotqa(hb,hg,hc)+mdecay(hb,j)*mfprodqa(j,hg,hc)
      mftotaq(hb,hg,hc)=mftotaq(hb,hg,hc)+mdecay(hb,j)*mfprodaq(j,hg,hc)
      mftotgq(hb,hg,hc)=mftotgq(hb,hg,hc)+mdecay(hb,j)*mfprodgq(j,hg,hc)
      mftotga(hb,hg,hc)=mftotga(hb,hg,hc)+mdecay(hb,j)*mfprodga(j,hg,hc)
      mftotag(hb,hg,hc)=mftotag(hb,hg,hc)+mdecay(hb,j)*mfprodag(j,hg,hc)
      mftotqg(hb,hg,hc)=mftotqg(hb,hg,hc)+mdecay(hb,j)*mfprodqg(j,hg,hc)
      enddo
      elseif (nwz == -1) then
      do j=jmin,2
      mitotqa(hb,hg,hc)=mitotqa(hb,hg,hc)+miprodqa(hb,hg,j)*mdecay(j,hc)
      mitotaq(hb,hg,hc)=mitotaq(hb,hg,hc)+miprodaq(hb,hg,j)*mdecay(j,hc)
      mitotgq(hb,hg,hc)=mitotgq(hb,hg,hc)+miprodgq(hb,hg,j)*mdecay(j,hc)
      mitotga(hb,hg,hc)=mitotga(hb,hg,hc)+miprodga(hb,hg,j)*mdecay(j,hc)
      mitotag(hb,hg,hc)=mitotag(hb,hg,hc)+miprodag(hb,hg,j)*mdecay(j,hc)
      mitotqg(hb,hg,hc)=mitotqg(hb,hg,hc)+miprodqg(hb,hg,j)*mdecay(j,hc)
      mftotqa(hb,hg,hc)=mftotqa(hb,hg,hc)+mfprodqa(hb,hg,j)*mdecay(j,hc)
      mftotaq(hb,hg,hc)=mftotaq(hb,hg,hc)+mfprodaq(hb,hg,j)*mdecay(j,hc)
      mftotgq(hb,hg,hc)=mftotgq(hb,hg,hc)+mfprodgq(hb,hg,j)*mdecay(j,hc)
      mftotga(hb,hg,hc)=mftotga(hb,hg,hc)+mfprodga(hb,hg,j)*mdecay(j,hc)
      mftotag(hb,hg,hc)=mftotag(hb,hg,hc)+mfprodag(hb,hg,j)*mdecay(j,hc)
      mftotqg(hb,hg,hc)=mftotqg(hb,hg,hc)+mfprodqg(hb,hg,j)*mdecay(j,hc)
      enddo
      endif
c--- note: our usual definition of the s-channel process does not include
c---       diagrams where an incoming gluon is attached to the heavy quark line
c---       (these are included as t-channel instead; hence the 'zip' here)
      qa=qa+aveqq*fac*(abs(mitotqa(hb,hg,hc))**2
     &                +abs(mftotqa(hb,hg,hc))**2)
      aq=aq+aveqq*fac*(abs(mitotaq(hb,hg,hc))**2
     &                +abs(mftotaq(hb,hg,hc))**2)
      gq=gq+aveqg*fac*(abs(mitotgq(hb,hg,hc))**2
     &                +zip*abs(mftotgq(hb,hg,hc))**2)
      ga=ga+aveqg*fac*(abs(mitotga(hb,hg,hc))**2
     &                +zip*abs(mftotga(hb,hg,hc))**2)
      ag=ag+aveqg*fac*(abs(mitotag(hb,hg,hc))**2
     &                +zip*abs(mftotag(hb,hg,hc))**2)
      qg=qg+aveqg*fac*(abs(mitotqg(hb,hg,hc))**2
     &                +zip*abs(mftotqg(hb,hg,hc))**2)
      enddo
      enddo
      enddo

C---fill elements
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp

c--- Q-Qbar
      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=Vsq(j,k)*qa
c--- Qbar-Q
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=Vsq(j,k)*aq
c--- g-Q
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)
     &             +Vsq(-4,k)+Vsq(-5,k))*gq
c--- g-Qbar
      elseif ((j == 0) .and. (k < 0)) then
        msq(j,k)=(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)
     &           +Vsq(+4,k)+Vsq(+5,k))*ga
c--- Q-g
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)
     &             +Vsq(j,-4)+Vsq(j,-5))*qg
c--- Qbar-g
      elseif ((j < 0) .and. (k == 0)) then
        msq(j,k)=(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)
     &           +Vsq(j,+4)+Vsq(j,+5))*ag
      endif

      enddo
      enddo

      return
      end
