      subroutine hbbdecay_v(p,ib,ibb,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell, June 2012                                 *
*                                                                      *
*     virtual matrix element contribution for the process of           *
*     Higgs decay  H --> b(ib)+b~(ibb)                                 *
*     with bottom mass included                                        *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'scale.f'
      include 'includect.f'
      integer:: ib,ibb
      real(dp):: p(mxpart,4),s56,msq,
     & ddilog,besq,beta,BigL,ratio,ctff_mqq,ff_mqq,xl56,mbar56
      
      s56=two*(p(ib,4)*p(ibb,4)-p(ib,1)*p(ibb,1)
     &        -p(ib,2)*p(ibb,2)-p(ib,3)*p(ibb,3))+two*mb**2
      
      if (mb < 1.e-6) then
        write(6,*) 'warning: mb=0 in input file'
        stop
      endif
      
c      msq=xn*gwsq*mbsq/(four*wmass**2)
c     & *two*(s56-four*mb**2)*Cf*ason2pi

      xl56=log(s56/musq)
      mbar56=mb/sqrt(s56)
c--- include two dipoles
      ctff_mqq=two*ff_mqq(1._dp,xl56,mbar56,1)
      
c--- Formula taken from Braaten and Leveille, PR D22, 715, 1980, Eq. (6)
C--- Last line is Hbb coupling counterterm of Eq. (13)
C--- No extra finite terms since Gamma(1-e)*Gamma(1+e)=1+O(e^2)
      besq=one-four*mb**2/s56
      beta=sqrt(besq)
      ratio=(one-beta)/(one+beta)
      BigL=-log(ratio)
      msq=
     & three*(epinv+log(musq/s56))
     & +((one+besq)/beta*BigL-two)*(epinv+log(musq/s56))
     & +two-log(mb**2/s56)+(two/beta-two*beta)*BigL
     & +(one+besq)/beta*(half*BigL**2-two*BigL*log(beta)
     & +two*ddilog(ratio)+two/three*pisq)
     & -3._dp*(epinv+log(musq/s56)-log(mb**2/s56)+4._dp/3._dp)
     & +ctff_mqq

c--- overall factor from Eq. (6) relative to LO width
      msq=CF*ason2pi*msq

c      if (includect) then
c        write(6,*) 'Counter-term not yet included in hbbdecay_v.f'
c        stop
c      endif

      return
      end
      
