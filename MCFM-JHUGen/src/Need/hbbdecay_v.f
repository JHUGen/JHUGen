      subroutine hbbdecay_v(p,ib,ibb,msq)
************************************************************************
*     Author: J.M. Campbell, June 2012                                 *
*                                                                      *
*     virtual matrix element contribution for the process of           *
*     Higgs decay  H --> b(ib)+b~(ibb)                                 *
*     with bottom mass included                                        *
************************************************************************
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'scale.f'
      include 'includect.f'
      integer ib,ibb
      double precision p(mxpart,4),s56,msq,
     & ddilog,besq,beta,BigL,ratio
      
      s56=2d0*(p(ib,4)*p(ibb,4)-p(ib,1)*p(ibb,1)
     &        -p(ib,2)*p(ibb,2)-p(ib,3)*p(ibb,3))+2d0*mb**2
      
      msq=xn*gwsq*mbsq/(4d0*wmass**2)
     & *2d0*(s56-4d0*mb**2)*Cf*ason2pi
      
c----Formula taken from Braaten and Leveille, PR D22, 715, 1980
      besq=1d0-4d0*mb**2/s56
      beta=sqrt(besq)
      ratio=(1d0-beta)/(1d0+beta)
      BigL=-log(ratio)
      msq=msq*(3d0*(epinv+log(musq/s56))
     & +((1d0+besq)/beta*BigL-2d0)*(epinv+log(musq/s56))
     & +2d0-log(mb**2/s56)+(2d0/beta-2*beta)*BigL
     & +(1d0+besq)/beta*(0.5d0*BigL**2-2d0*BigL*log(beta)
     & +2d0*ddilog(ratio)+2d0/3d0*pisq)
     & -3d0*(epinv+log(musq/s56)-log(mb**2/s56)+4d0/3d0))
C-----Last line is counterterm
C-----Gamma(1-e)*Gamma(1+e)=1+O(e^2)

      if (includect) then
        write(6,*) 'Counter-term not yet included in hbbdecay_v.f'
        stop
      endif

      return
      end
      
