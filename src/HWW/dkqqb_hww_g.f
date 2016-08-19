      subroutine dkqqb_hww_g(p,msq)
      implicit none
C----Author R.K. Ellis, May 2012
c----NLO matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c----with radiation in decay of W^- quark
c---- g(-p1)+g(-p2)-->H -->  W^+ (nu(p3)+e^+(p4))
c----                      + W^- (q(p5)+qbar(p6)+g(p7)) 
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),
     & s12,s34,s35,s46,s567,s56,s57,s67,s,Asq,decay,gg

c--- statement function
      s(j,k)=2d0*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     & -p(j,2)*p(k,2)-p(j,3)*p(k,3))

c---set msq=0 to initialize
      msq(:,:)=0d0

      s12=s(1,2)
      s34=s(3,4)
      s35=s(3,5)
      s46=s(4,6)
      s56=s(5,6)
      s57=s(5,7)
      s67=s(6,7)
      s567=s56+s57+s67
      decay=gwsq**3*wmass**2*s35*s46
c---calculate propagators
      decay=decay/((s34-wmass**2)**2+(wmass*wwidth)**2)
      decay=decay/((s567-wmass**2)**2+(wmass*wwidth)**2)
      decay=decay/((s12-hmass**2)**2+(hmass*hwidth)**2)
C-----Add factors for decay of W for all cases in which this routine is called
      decay=2d0*xn*decay
      decay=decay*2d0*gsq*CF/(s57*s67)
     & *(((s35+s(3,7))*(s56+s67)-s57*s(3,6))/s35     
     &  +((s46+s(4,7))*(s56+s57)-s67*s(4,5))/s46) 

      Asq=(as/(3d0*pi))**2/vevsq
      gg=0.5d0*V*Asq*s12**2

      msq(0,0)=avegg*gg*decay

      return
      end
