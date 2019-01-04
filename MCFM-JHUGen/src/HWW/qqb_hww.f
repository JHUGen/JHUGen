      subroutine qqb_hww(p,msq)
      implicit none
c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H -->  W^+ (nu(p3)+e^+(p4))+W^- (e^-(p5)+nubar(p6))
c---
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'process.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),s,s12
      double precision decay,gg,Asq
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

c---set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      s12=s(1,2)

      decay=gwsq**3*wmass**2*s(3,5)*s(4,6)
      decay=decay/((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
      decay=decay/((s(5,6)-wmass**2)**2+(wmass*wwidth)**2)
      decay=decay/((s12-hmass**2)**2+(hmass*hwidth)**2)
      if ((case .eq. 'HWW2lq') .or. (case .eq. 'HWWdkW')) then
        decay=2d0*xn*decay
      endif
      Asq=(as/(3d0*pi))**2/vevsq
      gg=0.5d0*V*Asq*s12**2

c---calculate propagators
      msq(0,0)=avegg*gg*decay

      return
      end
