c--- note that p has mxpart entries, pjet has 20
      function deltarj(i,j,p,pjet)
      implicit none
      include 'types.f'
      real(dp):: deltarj
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),pjet(20,4),phi1,phi2,etarap,dphi,
     & etarap20
      integer:: i,j
      
      phi1=atan2(p(i,1),p(i,2))
      phi2=atan2(pjet(j,1),pjet(j,2))
      dphi=phi1-phi2
      if (dphi > pi) dphi=twopi-dphi
      if (dphi < -pi) dphi=twopi+dphi
      deltarj=(etarap(i,p)-etarap20(j,pjet))**2+dphi**2
      deltarj=sqrt(deltarj)
      return
      end

c--- note that p has 20 entries
      function deltarq(i,j,p)
      implicit none
      include 'types.f'
      real(dp):: deltarq
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(20,4),phi1,phi2,etarap20,dphi
      integer:: i,j

      phi1=atan2(p(i,1),p(i,2))
      phi2=atan2(p(j,1),p(j,2))
      dphi=phi1-phi2
      if (dphi > pi) dphi=twopi-dphi
      if (dphi < -pi) dphi=twopi+dphi
      deltarq=(etarap20(i,p)-etarap20(j,p))**2+dphi**2
      deltarq=sqrt(deltarq)

      return
      end

c--- note that p has 20 entries
      function etarap20(j,p)
      implicit none
      include 'types.f'
      real(dp):: etarap20
      
C---returns the value of the pseudorapidity
      integer:: j
      real(dp):: p(20,4)
      etarap20=sqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)
      etarap20=(etarap20+p(j,3))/(etarap20-p(j,3))
      if (etarap20 < 1.e-13_dp) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      etarap20=100._dp
      else
      etarap20=0.5_dp*log(etarap20)
      endif
      return
      end

