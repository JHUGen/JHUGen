      subroutine hbbdecay_g(p,ib,ibb,ig,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell, June 2012                                 *
*                                                                      *
*     matrix element squared for the process of                        *
*     Higgs decay  H --> b(ib)+b~(ibb)+g(ig)                           *
*     with bottom mass included                                        *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      integer:: ib,ibb,ig,j,k
      real(dp):: p(mxpart,4),s,s56,s57,s67,msq,msqhbbg
      
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))
      s56=s(ib,ibb)+2._dp*mb**2
      s57=s(ib,ig)
      s67=s(ibb,ig)
      
      msq=msqhbbg(s56,s57,s67)
      return
      end
      
      function msqhbbg(s56,s57,s67)
      implicit none
      include 'types.f'
      real(dp):: msqhbbg
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'kpart.f'
      include 'msbarmasses.f'
      include 'ewcouple.f'
      include 'first.f'
      include 'hbbparams.f'
      real(dp):: s56,s57,s67

      msqhbbg=xn*gwsq*mb_eff**2/(4._dp*wmass**2)*gsq*Cf
     &*(2._dp*(s56-4._dp*mb**2)
     &*(-4._dp*mb**2/s57**2-4._dp*mb**2/s67**2
     &  +4._dp*(s56-2._dp*mb**2)/s57/s67)
     & +4._dp*(2._dp*s56+s57+s67-6._dp*mb**2)/s57
     & +4._dp*(2._dp*s56+s57+s67-6._dp*mb**2)/s67
     & -8._dp*mb**2*(s67/s57**2+s57/s67**2))

      return
      end
