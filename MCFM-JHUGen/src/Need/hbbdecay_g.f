      subroutine hbbdecay_g(p,ib,ibb,ig,msq)
************************************************************************
*     Author: J.M. Campbell, June 2012                                 *
*                                                                      *
*     matrix element squared for the process of                        *
*     Higgs decay  H --> b(ib)+b~(ibb)+g(ig)                           *
*     with bottom mass included                                        *
************************************************************************
      implicit none
      include 'constants.f'
      include 'masses.f'
      integer ib,ibb,ig,j,k
      double precision p(mxpart,4),s,s56,s57,s67,msq,msqhbbg
      
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))
      s56=s(ib,ibb)+2d0*mb**2
      s57=s(ib,ig)
      s67=s(ibb,ig)
      
      msq=msqhbbg(s56,s57,s67)
      return
      end
      
      double precision function msqhbbg(s56,s57,s67)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'part.f'
      include 'msbarmasses.f'
      include 'ewcouple.f'
      include 'first.f'
      double precision mb_eff,massfrun,s56,s57,s67
      save mb_eff

      if (first) then
c--- run mb to appropriate scale
        if (part .eq. 'lord') then
          mb_eff=massfrun(mb_msbar,hmass,amz,1)
        else
          mb_eff=massfrun(mb_msbar,hmass,amz,2)
        endif
        first=.false.
      endif

      msqhbbg=xn*gwsq*mb_eff**2/(4d0*wmass**2)*gsq*Cf
     &*(2d0*(s56-4d0*mb**2)
     &*(-4d0*mb**2/s57**2-4d0*mb**2/s67**2+4d0*(s56-2d0*mb**2)/s57/s67)
     & +4d0*(2d0*s56+s57+s67-2d0*mb**2)/s57
     & +4d0*(2d0*s56+s57+s67-2d0*mb**2)/s67
     & -8d0*mb**2*(s67/s57**2+s57/s67**2))

      return
      end
