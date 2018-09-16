      subroutine hbbdecay(p,ib,ibb,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell, June 2012                                 *
*                                                                      *
*     matrix element squared for the process of                        *
*     Higgs decay  H --> b(ib)+b~(ibb)                                 *
*     with bottom mass included                                        *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      integer:: ib,ibb
      real(dp):: p(mxpart,4),s56,msq,msqhbb
      
      s56=2._dp*(p(ib,4)*p(ibb,4)-p(ib,1)*p(ibb,1)
     &        -p(ib,2)*p(ibb,2)-p(ib,3)*p(ibb,3))+2._dp*mb**2
      
      msq=msqhbb(s56)
      
      return
      end
      
      
      
      function msqhbb(s)
      implicit none
      include 'types.f'
      real(dp):: msqhbb
      real(dp):: s
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'couple.f'
      include 'kpart.f'
      include 'msbarmasses.f'
      include 'ewcouple.f'
      include 'first.f'
      include 'hbbparams.f'

      msqhbb=xn*gwsq*mb_eff**2/(4._dp*wmass**2)*2._dp*(s-4._dp*mb**2)

      return
      end
      
      
