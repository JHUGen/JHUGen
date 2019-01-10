      subroutine hbbdecay(p,ib,ibb,msq)
************************************************************************
*     Author: J.M. Campbell, June 2012                                 *
*                                                                      *
*     matrix element squared for the process of                        *
*     Higgs decay  H --> b(ib)+b~(ibb)                                 *
*     with bottom mass included                                        *
************************************************************************
      implicit none
      include 'constants.f'
      include 'masses.f'
      integer ib,ibb
      double precision p(mxpart,4),s56,msq,msqhbb
      
      s56=2d0*(p(ib,4)*p(ibb,4)-p(ib,1)*p(ibb,1)
     &        -p(ib,2)*p(ibb,2)-p(ib,3)*p(ibb,3))+2d0*mb**2
      
      msq=msqhbb(s56)
      
      return
      end
      
      
      
      double precision function msqhbb(s)
      implicit none
      double precision s
      include 'constants.f'
      include 'masses.f'
      include 'couple.f'
      include 'part.f'
      include 'msbarmasses.f'
      include 'ewcouple.f'
      include 'first.f'
      double precision mb_eff,massfrun
      save mb_eff
!$omp threadprivate(mb_eff)

      if (first) then
c--- run mb to appropriate scale
        if (part .eq. 'lord') then
          mb_eff=massfrun(mb_msbar,hmass,amz,1)
        else
          mb_eff=massfrun(mb_msbar,hmass,amz,2)
        endif  
        first=.false.
      endif

      msqhbb=xn*gwsq*mb_eff**2/(4d0*wmass**2)*2d0*(s-4d0*mb**2)

      return
      end
      
      
