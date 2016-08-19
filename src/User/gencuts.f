      logical function gencuts(pjet,njets)
************************************************************************
*   Author: J.M. Campbell, 26th November 2013                          *
*                                                                      *
*   This routine calls specific cutting routines if requested,         *
*   otherwise just uses the cuts specified in the input file           *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'process.f'
      include 'runstring.f'
      include 'first.f'
      logical failed,gencuts_Zt,gencuts_WZjj,gencuts_input
      logical, save :: makeVBScuts,makeATLAS_sscuts,
     & makeCMS_hzz,makeCMS_hzz_vbf
      integer njets
      double precision pjet(mxpart,4)
!$omp threadprivate(makeVBScuts,makeATLAS_sscuts)
!$omp threadprivate(makeCMS_hzz,makeCMS_hzz_vbf)
      
      if (first) then
        first=.false.
        makeVBScuts=(index(runstring,'VBS') .gt. 0)
        makeATLAS_sscuts=(index(runstring,'ATLAS_ss') .gt. 0)
        makeCMS_hzz=(index(runstring,'CMS_hzz') .gt. 0)
        makeCMS_hzz_vbf=(index(runstring,'CMS_hzz_vbf') .gt. 0)
      endif
      
      gencuts=.false.
      
      if (makeVBScuts) then
        call VBS(pjet,failed)
        if (failed) gencuts=.true.
        return
      endif

      if (makeATLAS_sscuts) then
        call ATLAS_ss(pjet,failed)
        if (failed) gencuts=.true.
        return
      endif

      if (makeCMS_hzz_vbf) then
        call CMS_hzz_vbf(pjet,failed)
        if (failed) gencuts=.true.
        return
      endif
      
      if (makeCMS_hzz) then
        call CMS_hzz(pjet,failed)
        if (failed) gencuts=.true.
        return
      endif
      
c--- Default: use the cuts from the input file
      gencuts=gencuts_input(pjet,njets)
        
      return


      end
 



 
 
c      if ( (case .eq. 'HZZ_4l')
c     & .or.(case .eq. 'HZZ_tb')
c     & .or.(case .eq. 'HZZint')
c     & .or.(case .eq. 'HZZH+i')
c     & .or.(case .eq. 'ggZZ4l') 
c     & .or.(case .eq. 'HZZqgI')) then 
c        call CMS_hzz(pjet,failed)
c        if (failed) gencuts=.true.
c        return
c      endif

c      if ( (case .eq. 'HWW_4l')
c     & .or.(case .eq. 'HWW_tb')
c     & .or.(case .eq. 'HWWint')
c     & .or.(case .eq. 'HWWH+i')
c     & .or.(case .eq. 'ggWW4l') 
c     & .or.(case .eq. 'WWqqbr')) then 
c        call ATLAS_hww2013(pjet,failed)
c        if (failed) gencuts=.true.
c        return
c      endif

c -- CMS FCNC cuts for Z_tdkj
c      if ( case .eq. 'Z_tdkj' .and. runstring(1:3) .eq. 'CMS' ) then
c         gencuts=gencuts_Zt(pjet,njets)
c         return
c      endif

c      if ( case .eq. 'WpmZjj' .and. runstring(1:3) .eq. 'CMS' ) then
c         gencuts=gencuts_WZjj(pjet,njets)
c         return
c      endif

 
