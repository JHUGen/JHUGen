      subroutine masscuts(p,*)
      implicit none
      include 'constants.f'
      include 'limits.f'
      include 'process.f'
      include 'cutoff.f'
      include 'interference.f'
      include 'nqcdjets.f'
      include 'nproc.f'
      logical first,VBSprocess
      double precision p(mxpart,4),s34,s56,s36,s45,s3456,s78
      data first/.true./
      save first,VBSprocess
!$omp threadprivate(first,VBSprocess)

      if (first) then
      first=.false.
c--- do not allow a cut on m56 for W/Z+gamma processes, tau pairs, or DM
      if ( (case .eq. 'Wgamma') .or. (case .eq. 'Zgamma')
     & .or.(case .eq. 'tautau') .or. (case .eq. 'dm_jet')
     & .or.(case .eq. 'dm_gam') ) then
        bbsqmin=0d0
        bbsqmax=81d8
      endif

c--- do not allow a cut on m34 for direct photon process, or tau pairs
      if ( (case .eq. 'dirgam')  .or. (case .eq. 'tautau')) then
        wsqmin=0d0
        wsqmax=81d8
      endif

c---- do not allow cuts on m34 if doing gamma gamma (will be done elsewhere)
      if(case .eq. 'gamgam') then
         return
      endif

c---- check to see whether this is a VBS process
      if ( ((nproc .ge. 220) .and. (nproc .le. 229))
     & .or. (nproc .eq. 2201) .or. (nproc .eq. 2221)
     & .or. (nproc .eq. 2231) .or. (nproc .eq. 2241)
     & .or. (nproc .eq. 2251) .or. (nproc .eq. 2281)
     & .or. (nproc .eq. 2291) ) then
        VBSprocess=.true.
      else
        VBSprocess=.false.
      endif

      write(6,*)
      write(6,*) '****************** Basic mass cuts *****************'
      write(6,*) '*                                                  *'
      write(6,99) dsqrt(wsqmin),'m34',dsqrt(wsqmax)
      if ((nqcdjets .lt. 2) .or. (VBSprocess)) then
      write(6,99) dsqrt(bbsqmin),'m56',dsqrt(bbsqmax)
      else
      write(6,98) dsqrt(bbsqmin),'m(jet1,jet2)',dsqrt(bbsqmax)
      endif
      if (interference) then
      write(6,99) dsqrt(wsqmin),'m36',dsqrt(wsqmax)
      write(6,99) dsqrt(wsqmin),'m45',dsqrt(wsqmax)
      endif
      write(6,97) m3456min,'m3456',m3456max
      if (VBSprocess) then
       write(6,*)'*               m(jet1,jet2) > 100 GeV             *'
      endif
      write(6,*) '****************************************************'
      endif

c--- only apply cuts on s34 if vectors 3 and 4 are defined
      if ((abs(p(3,4)) .gt. 1d-8) .and. (abs(p(4,4)) .gt. 1d-8)) then
        s34=+(p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     .      -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2
c--- do not accept s34<cutoff either
        if ((s34 .lt. max(wsqmin,cutoff)).or.(s34 .gt. wsqmax)) return 1
      endif

c--- only apply cuts on s56 if vectors 5 and 6 are defined
      if ((abs(p(5,4)) .gt. 1d-8) .and. (abs(p(6,4)) .gt. 1d-8)) then
        s56=+(p(5,4)+p(6,4))**2-(p(5,1)+p(6,1))**2
     .      -(p(5,2)+p(6,2))**2-(p(5,3)+p(6,3))**2
c--- do not accept s56<cutoff either
        if ((s56 .lt. max(bbsqmin,cutoff)).or.(s56.gt.bbsqmax)) return 1
      endif

      if (interference) then
      s45=+(p(4,4)+p(5,4))**2-(p(4,1)+p(5,1))**2
     .    -(p(4,2)+p(5,2))**2-(p(4,3)+p(5,3))**2
      s36=+(p(3,4)+p(6,4))**2-(p(3,1)+p(6,1))**2
     .    -(p(3,2)+p(6,2))**2-(p(3,3)+p(6,3))**2
        if (wsqmin .ne. bbsqmin) then
        write(6,*) 'masscuts: min. cuts must be equal for interference'
        stop
        endif
        if ((s45 .lt. bbsqmin) .or. (s45 .gt. bbsqmax)) return 1
        if ((s36 .lt. bbsqmin) .or. (s36 .gt. bbsqmax)) return 1
      endif

      s3456=+(p(3,4)+p(4,4)+p(5,4)+p(6,4))**2
     &      -(p(3,1)+p(4,1)+p(5,1)+p(6,1))**2
     &      -(p(3,2)+p(4,2)+p(5,2)+p(6,2))**2
     &      -(p(3,3)+p(4,3)+p(5,3)+p(6,3))**2
      if (s3456 .lt. m3456min**2) return 1
      if (s3456 .gt. m3456max**2) return 1

      if (VBSprocess) then
        s78=+(p(7,4)+p(8,4))**2-(p(7,1)+p(8,1))**2
     &      -(p(7,2)+p(8,2))**2-(p(7,3)+p(8,3))**2
        if (s78 .lt. 1d2) return 1
      endif

   97 format(' *         ',f8.2,'  <   ',a5,'  < ',f8.2,'          *')
   98 format(' *      ',f8.2,'  <   ',a12,'  < ',f8.2,'      *')
   99 format(' *          ',f8.2,'  <   ',a3,'  < ',f8.2,'           *')

      return
      end

