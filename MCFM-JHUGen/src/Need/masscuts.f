      subroutine masscuts(p,*)
      implicit none
      include 'omp_lib.h'
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'limits.f'
      include 'kprocess.f'
      include 'cutoff.f'
      include 'interference.f'
      include 'nqcdjets.f'
      include 'nproc.f'
      include 'verbose.f'
      include 'mpicommon.f'
      logical:: first,VBSprocess
      real(dp):: p(mxpart,4),s34,s56,s36,s45,s3456,s78
      data first/.true./
      save first,VBSprocess
!$omp threadprivate(first,VBSprocess)

      if (first) then
      first=.false.
c--- do not allow a cut on m56 for W/Z+gamma processes, tau pairs, or DM
      if ( (kcase==kWgamma) .or. (kcase==kZgamma)
     & .or.(kcase==ktautau) .or. (kcase==kdm_jet)
     & .or.(kcase==kdm_gam) ) then
        bbsqmin=0._dp
        bbsqmax=81d8
      endif

c--- do not allow a cut on m34 for direct photon process, or tau pairs
      if ((kcase==kdirgam) .or. (kcase==ktautau)) then
        wsqmin=0._dp
        wsqmax=81d8
      endif

c---- do not allow cuts on m34 if doing gamma gamma (will be done elsewhere) 
      if((kcase==kgamgam) .or. (kcase==kgg2gam) .or. (kcase==kgmgmjt)) then 
         return 
      endif
     
c---- check to see whether this is a VBS process
      if ( ((nproc >= 220) .and. (nproc <= 229))
     & .or. (nproc == 2201) .or. (nproc == 2221)
     & .or. (nproc == 2231) .or. (nproc == 2241)
     & .or. (nproc == 2251) .or. (nproc == 2281)
     & .or. (nproc == 2291) ) then
        VBSprocess=.true.
      else
        VBSprocess=.false.
      endif

!$omp master
      if (rank.eq.0) then
      write(6,*)
      write(6,*) '****************** Basic mass cuts *****************'
      write(6,*) '*                                                  *'
      write(6,99) sqrt(wsqmin),'m34',sqrt(wsqmax)
      if ((nqcdjets < 2) .or. (VBSprocess)) then
      write(6,99) sqrt(bbsqmin),'m56',sqrt(bbsqmax)
      else
      write(6,98) sqrt(bbsqmin),'m(jet1,jet2)',sqrt(bbsqmax)
      endif
      if (interference) then
      write(6,99) sqrt(wsqmin),'m36',sqrt(wsqmax)
      write(6,99) sqrt(wsqmin),'m45',sqrt(wsqmax)
      endif
      write(6,97) m3456min,'m3456',m3456max
      if (VBSprocess) then
       write(6,*)'*               m(jet1,jet2) > 100 GeV             *'
      endif
      write(6,*) '****************************************************'
      endif
!$omp end master
      endif

c--- only apply cuts on s34 if vectors 3 and 4 are defined
      if ((abs(p(3,4)) > 1.e-8_dp) .and. (abs(p(4,4)) > 1.e-8_dp)) then
        s34=+(p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     &      -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2
c--- do not accept s34<cutoff either
        if ((s34 < max(wsqmin,cutoff)).or.(s34 > wsqmax)) return 1
      endif
         
c--- only apply cuts on s56 if vectors 5 and 6 are defined
      if ((abs(p(5,4)) > 1.e-8_dp) .and. (abs(p(6,4)) > 1.e-8_dp)) then
        s56=+(p(5,4)+p(6,4))**2-(p(5,1)+p(6,1))**2
     &      -(p(5,2)+p(6,2))**2-(p(5,3)+p(6,3))**2
c--- do not accept s56<cutoff either
c        if ((s56 < max(bbsqmin,cutoff)).or.(s56>bbsqmax)) return 1
        if ((s56 < bbsqmin).or.(s56>bbsqmax)) return 1
      endif
     
      if (interference) then
      s45=+(p(4,4)+p(5,4))**2-(p(4,1)+p(5,1))**2
     &    -(p(4,2)+p(5,2))**2-(p(4,3)+p(5,3))**2
      s36=+(p(3,4)+p(6,4))**2-(p(3,1)+p(6,1))**2
     &    -(p(3,2)+p(6,2))**2-(p(3,3)+p(6,3))**2
        if (wsqmin .ne. bbsqmin) then
        write(6,*) 'masscuts: min. cuts must be equal for interference'
        stop
        endif
        if ((s45 < bbsqmin) .or. (s45 > bbsqmax)) return 1
        if ((s36 < bbsqmin) .or. (s36 > bbsqmax)) return 1
      endif

      s3456=+(p(3,4)+p(4,4)+p(5,4)+p(6,4))**2
     &      -(p(3,1)+p(4,1)+p(5,1)+p(6,1))**2
     &      -(p(3,2)+p(4,2)+p(5,2)+p(6,2))**2
     &      -(p(3,3)+p(4,3)+p(5,3)+p(6,3))**2
      if (s3456 < m3456min**2) return 1
      if (s3456 > m3456max**2) return 1

      if (VBSprocess) then
        s78=+(p(7,4)+p(8,4))**2-(p(7,1)+p(8,1))**2
     &      -(p(7,2)+p(8,2))**2-(p(7,3)+p(8,3))**2
        if (s78 < 1d4) return 1
      endif

   97 format(' *          ',f8.2,'  <  ',a5,' < ',f8.2,'           *')
   98 format(' *      ',f8.2,'  <   ',a12,'  < ',f8.2,'      *')
   99 format(' *          ',f8.2,'  <   ',a3,'  < ',f8.2,'           *')

      return
      end

