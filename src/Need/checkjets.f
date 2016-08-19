      subroutine checkjets(jetsfound,qfinal,isub,failed)
c--- performs checks on the jets that are found by the clustering algorithm,
c--- to ensure that the correct number of heavy quark jets is found,
c--- in the right invariant mass range
c--- given integer 'jetsfound' jets with momenta 'qfinal', 'isub'
c--- returns logical 'failed'
      implicit none
      include 'bbproc.f'
      include 'constants.f'
      include 'jetlabel.f'
      include 'limits.f'
      include 'process.f'
      include 'removebr.f'
      include 'npart.f'
      include 'nproc.f'
      double precision qfinal(mxpart,4),m56,m57,m67
      logical failed
      integer countb,jetsfound,nbq,nba,isub
      
      failed=.false.

c--- check that particle 5 is a b for H+b, W+c, Z+Q and W+b+jet processes
c--- also check for t-channel single top when the BR is not removed
      if ( (case .eq. 'H_1jet')
     . .or.(case .eq. 'W_cjet') .or. (case .eq. 'Wcjet0')
     . .or.(case .eq. 'gQ__ZQ') .or. (case .eq. 'W_bjet')
     . .or.(case .eq. 'Z_bjet')
     . .or.((case.eq. 'bq_tpq') .and. (removebr .eqv. .false.))
     . .or.((case.eq. 'ttdkay') .and. (removebr .eqv. .false.))
     . .or.((case.eq. 'W_twdk') .and. (removebr .eqv. .false.))
     . .or.((case.eq. 'Wtdkay') .and. (removebr .eqv. .false.)) ) then
        countb=0
        if ((jetsfound .ge. 1) .and. ((jetlabel(1) .eq. 'bq')
     .    .or. (jetlabel(1) .eq. 'ba'))) countb=1
        if ((jetsfound .ge. 2) .and. ((jetlabel(2) .eq. 'bq')
     .    .or. (jetlabel(2) .eq. 'ba'))) countb=countb+1
        if ((jetsfound .ge. 3) .and. ((jetlabel(3) .eq. 'bq')
     .    .or. (jetlabel(3) .eq. 'ba'))) countb=countb+1
        if ((jetsfound .eq. 1) .and. (countb .eq. 0)) failed=.true.
        if (   (     (case .eq. 'bq_tpq') .or. (case .eq. 'ttdkay')
     .          .or. (case .eq. 'W_bjet') .or. (case .eq. 'Z_bjet') )
     .   .and. (countb .lt. 1) ) failed=.true.
        if ((nproc .eq. 132) .and. (jetsfound .eq. 2)
     .      .and. (countb .ne. 1)) failed=.true.
        if ((nproc .eq. 133) .and. (jetsfound .eq. 2)
     .      .and. (countb .ne. 2)) failed=.true.
        if ((nproc .eq. 342) .and. (jetsfound .eq. 2)
     .      .and. (countb .ne. 1) .and. (isub .eq. 0)) failed=.true.
      endif
            
c--- check that 5 and 6 are b and b-bar (if appropriate)
      if (bbproc) then
        call getbs(qfinal,nbq,nba)
        if ((nbq .eq. 0) .or. (nba .eq. 0)) failed=.true.    
      endif

c--- perform m56 mass cut if there are 2 or more jets found
c--- and there are at least 6 particles in the final state
      if ((npart .gt. 3) .and. (jetsfound .ge. 2)) then
        m56=(qfinal(5,4)+qfinal(6,4))**2
     .     -(qfinal(5,1)+qfinal(6,1))**2
     .     -(qfinal(5,2)+qfinal(6,2))**2
     .     -(qfinal(5,3)+qfinal(6,3))**2
        if (jetsfound .ge. 3) then
        m57=(qfinal(5,4)+qfinal(7,4))**2
     .     -(qfinal(5,1)+qfinal(7,1))**2
     .     -(qfinal(5,2)+qfinal(7,2))**2
     .     -(qfinal(5,3)+qfinal(7,3))**2
        m67=(qfinal(6,4)+qfinal(7,4))**2
     .     -(qfinal(6,1)+qfinal(7,1))**2
     .     -(qfinal(6,2)+qfinal(7,2))**2
     .     -(qfinal(6,3)+qfinal(7,3))**2
        m56=max(m56,max(m57,m67))
        endif
        if ((m56 .lt. bbsqmin) .or. (m56 .gt. bbsqmax)) then
          failed=.true.
        endif
      endif
                  
      return
      end
