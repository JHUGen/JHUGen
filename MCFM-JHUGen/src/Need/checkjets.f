      subroutine checkjets(jetsfound,qfinal,isub,failed)
      implicit none
      include 'types.f'
c--- performs checks on the jets that are found by the clustering algorithm,
c--- to ensure that the correct number of heavy quark jets is found,
c--- in the right invariant mass range
c--- given integer:: 'jetsfound' jets with momenta 'qfinal', 'isub'
c--- returns logical:: 'failed'

      include 'bbproc.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'jetlabel.f'
      include 'limits.f'
      include 'kprocess.f'
      include 'removebr.f'
      include 'npart.f'
      include 'nproc.f'
      real(dp):: qfinal(mxpart,4),m56,m57,m67
      logical:: failed
      integer:: countb,jetsfound,nbq,nba,isub

      failed=.false.

c--- check that particle 5 is a b for H+b, W+c, Z+Q and W+b+jet processes
c--- also check for t-channel single top when the BR is not removed
      if ( (kcase==kH_1jet)
     & .or.(kcase==kW_cjet) .or. (kcase==kWcjet0)
     & .or.(kcase==kgQ__ZQ) .or. (kcase==kW_bjet)
     & .or.(kcase==kZ_bjet)
     & .or.((kcase==kbq_tpq) .and. (removebr .eqv. .false.))
     & .or.((kcase==kttdkay) .and. (removebr .eqv. .false.))
     & .or.((kcase==kW_twdk) .and. (removebr .eqv. .false.))
     & .or.((kcase==kWtdkay) .and. (removebr .eqv. .false.)) ) then
        countb=0
        if ((jetsfound >= 1) .and. ((jetlabel(1) == 'bq')
     &    .or. (jetlabel(1) == 'ba'))) countb=1
        if ((jetsfound >= 2) .and. ((jetlabel(2) == 'bq')
     &    .or. (jetlabel(2) == 'ba'))) countb=countb+1
        if ((jetsfound >= 3) .and. ((jetlabel(3) == 'bq')
     &    .or. (jetlabel(3) == 'ba'))) countb=countb+1
        if ((jetsfound == 1) .and. (countb == 0)) failed=.true.
        if (   (     (kcase==kbq_tpq) .or. (kcase==kttdkay)
     &          .or. (kcase==kW_bjet) .or. (kcase==kZ_bjet) )
     &   .and. (countb < 1) ) failed=.true.
        if ((nproc == 132) .and. (jetsfound == 2)
     &      .and. (countb .ne. 1)) failed=.true.
        if ((nproc == 133) .and. (jetsfound == 2)
     &      .and. (countb .ne. 2)) failed=.true.
        if ((nproc == 342) .and. (jetsfound == 2)
     &      .and. (countb .ne. 1) .and. (isub == 0)) failed=.true.
      endif

c--- check that 5 and 6 are b and b-bar (if appropriate)
      if (bbproc) then
        call getbs(qfinal,nbq,nba)
        if ((nbq == 0) .or. (nba == 0)) failed=.true.
      endif

c--- perform m56 mass cut if there are 2 or more jets found
c--- and there are at least 6 particles in the final state
      if ((npart > 3) .and. (jetsfound >= 2)) then
        m56=(qfinal(5,4)+qfinal(6,4))**2
     &     -(qfinal(5,1)+qfinal(6,1))**2
     &     -(qfinal(5,2)+qfinal(6,2))**2
     &     -(qfinal(5,3)+qfinal(6,3))**2
        if (jetsfound >= 3) then
        m57=(qfinal(5,4)+qfinal(7,4))**2
     &     -(qfinal(5,1)+qfinal(7,1))**2
     &     -(qfinal(5,2)+qfinal(7,2))**2
     &     -(qfinal(5,3)+qfinal(7,3))**2
        m67=(qfinal(6,4)+qfinal(7,4))**2
     &     -(qfinal(6,1)+qfinal(7,1))**2
     &     -(qfinal(6,2)+qfinal(7,2))**2
     &     -(qfinal(6,3)+qfinal(7,3))**2
        m56=max(m56,max(m57,m67))
        endif
c        if ((m56 < bbsqmin) .or. (m56 > bbsqmax)) then
c          failed=.true.
c        endif
      endif

      return
      end
