      subroutine wbfcuts_jeppe(p,maxparts,passed)
      implicit none
      include 'types.f'
c--- given the event momenta in p (with maxparts entries),
c--- performs some generic WBF cuts
c---  a point that fails the cuts returns passed=.false.
c--- Note: this implements a modified version of Eq. (3.2) of CEZ
c---       paper, where any pair of three jets can satisy the criteria
c---       (as proposed by J. Andersen)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'plabel.f'
      logical:: passed
      integer:: j,k,maxparts,j1,j2
      real(dp):: etarap,p(mxpart,4),etaj1,etaj2

      passed=.true.

***************************** START CUTS *******************************

c--- check for pairs of jets
      do j=3,maxparts
        j1=0
        if (  (plabel(j) == 'pp') .or. (plabel(j) == 'qj')
     &   .or. (plabel(j) == 'bq') .or. (plabel(j) == 'ba')) j1=j
        do k=j+1,maxparts
        j2=0
          if (  (plabel(k) == 'pp') .or. (plabel(k) == 'qj')
     &   .  or. (plabel(k) == 'bq') .or. (plabel(k) == 'ba')) j2=k
c--- only look at jet-jet combinations
        if (( j1 == 0) .or. (j2 == 0)) goto 77

c--- now do cuts
          etaj1=etarap(j1,p)
          etaj2=etarap(j2,p)

c--- ensure a rapidity gap of at least 4.2 between the jets
c--- and ensure also that they jets lie in opposite hemispheres
      if ((abs(etaj1-etaj2) > 4.2_dp) .and. (etaj1*etaj2 <= 0._dp))
     &  goto 999


   77   continue

      enddo
      enddo


********************** END OF CUT CROSS SECTIONS ***********************

c--- if we get here, the cuts have failed
      passed=.false.

c--- this is the alternate return if the cuts are passed
  999 continue

      return

      end


