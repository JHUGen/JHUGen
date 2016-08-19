      subroutine wbfcuts_jeppe(p,maxparts,passed)
c--- given the event momenta in p (with maxparts entries),
c--- performs some generic WBF cuts 
c---  a point that fails the cuts returns passed=.false.
c--- Note: this implements a modified version of Eq. (3.2) of CEZ
c---       paper, where any pair of three jets can satisy the criteria
c---       (as proposed by J. Andersen)
      implicit none
      include 'constants.f'
      include 'plabel.f'
      logical passed
      integer j,k,maxparts,j1,j2
      double precision etarap,p(mxpart,4),etaj1,etaj2
      
      passed=.true.
      
***************************** START CUTS *******************************

c--- check for pairs of jets
      do j=3,maxparts
        j1=0
        if (  (plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'qj')
     .   .or. (plabel(j) .eq. 'bq') .or. (plabel(j) .eq. 'ba')) j1=j
        do k=j+1,maxparts
        j2=0
          if (  (plabel(k) .eq. 'pp') .or. (plabel(k) .eq. 'qj')
     .   .  or. (plabel(k) .eq. 'bq') .or. (plabel(k) .eq. 'ba')) j2=k
c--- only look at jet-jet combinations
        if (( j1 .eq. 0) .or. (j2 .eq. 0)) goto 77

c--- now do cuts
          etaj1=etarap(j1,p)
          etaj2=etarap(j2,p)
      
c--- ensure a rapidity gap of at least 4.2 between the jets
c--- and ensure also that they jets lie in opposite hemispheres
      if ((abs(etaj1-etaj2) .gt. 4.2d0) .and. (etaj1*etaj2 .le. 0d0))
     .  goto 999      
      
      
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
      
     
