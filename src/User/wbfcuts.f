      subroutine wbfcuts(p,maxparts,passed)
c--- given the event momenta in p (with maxparts entries),
c--- performs some generic WBF cuts (a la Del Duca et al.)
c---  a point that fails the cuts returns passed=.false.
c--- Note: implements Eq. (3.2) of CEZ paper
      implicit none
      include 'constants.f'
      include 'plabel.f'
      logical passed
      integer j,maxparts,found,j1,j2
      double precision pt,etarap,p(mxpart,4),ptj,ptj1,ptj2,
     . etaj1,etaj2
      
      passed=.false.
      
***************************** START CUTS *******************************

      found=0
      ptj1=0d0
      ptj2=0d0

c--- check for 2 tagging jets (highest pt)
      do j=3,maxparts
        if (  (plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'qj')
     .   .or. (plabel(j) .eq. 'bq') .or. (plabel(j) .eq. 'ba')) then
          ptj=pt(j,p)
          if (ptj .gt. min(ptj1,ptj2)) then
            found=found+1
            if (found .eq. 1) then
              ptj1=ptj
              j1=j
            endif
            if (found .ge. 2) then
              if (ptj .gt. ptj1) then
                 ptj2=ptj1
                 ptj1=ptj
                 j2=j1
                 j1=j
              else
                 ptj2=ptj
                 j2=j
              endif
            endif
          endif
        endif
      enddo
      
      if (found .lt. 2) goto 999
                   
      etaj1=etarap(j1,p)
      etaj2=etarap(j2,p)
      
c--- ensure a rapidity gap of at least 4.2 between the tagging jets
      if (abs(etaj1-etaj2) .lt. 4.2d0) goto 999      
      
c--- ensure the tagging jets lie in opposite hemispheres
      if (etaj1*etaj2 .ge. 0d0) goto 999      

c      mj1j2=dsqrt(max(0d0,(p(j1,4)+p(j2,4))**2
c     . -(p(j1,1)+p(j2,1))**2-(p(j1,2)+p(j2,2))**2-(p(j1,3)+p(j2,3))**2))
c--- ensure the tagging jets have an invariant mass larger than 600 GeV
c--- NOT USED ANY LONGER
c      if (mj1j2 .lt. 600d0) goto 999      

********************** END OF CUT CROSS SECTIONS ***********************

      passed=.true.

c--- this is the alternate return if the cuts are failed    
  999 continue
      
      return
      
      end
      
     
