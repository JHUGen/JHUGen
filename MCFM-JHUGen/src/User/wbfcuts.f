      subroutine wbfcuts(p,maxparts,passed)
      implicit none
      include 'types.f'
c--- given the event momenta in p (with maxparts entries),
c--- performs some generic WBF cuts (a la Del Duca et al.)
c---  a point that fails the cuts returns passed=.false.
c--- Note: implements Eq. (3.2) of CEZ paper

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'plabel.f'
      logical:: passed
      integer:: j,maxparts,found,j1,j2
      real(dp):: pt,etarap,p(mxpart,4),ptj,ptj1,ptj2,
     & etaj1,etaj2

      passed=.false.

***************************** START CUTS *******************************

      found=0
      ptj1=zip
      ptj2=zip

c--- check for 2 tagging jets (highest pt)
      do j=3,maxparts
        if (  (plabel(j) == 'pp') .or. (plabel(j) == 'qj')
     &   .or. (plabel(j) == 'bq') .or. (plabel(j) == 'ba')) then
          ptj=pt(j,p)
          if (ptj > min(ptj1,ptj2)) then
            found=found+1
            if (found == 1) then
              ptj1=ptj
              j1=j
            endif
            if (found >= 2) then
              if (ptj > ptj1) then
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

      if (found < 2) goto 999

      etaj1=etarap(j1,p)
      etaj2=etarap(j2,p)

c--- ensure a rapidity gap of at least 4.2 between the tagging jets
      if (abs(etaj1-etaj2) < 4.2_dp) goto 999

c--- ensure the tagging jets lie in opposite hemispheres
      if (etaj1*etaj2 >= zip) goto 999

c      mj1j2=sqrt(max(zip,(p(j1,4)+p(j2,4))**2
c     & -(p(j1,1)+p(j2,1))**2-(p(j1,2)+p(j2,2))**2-(p(j1,3)+p(j2,3))**2))
c--- ensure the tagging jets have an invariant mass larger than 600 GeV
c--- NOT USED ANY LONGER
c      if (mj1j2 < 600._dp) goto 999

********************** END OF CUT CROSS SECTIONS ***********************

      passed=.true.

c--- this is the alternate return if the cuts are failed
  999 continue

      return

      end


