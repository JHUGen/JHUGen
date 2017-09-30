      subroutine ttggHdriver(q,ampsq1,ampsq2,ampsq0)
      implicit none
      include 'constants.f'
      double precision q(mxpart,4),ampsq1,ampsq2,ampsq0
      integer h1,h2,h3,h4
      double complex
     & ampAB(2,2,2,2),ampBA(2,2,2,2),ampQED(2,2,2,2)
      call ttggHamp(q,3,4,1,2,1,1,ampAB)
      call ttggHamp(q,3,4,2,1,1,1,ampBA)
      ampsq1=0d0
      ampsq2=0d0
      ampsq0=0d0
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      ampQED(h1,h2,h3,h4)=ampAB(h1,h2,h3,h4)+ampBA(h1,h2,h4,h3)
      ampsq1=ampsq1
     &  +dble(ampAB(h1,h2,h3,h4)*Dconjg(ampAB(h1,h2,h3,h4)))
      ampsq2=ampsq2
     &  +dble(ampBA(h1,h2,h4,h3)*Dconjg(ampBA(h1,h2,h4,h3)))
      ampsq0=ampsq0
     & -dble(ampQED(h1,h2,h3,h4)*Dconjg(ampQED(h1,h2,h3,h4)))/xn**2
      enddo
      enddo
      enddo
      enddo
      return
      end
