      subroutine HAQggnew(p1,p2,p3,p4,ampsq,ampsq_ba,ampsq_ab,ampsq_sym)
c--- Note that the ab/ba ordering of the color-ordered matrix elements
c--- returned by this routine is OPOPSITE to that in hqqggdfm.f;
c--- the correct ordering is important for real subtraction terms
      implicit none
      include 'constants.f'
      integer p1,p2,p3,p4,j1,j2,j3
      double complex ab(2,2,2),ba(2,2,2)
      double precision ampsq,ampsq_ab,ampsq_ba,ampsq_sym

      call Amplo_AQgg(p1,p2,p3,p4,ab,ba)

c--- calculate the matrix element as the sum of two colour orderings,
c--- plus a colour-suppressed QED-like piece which is symmetric
c--- in the ordering of the two gluons
      ampsq_ab=0d0
      ampsq_ba=0d0
      ampsq_sym=0d0
      do j1=1,2
      do j2=1,2
      do j3=1,2
      ampsq_ab=ampsq_ab
     .  +cf*xn**2/2d0*cdabs(ab(j1,j2,j3))**2
      ampsq_ba=ampsq_ba
     .  +cf*xn**2/2d0*cdabs(ba(j1,j2,j3))**2
      ampsq_sym=ampsq_sym
     .  -cf/2d0*cdabs(ab(j1,j2,j3)+ba(j1,j2,j3))**2
      enddo
      enddo
      enddo

      ampsq=ampsq_ab+ampsq_ba+ampsq_sym

      return
      end



