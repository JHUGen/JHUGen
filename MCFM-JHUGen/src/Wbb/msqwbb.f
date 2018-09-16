      double precision function msqwbb(i1,i2,i5,i6)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
C---This is the matrix element squared for 
c    |q_L(p1)+q_L(p5) --> q_L(p2)+q_L(pfive)+W(l(p6)+antilepton(p7))|^2
c   +|q_L(p1)+q_R(p5) --> q_L(p2)+q_R(pfive)+W(l(p6)+antilepton(p7))|^2
c---with couplings for W included
      integer i1,i2,i5,i6
      double complex aqqb_wbb
      double precision faclo
      faclo=V*gsq**2*gwsq**2*aveqq
      msqwbb=faclo
     &*(cdabs(aqqb_wbb(i1,i2,i5,i6,3,4))**2
     & +cdabs(aqqb_wbb(i1,i2,i6,i5,3,4))**2)
      return 
      end
