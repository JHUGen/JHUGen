      double complex function t2(j1,j2,j3,j4)
      implicit none
      include 'constants.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4
      t2=zb(j1,j2)*za(j4,j2)+zb(j1,j3)*za(j4,j3)
      return
      end

