
      subroutine store_zdip(nd,z)
      implicit none
      include 'types.f'
      
      include 'z_dip.f'
      real(dp):: z
      integer:: nd

      z_dip(nd)=z

      return 
      end subroutine
