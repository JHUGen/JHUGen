
      subroutine store_zdip(nd,z)
      implicit none
      include 'z_dip.f'
      double precision z
      integer nd

      z_dip(nd)=z

      return 
      end subroutine
