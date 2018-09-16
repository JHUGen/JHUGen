      program mcfm
      implicit none
      include 'types.f'
      include 'mpicommon.f'
      real(dp):: r,er
c
      rank=0
      size=1

      call mcfmsub(r,er)
c
      end
