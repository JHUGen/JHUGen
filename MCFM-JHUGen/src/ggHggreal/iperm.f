      subroutine iperm(IHEL,PERM,IHELX,NGLUONS)
      implicit none
      include 'types.f'

      integer::J,NGLUONS,PERM(NGLUONS),IHEL(NGLUONS),IHELX(NGLUONS)
C permutes helicities to match momenta
      do j=1,ngluons
      ihelx(j)=IHEL(PERM(j))
      enddo
      return
      end
