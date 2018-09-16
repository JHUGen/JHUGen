c--- Provide dummy routines so that the calls for n-tuple handling
c--- and event unweighting still exist even when the program is not
c--- compiled against CERNLIB.
c--- We could write our own routines here, but for now we just
c--- return an error message.

      subroutine dsw_error
      implicit none
      include 'types.f'


      write(6,*) 'This version of MCFM has not been compiled against'
      write(6,*) 'CERNLIB, so output of n-tuples and event unweighting'
      write(6,*) 'are not available.'
      write(6,*)
      stop
      end
c

      subroutine dswhbook(n,titlex,dx,xmin,xmax)
      implicit none
      include 'types.f'

      integer:: n
      character titlex*8
      real(dp):: dx,xmin,xmax

      call dsw_error

      return
      end
c

      subroutine dswhfill(n,var,wgt)
      implicit none
      include 'types.f'

      integer:: n
      real(dp):: var,wgt

      call dsw_error

      return
      end

c

      subroutine NTfinalize
      implicit none
      include 'types.f'


      call dsw_error

      return
      end

c
      subroutine bookfill(tag,p,wt)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'maxwt.f'

      integer tag
      real(dp):: p(mxpart,4)
      real(dp):: wt

      call dsw_error

      return
      end

c
      subroutine dswntuplebook
      implicit none
      include 'types.f'


      call dsw_error

      return
      end

c
      subroutine dswntuplefill(p,wt)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'

      real(dp):: p(mxpart,4)
      real(dp):: wt

      call dsw_error

      return
      end

c
!      subroutine SORTZV(A,INDEX,N,MODE,NWAY,NSORT)
!      implicit none
!      include 'types.f'
!
!      include 'constants.f'
!      include 'nf.f'
!      include 'mxpart.f'
!      include 'cplx.h'
!      include 'eventbuffer.f'
!      real*4 A(buffersize)
!      integer:: INDEX(buffersize),N,MODE,NWAY,NSORT
!
!      call dsw_error
!
!      return
!      end

c
