c--- Provide dummy routines so that the calls for n-tuple handling
c--- and event unweighting still exist even when the program is not
c--- compiled against CERNLIB.
c--- We could write our own routines here, but for now we just
c--- return an error message.

      subroutine dsw_error
      implicit none

      write(6,*) 'This version of MCFM has not been compiled against'
      write(6,*) 'CERNLIB, so output of n-tuples and event unweighting'
      write(6,*) 'are not available.'
      write(6,*)
      stop
      end
c      

      subroutine dswhbook(n,titlex,dx,xmin,xmax)
      implicit none
      integer n
      character titlex*8
      real*8 dx,xmin,xmax

      call dsw_error

      return
      end
c

      subroutine dswhfill(n,var,wgt)
      implicit none
      integer n
      real*8 var,wgt

      call dsw_error

      return
      end

c

      subroutine NTfinalize
      implicit none

      call dsw_error

      return
      end

c
      subroutine bookfill(tag,p,wt)
      implicit none
      include 'constants.f'
      include 'maxwt.f'

      character tag*4
      double precision p(mxpart,4)
      double precision wt 

      call dsw_error

      return
      end

c
      subroutine dswntuplebook
      implicit none

      call dsw_error

      return
      end

c
      subroutine dswntuplefill(p,wt)
      implicit none
      include 'constants.f'

      double precision p(mxpart,4)
      double precision wt 

      call dsw_error

      return
      end

c
      subroutine SORTZV(A,INDEX,N,MODE,NWAY,NSORT)
      implicit none
      include 'constants.f'
      include 'eventbuffer.f'
      real*4 A(buffersize)
      integer INDEX(buffersize),N,MODE,NWAY,NSORT
      
      call dsw_error

      return
      end

c
