      subroutine genrad(p,i1,i2,i7,r1,r2,phit,wt,*)
c----this is only a switchyard routine.
c----i1 is the initial state vector, except for final-final case
c----i7 is the generated vector
      implicit none
      include 'constants.f'
      integer i1,i2,i7,j
      double precision p(mxpart,4),phit,wt
      double precision r1,r2

c---change sign of incoming momenta
      wt=0d0
      do j=1,4
      p(1,j)=-p(1,j)
      p(2,j)=-p(2,j)
      enddo

c----initial-initial case
      if ((i1 .le. 2) .and. (i2 .le. 2)) then
      call genrii(p,i1,i2,i7,r1,r2,phit,wt,*999)
c******************************************************************
c----final-final case
      elseif ((i1 .gt. 2) .and. (i2 .gt. 2)) then
      call genrff(p,i1,i2,i7,r1,r2,phit,wt,*999)
c------initial final
      elseif ((i1 .le. 2) .and. (i2 .gt. 2)) then
      call genrif(p,i1,i2,i7,r1,r2,phit,wt,*999)
c---protection 
      else
      write(6,*) 'Unimplemented case'
      stop
      endif

c---reverse signs of incoming momenta so that they are all incoming
      do j=1,4
      p(1,j)=-p(1,j)
      p(2,j)=-p(2,j)
      enddo

      return
 999  continue
      wt=0
      return 1

      end

