      subroutine h2q3g(p1,p2,p3,p4,p5,Hqaggg)
      implicit none
C-----calculates the matrix element squared
C-----for q(p1)+q~(p2)-->g(p3)+g(p4)+g(p5)+H
c----- using the results of S. Badger and
C----   V.~Del Duca, A.~Frizzo and F.~Maltoni,
C----   %``Higgs boson production in association with three jets,''
C----   JHEP {\bf 0405}, 064 (2004)
C----   [arXiv:hep-ph/0404013].
C-----has to be preceded by a call to spinoru to set up the za,zb
      include 'constants.f'
      include 'zprods_com.f'
      double precision Hqaggg,sign,factor
      double complex qamp(6,2,2,2,2),temp(6,2,2,2,2)
c      double complex otmp(6,2,2,2,2)
      double complex A2q3g_mpppp,A2q3g_mmmpp, A2q3g_mmpmp, A2q3g_mpmmp,
     & na2q3g_mmmpp,na2q3g_mmpmp,na2q3g_mpmmp,na2q3g_mpppp
      double complex a0hqbqggg_mp_ppp,a0hqbqggg_mp_mpp,a0hqbqggg_mpppm,
     & a0hqbqggg_mp_pmp


      double precision DJK(6,6),xa,xb,xc,xd
c--- Full result
      parameter(xa=xn*cf**2,xb=-cf/2d0,xc=0.25d0/xn,xd=(xn**2+1d0)*xc)
c--- Leading in colour
c      parameter(xa=xn**3/4d0,xb=0d0,xc=0d0,xd=0d0)
c--- 1/N^2 suppressed
c      parameter(xa=-xn/2d0,xb=-xn/4d0,xc=0d0,xd=xn/4d0)
c--- 1/N^4 suppressed
c      parameter(xa=0.25d0/xn,xb=0.25d0/xn,xc=0.25d0/xn,xd=0.25d0/xn)
c--- Note that the definition of matrix DJK here differs from the one
c--- in (B.22) of the paper by an overall factor of Cf which is
c--- restored in the sum below
      integer j,k,h1,h3,h4,h5,p1,p2,p3,p4,p5,n(5)
      integer, parameter, dimension(6) :: i3=(/3,3,4,4,5,5/)
      integer, parameter, dimension(6) :: i4=(/4,5,3,5,3,4/)
      integer, parameter, dimension(6) :: i5=(/5,4,5,3,4,3/)

      DATA  (DJK(J,1),J=1,6)/xa,xb,xb,xc,xc,xd/
      DATA  (DJK(J,2),J=1,6)/xb,xa,xc,xd,xb,xc/
      DATA  (DJK(J,3),J=1,6)/xb,xc,xa,xb,xd,xc/
      DATA  (DJK(J,4),J=1,6)/xc,xd,xb,xa,xc,xb/
      DATA  (DJK(J,5),J=1,6)/xc,xb,xd,xc,xa,xb/
      DATA  (DJK(J,6),J=1,6)/xd,xc,xc,xb,xb,xa/

      save djk
!$omp threadprivate(djk)

      n(1)=p1
      n(2)=p2
      n(3)=p3
      n(4)=p4
      n(5)=p5

c--- To be used when relating amplitudes by symmetry:
c--- additional sign for the crossed amplitudes with initial gluon,
c--- but no sign for two gluons
      if ((p4 .eq. 4) .and. (p3 .lt. 3)) then
        sign=-1d0
      else
        sign=+1d0
      endif

C----definition of helicities is for outgoing lines
C----labelling is as follows
C     temp(j,h1,h3,h4,h5) since h2 can be obtained from h1

      do j=1,6

c      otmp(j,2,2,2,2)=A2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j)),
c     . za,zb)
      temp(j,2,2,2,2)=-nA2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j))
     & ,za,zb)

C--------
c      write(6,*) '2q3g:o:+- +++',
c     & A2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j)),za,zb)
c      write(6,*) '2q3g:n:+- +++',
c     & -nA2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j)),za,zb)
C--------

c      otmp(j,2,1,1,1)=A2q3g_mpppp(n(2),n(1),n(i5(j)),n(i4(j)),n(i3(j)),
c     . zb,za)
      temp(j,2,1,1,1)=-nA2q3g_mpppp(n(2),n(1),n(i5(j)),n(i4(j)),n(i3(j))
     & ,zb,za)

c      otmp(j,2,2,2,1) = A2q3g_mmmpp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
c     . n(2),zb,za)
      temp(j,2,2,2,1) =-nA2q3g_mmmpp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
     & n(2),zb,za)

C----------
c      write(6,*) '2q3g:o:+- ++-',
c     & A2q3g_mmmpp(n(2),n(i3(j)),n(i4(j)),n(i5(j)),n(1),zb,za)
c      write(6,*) '2q3g:n:+- ++-',
c     & -nA2q3g_mmmpp(n(2),n(i3(j)),n(i4(j)),n(i5(j)),n(1),zb,za)
c------------

c      otmp(j,2,2,1,1) = A2q3g_mmmpp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
c     . n(1),za,zb)
      temp(j,2,2,1,1) =-nA2q3g_mmmpp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     & n(1),za,zb)

c      otmp(j,2,2,1,2) = A2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j))
c     . ,n(2),zb,za)
      temp(j,2,2,1,2) = -nA2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
     & n(2),zb,za)

C--------
c      write(6,*) '2q3g:o:+- +-+',
c     & A2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),n(2),zb,za)
c      write(6,*) '2q3g:n:+- +-+',
c     & -nA2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),n(2),zb,za)
C--------

c      otmp(j,2,1,2,1) = A2q3g_mmpmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
c     . n(1),za,zb)
      temp(j,2,1,2,1) = -nA2q3g_mmpmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     & n(1),za,zb)

c      otmp(j,2,1,2,2) = A2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
c     . n(2),zb,za)
      temp(j,2,1,2,2) = -nA2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
     & n(2),zb,za)

c--------
c      write(6,*) '2q3g:o:+- -++',
c     & A2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),n(2),zb,za)
c      write(6,*) '2q3g:n:+- -++',
c     & -nA2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),n(2),zb,za)
c---------

c      otmp(j,2,1,1,2) = A2q3g_mpmmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
c     .n(1),za,zb)
      temp(j,2,1,1,2) = -nA2q3g_mpmmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     & n(1),za,zb)

c       pause

c      temp(j,1,1,1,1)=A2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j)),
c     . zb,za)
c      temp(j,1,2,2,2)=A2q3g_mpppp(n(2),n(1),n(i5(j)),n(i4(j)),n(i3(j)),
c     . za,zb)

c      temp(j,1,2,1,1) = A2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
c     . n(2),za,zb)
c      temp(j,1,2,2,1) = A2q3g_mpmmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
c     . n(1),zb,za)

c      temp(j,1,1,1,2) = A2q3g_mmmpp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
c     . n(2),za,zb)
c      temp(j,1,1,2,2) = A2q3g_mmmpp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
c     . n(1),zb,za)

c      temp(j,1,1,2,1)= A2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
c     . n(2),za,zb)
c      temp(j,1,2,1,2)= A2q3g_mmpmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
c     . n(1),zb,za)

c--- fastest to obtain remaining amplitudes by symmetry
c      otmp(j,1,1,1,1)=-sign*dconjg(otmp(j,2,2,2,2))
c      otmp(j,1,2,2,2)=-sign*dconjg(otmp(j,2,1,1,1))
c      otmp(j,1,2,1,1)=-sign*dconjg(otmp(j,2,1,2,2))
c      otmp(j,1,2,2,1)=-sign*dconjg(otmp(j,2,1,1,2))
c      otmp(j,1,1,1,2)=-sign*dconjg(otmp(j,2,2,2,1))
c      otmp(j,1,1,2,2)=-sign*dconjg(otmp(j,2,2,1,1))
c      otmp(j,1,1,2,1)=-sign*dconjg(otmp(j,2,2,1,2))
c      otmp(j,1,2,1,2)=-sign*dconjg(otmp(j,2,1,2,1))

      temp(j,1,1,1,1)=-sign*dconjg(temp(j,2,2,2,2))
      temp(j,1,2,2,2)=-sign*dconjg(temp(j,2,1,1,1))
      temp(j,1,2,1,1)=-sign*dconjg(temp(j,2,1,2,2))
      temp(j,1,2,2,1)=-sign*dconjg(temp(j,2,1,1,2))
      temp(j,1,1,1,2)=-sign*dconjg(temp(j,2,2,2,1))
      temp(j,1,1,2,2)=-sign*dconjg(temp(j,2,2,1,1))
      temp(j,1,1,2,1)=-sign*dconjg(temp(j,2,2,1,2))
      temp(j,1,2,1,2)=-sign*dconjg(temp(j,2,1,2,1))

      enddo

c      do j=1,6
c      do h1=1,2
c      do h3=1,2
c      do h4=1,2
c      do h5=1,2
c      write(6,*) j,h1,h3,h4,h5,temp(j,h1,h3,h4,h5)/
c     & otmp(j,h1,h3,h4,h5)
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo
c      pause


C----At this stage we have setup the amplitudes but failed
C----to assign the helicities properly. So we now reshuffle
C----to get these right.
      do h1=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
      do j=1,6
      if (j.eq. 1) qamp(j,h1,h3,h4,h5)=temp(j,h1,h3,h4,h5)
      if (j.eq. 2) qamp(j,h1,h3,h4,h5)=temp(j,h1,h3,h5,h4)
      if (j.eq. 3) qamp(j,h1,h3,h4,h5)=temp(j,h1,h4,h3,h5)
      if (j.eq. 4) qamp(j,h1,h3,h4,h5)=temp(j,h1,h4,h5,h3)
      if (j.eq. 5) qamp(j,h1,h3,h4,h5)=temp(j,h1,h5,h3,h4)
      if (j.eq. 6) qamp(j,h1,h3,h4,h5)=temp(j,h1,h5,h4,h3)
      enddo
      enddo
      enddo
      enddo
      enddo

c--- now perform the sum with the appropriate weights,
c--- c.f. Eq. (B.20); note that the factor of (gsq)**3
c--- is included in the wrapping routine

c--- NB: use symmetry to slightly improve speed, summing over
c---     diagonal and above in matrix
      Hqaggg=0d0
      do j=1,6
      do k=j,6
      do h1=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
      if (j .eq. k) then
      factor=one
      else
      factor=two
      endif
      Hqaggg=Hqaggg+2d0*Cf*djk(j,k)
     . *dble(qamp(j,h1,h3,h4,h5)*dconjg(qamp(k,h1,h3,h4,h5)))*factor
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
