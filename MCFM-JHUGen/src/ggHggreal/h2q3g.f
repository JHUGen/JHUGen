      subroutine h2q3g(p1,p2,p3,p4,p5,Hqaggg)
      implicit none
      include 'types.f'

!-----calculates the matrix element squared
!-----for q(p1)+q~(p2)-->g(p3)+g(p4)+g(p5)+H
!----- using the results of S. Badger and
!----   V.~Del Duca, A.~Frizzo and F.~Maltoni,
!----   %``Higgs boson production in association with three jets,''
!----   JHEP {\bf 0405}, 064 (2004)
!----   [arXiv:hep-ph/0404013].
!-----has to be preceded by a call to spinoru to set up the za,zb
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      real(dp)::Hqaggg,sign,factor
      complex(dp)::qamp(6,2,2,2,2),temp(6,2,2,2,2)
!      complex(dp)::otmp(6,2,2,2,2)
      complex(dp)::A2q3g_mpppp,A2q3g_mmmpp,A2q3g_mmpmp,
     & A2q3g_mpmmp,
     & na2q3g_mmmpp,na2q3g_mmpmp,na2q3g_mpmmp,na2q3g_mpppp
      complex(dp)::a0hqbqggg_mp_ppp,a0hqbqggg_mp_mpp,
     & a0hqbqggg_mpppm,a0hqbqggg_mp_pmp


      real(dp)::xa,xb,xc,xd
!--- Full result
      parameter(xa=xn*cf**2,xb=-cf/two,xc=0.25_dp/xn,xd=(xn**2+one)*xc)

!--- Leading in colour
c    parameter(xa=xn**3/four,xb=zip,xc=zip,xd=zip)
!--- 1/N^2 suppressed
c    parameter(xa=-xn/two,xb=-xn/four,xc=zip,xd=xn/four)
!--- 1/N^4 suppressed
c   parameter(xa=one/four/xn,xb=one/four/xn,xc=one/four/xn,xd=one/four/xn)
!--- Note that the definition of matrix DJK here differs from the one
!--- in (B.22) of the paper by an overall factor of Cf which is
!--- restored in the sum below
      integer::j,k,h1,h3,h4,h5,p1,p2,p3,p4,p5,n(5)
      integer, parameter, dimension(6) ::i3=(/3,3,4,4,5,5/)
      integer, parameter, dimension(6) ::i4=(/4,5,3,5,3,4/)
      integer, parameter, dimension(6) ::i5=(/5,4,5,3,4,3/)

      real(dp), parameter ::DJK(6,6)=reshape(
     & (/xa,xb,xb,xc,xc,xd,
     &   xb,xa,xc,xd,xb,xc,
     &   xb,xc,xa,xb,xd,xc,
     &   xc,xd,xb,xa,xc,xb,
     &   xc,xb,xd,xc,xa,xb,
     &   xd,xc,xc,xb,xb,xa/),(/6,6/))
      include 'cplx.h'

      n(1)=p1
      n(2)=p2
      n(3)=p3
      n(4)=p4
      n(5)=p5

!--- To be used when relating amplitudes by symmetry:
!--- additional sign for the crossed amplitudes with initial gluon,
!--- but no sign for two gluons
      if ((p4 == 4) .and. (p3 < 3)) then
        sign=-one
      else
        sign=+one
      endif

!----definition of helicities is for outgoing lines
!----labelling is as follows
!     temp(j,h1,h3,h4,h5) since h2 can be obtained from h1

      do j=1,6

!      otmp(j,2,2,2,2)=A2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j)),
!     & za,zb)
      temp(j,2,2,2,2)=-nA2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j))
     & ,za,zb)

!--------
!      write(6,*) '2q3g:o:+- +++',
!     & A2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j)),za,zb)
!      write(6,*) '2q3g:n:+- +++',
!     & -nA2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j)),za,zb)
!--------

!      otmp(j,2,1,1,1)=A2q3g_mpppp(n(2),n(1),n(i5(j)),n(i4(j)),n(i3(j)),
!     & zb,za)
      temp(j,2,1,1,1)=-nA2q3g_mpppp(n(2),n(1),n(i5(j)),n(i4(j)),n(i3(j))
     & ,zb,za)

!      otmp(j,2,2,2,1) = A2q3g_mmmpp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
!     & n(2),zb,za)
      temp(j,2,2,2,1) =-nA2q3g_mmmpp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
     & n(2),zb,za)

!----------
!      write(6,*) '2q3g:o:+- ++-',
!     & A2q3g_mmmpp(n(2),n(i3(j)),n(i4(j)),n(i5(j)),n(1),zb,za)
!      write(6,*) '2q3g:n:+- ++-',
!     & -nA2q3g_mmmpp(n(2),n(i3(j)),n(i4(j)),n(i5(j)),n(1),zb,za)
!------------

!      otmp(j,2,2,1,1) = A2q3g_mmmpp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
!     & n(1),za,zb)
      temp(j,2,2,1,1) =-nA2q3g_mmmpp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     & n(1),za,zb)

!      otmp(j,2,2,1,2) = A2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j))
!     & ,n(2),zb,za)
      temp(j,2,2,1,2) = -nA2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
     & n(2),zb,za)

!--------
!      write(6,*) '2q3g:o:+- +-+',
!     & A2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),n(2),zb,za)
!      write(6,*) '2q3g:n:+- +-+',
!     & -nA2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),n(2),zb,za)
!--------

!      otmp(j,2,1,2,1) = A2q3g_mmpmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
!     & n(1),za,zb)
      temp(j,2,1,2,1) = -nA2q3g_mmpmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     & n(1),za,zb)

!      otmp(j,2,1,2,2) = A2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
!     & n(2),zb,za)
      temp(j,2,1,2,2) = -nA2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
     & n(2),zb,za)

!--------
!      write(6,*) '2q3g:o:+- -++',
!     & A2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),n(2),zb,za)
!      write(6,*) '2q3g:n:+- -++',
!     & -nA2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),n(2),zb,za)
!---------

!      otmp(j,2,1,1,2) = A2q3g_mpmmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
!     .n(1),za,zb)
      temp(j,2,1,1,2) = -nA2q3g_mpmmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     & n(1),za,zb)

!       pause

!      temp(j,1,1,1,1)=A2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j)),
!     & zb,za)
!      temp(j,1,2,2,2)=A2q3g_mpppp(n(2),n(1),n(i5(j)),n(i4(j)),n(i3(j)),
!     & za,zb)

!      temp(j,1,2,1,1) = A2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
!     & n(2),za,zb)
!      temp(j,1,2,2,1) = A2q3g_mpmmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
!     & n(1),zb,za)

!      temp(j,1,1,1,2) = A2q3g_mmmpp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
!     & n(2),za,zb)
!      temp(j,1,1,2,2) = A2q3g_mmmpp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
!     & n(1),zb,za)

!      temp(j,1,1,2,1)= A2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
!     & n(2),za,zb)
!      temp(j,1,2,1,2)= A2q3g_mmpmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
!     & n(1),zb,za)

!--- fastest to obtain remaining amplitudes by symmetry
!      otmp(j,1,1,1,1)=-sign*conjg(otmp(j,2,2,2,2))
!      otmp(j,1,2,2,2)=-sign*conjg(otmp(j,2,1,1,1))
!      otmp(j,1,2,1,1)=-sign*conjg(otmp(j,2,1,2,2))
!      otmp(j,1,2,2,1)=-sign*conjg(otmp(j,2,1,1,2))
!      otmp(j,1,1,1,2)=-sign*conjg(otmp(j,2,2,2,1))
!      otmp(j,1,1,2,2)=-sign*conjg(otmp(j,2,2,1,1))
!      otmp(j,1,1,2,1)=-sign*conjg(otmp(j,2,2,1,2))
!      otmp(j,1,2,1,2)=-sign*conjg(otmp(j,2,1,2,1))

      temp(j,1,1,1,1)=-sign*conjg(temp(j,2,2,2,2))
      temp(j,1,2,2,2)=-sign*conjg(temp(j,2,1,1,1))
      temp(j,1,2,1,1)=-sign*conjg(temp(j,2,1,2,2))
      temp(j,1,2,2,1)=-sign*conjg(temp(j,2,1,1,2))
      temp(j,1,1,1,2)=-sign*conjg(temp(j,2,2,2,1))
      temp(j,1,1,2,2)=-sign*conjg(temp(j,2,2,1,1))
      temp(j,1,1,2,1)=-sign*conjg(temp(j,2,2,1,2))
      temp(j,1,2,1,2)=-sign*conjg(temp(j,2,1,2,1))

      enddo

!      do j=1,6
!      do h1=1,2
!      do h3=1,2
!      do h4=1,2
!      do h5=1,2
!      write(6,*) j,h1,h3,h4,h5,temp(j,h1,h3,h4,h5)/
!     & otmp(j,h1,h3,h4,h5)
!      enddo
!      enddo
!      enddo
!      enddo
!      enddo
!      pause


!----At this stage we have setup the amplitudes but failed
!----to assign the helicities properly. So we now reshuffle
!----to get these right.
      do h1=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
      do j=1,6
      if (j== 1) qamp(j,h1,h3,h4,h5)=temp(j,h1,h3,h4,h5)
      if (j== 2) qamp(j,h1,h3,h4,h5)=temp(j,h1,h3,h5,h4)
      if (j== 3) qamp(j,h1,h3,h4,h5)=temp(j,h1,h4,h3,h5)
      if (j== 4) qamp(j,h1,h3,h4,h5)=temp(j,h1,h4,h5,h3)
      if (j== 5) qamp(j,h1,h3,h4,h5)=temp(j,h1,h5,h3,h4)
      if (j== 6) qamp(j,h1,h3,h4,h5)=temp(j,h1,h5,h4,h3)
      enddo
      enddo
      enddo
      enddo
      enddo

!--- now perform the sum with the appropriate weights,
!--- c.f. Eq. (B.20); note that the factor of (gsq)**3
!--- is included in the wrapping routine

!--- NB: use symmetry to slightly improve speed, summing over
!---     diagonal and above in matrix
      Hqaggg=zip
      do j=1,6
      do k=j,6
      do h1=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
      if (j == k) then
      factor=one
      else
      factor=two
      endif
      Hqaggg=Hqaggg+two*Cf*djk(j,k)
     & *real(qamp(j,h1,h3,h4,h5)*conjg(qamp(k,h1,h3,h4,h5)))*factor
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
