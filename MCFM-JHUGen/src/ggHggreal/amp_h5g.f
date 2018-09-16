      function amp_h5g(J,JHEL)
      IMPLICIT NONE
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      complex(dp)::amp_h5g
c***************************************************************************
c     this is a wrapper to generate all the helicity configurations for
c     higgs+ 5 gluons
c
c     given the basis:
c     app=+++++
c     amp=-++++
c     amm=--+++
c     amplitudes calculated by alberto frizzo and stored in amp_h_frizzo.f
c***************************************************************************
C     Based on routine of the same name by Del Duca et al,
c----- using the results of S. Badger and
c V.~Del Duca, A.~Frizzo and F.~Maltoni,
c ``Higgs boson production in association with three jets,''
c JHEP {\bf 0405}, 064 (2004)
c [arXiv:hep-ph/0404013].
c%%CITATION = HEP-PH 0404013;%%


c
C
C ARGUMENTS
C
      integer::JHEL(5),J(5)
C
C LOCAL VARIABLES
C
      logical::found
      integer::i,k,iconf,id
C
C EXTERNAL VARIABLES
C
      complex(dp)::appppp,nAmpppp,nAmmppp
      integer::NHEL(5,16),ihel
      include 'hels.f'
C checking
C      complex(dp)::amp_h5g_old

C-------
C BEGIN
C-------

c-- identify which helicity combination is needed
      found=.false.
        k=1
      do while (.not.found)
         id=0
        do i=1,5
c          write (*,*) 'jhel(',i,')=',jhel(i),'  nhel(',i,',',k,')=',nhel(i,k)
         if(jhel(i)==nhel(i,k)) id=id+1
          enddo
          if(id==5) then
         iconf=k
           found=.true.
          endif
        k=k+1
      enddo

c-- start the various helicity cases:

      if(iconf==1) then ! +++++
          amp_h5g=Appppp(J(1),J(2),J(3),J(4),J(5))

        elseif(iconf==2) then ! -++++
          amp_h5g=nAmpppp(J(1),J(2),J(3),J(4),J(5))

        elseif(iconf==3) then ! +-+++
          amp_h5g=nAmpppp(J(2),J(3),J(4),J(5),J(1))

        elseif(iconf==4) then ! ++-++
          amp_h5g=nAmpppp(J(3),J(4),J(5),J(1),J(2))

        elseif(iconf==5) then ! +++-+
          amp_h5g=nAmpppp(J(4),J(5),J(1),J(2),J(3))

        elseif(iconf==6) then ! ++++-
          amp_h5g=nAmpppp(J(5),J(1),J(2),J(3),J(4))

c--
! GZ replace call to Frizzo Amp dfm_amm with nAmmppp which uses MHV
! amplitudes (hep-th/0411092)
        elseif(iconf== 7) then ! --+++
          amp_h5g=nAmmppp(J(1),J(2),J(3),J(4),J(5))

        elseif(iconf== 8) then ! -+-++
           amp_h5g=-nAmmppp(J(1),J(3),J(2),J(4),J(5))
     &             -nAmmppp(J(1),J(3),J(4),J(2),J(5))
     &             -nAmmppp(J(1),J(3),J(4),J(5),J(2))

        elseif(iconf== 9) then ! -++-+
           amp_h5g=-nAmmppp(J(4),J(1),J(5),J(2),J(3))
     &             -nAmmppp(J(4),J(1),J(2),J(5),J(3))
     &             -nAmmppp(J(4),J(1),J(2),J(3),J(5))

        elseif(iconf== 10) then ! -+++-
          amp_h5g=nAmmppp(J(5),J(1),J(2),J(3),J(4))

        elseif(iconf== 11) then ! +--++
          amp_h5g=nAmmppp(J(2),J(3),J(4),J(5),J(1))

        elseif(iconf== 12) then ! +-+-+
           amp_h5g=-nAmmppp(J(2),J(4),J(3),J(5),J(1))
     &             -nAmmppp(J(2),J(4),J(5),J(3),J(1))
     &             -nAmmppp(J(2),J(4),J(5),J(1),J(3))

        elseif(iconf== 13) then ! +-++-
           amp_h5g=-nAmmppp(J(5),J(2),J(1),J(3),J(4))
     &             -nAmmppp(J(5),J(2),J(3),J(1),J(4))
     &             -nAmmppp(J(5),J(2),J(3),J(4),J(1))

        elseif(iconf== 14) then ! ++--+
          amp_h5g=nAmmppp(J(3),J(4),J(5),J(1),J(2))

        elseif(iconf== 15) then ! ++-+-
           amp_h5g=-nAmmppp(J(3),J(5),J(4),J(1),J(2))
     &             -nAmmppp(J(3),J(5),J(1),J(4),J(2))
     &             -nAmmppp(J(3),J(5),J(1),J(2),J(4))

        elseif(iconf== 16) then ! +++--
          amp_h5g=nAmmppp(J(4),J(5),J(1),J(2),J(3))


      else

      write(*,*) 'unknown helicity configuration'
      endif

! Check with old results by Frizzo & Co.
!        if(iconf== 7) then ! --+++
!          amp_h5g_old=Amm(J(1),J(2),J(3),J(4),J(5))
!
!          !write(*,*) 'amp_h5g',amp_h5g
!          !amp_h5g=nAmmppp(J(1),J(2),J(3),J(4),J(5))
!          !write(*,*) 'amp_h5g',amp_h5g
!          !pause
!
!        elseif(iconf== 8) then ! -+-++
!           amp_h5g_old=-Amm(J(1),J(3),J(2),J(4),J(5))
!     &             -Amm(J(1),J(3),J(4),J(2),J(5))
!     &             -Amm(J(1),J(3),J(4),J(5),J(2))
!
!        elseif(iconf== 9) then ! -++-+
!           amp_h5g_old=-Amm(J(4),J(1),J(5),J(2),J(3))
!     &             -Amm(J(4),J(1),J(2),J(5),J(3))
!     &             -Amm(J(4),J(1),J(2),J(3),J(5))
!
!        elseif(iconf== 10) then ! -+++-
!          amp_h5g_old=Amm(J(5),J(1),J(2),J(3),J(4))
!
!        elseif(iconf== 11) then ! +--++
!          amp_h5g_old=Amm(J(2),J(3),J(4),J(5),J(1))
!
!        elseif(iconf== 12) then ! +-+-+
!           amp_h5g_old=-Amm(J(2),J(4),J(3),J(5),J(1))
!     &             -Amm(J(2),J(4),J(5),J(3),J(1))
!     &             -Amm(J(2),J(4),J(5),J(1),J(3))
!
!        elseif(iconf== 13) then ! +-++-
!           amp_h5g_old=-Amm(J(5),J(2),J(1),J(3),J(4))
!     &             -Amm(J(5),J(2),J(3),J(1),J(4))
!     &             -Amm(J(5),J(2),J(3),J(4),J(1))
!
!        elseif(iconf== 14) then ! ++--+
!          amp_h5g_old=Amm(J(3),J(4),J(5),J(1),J(2))
!
!        elseif(iconf== 15) then ! ++-+-
!           amp_h5g_old=-Amm(J(3),J(5),J(4),J(1),J(2))
!     &             -Amm(J(3),J(5),J(1),J(4),J(2))
!     &             -Amm(J(3),J(5),J(1),J(2),J(4))
!
!        elseif(iconf== 16) then ! +++--
!          amp_h5g_old=Amm(J(4),J(5),J(1),J(2),J(3))
!       endif
!
!       if (iconf >=7 .and. abs(amp_h5g_old-amp_h5g) > 0.00one) then
!          write(*,*) 'amp_h5g:', iconf, amp_h5g_old,amp_h5g
!       endif


      return
      end

