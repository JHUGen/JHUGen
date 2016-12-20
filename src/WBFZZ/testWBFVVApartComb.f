      subroutine testWBFVVApartComb(j1,j2,j7,j8,partSwapOk)
      implicit none
      include 'constants.f'
      include 'pid_pdg.f'
      integer j1,j2,j7,j8
      logical partSwapOk

      partSwapOk=.true.
      if (
     & (
     & (j1.eq.7 .or. j2.eq.7) .and. (
     &      (pid_pdg(7).ge.1 .and. pid_pdg(7).le.5)
     & .or. (pid_pdg(7) .eq. 12)
     & .or. (pid_pdg(7) .eq. 14)
     & .or. (pid_pdg(7) .eq. 16)
     & .or. (pid_pdg(7) .eq. 11)
     & .or. (pid_pdg(7) .eq. 13)
     & .or. (pid_pdg(7) .eq. 15)
     & )
     & ) .or. (
     & (j1.eq.8 .or. j2.eq.8) .and. (
     &      (pid_pdg(8).ge.1 .and. pid_pdg(8).le.5)
     & .or. (pid_pdg(8) .eq. 12)
     & .or. (pid_pdg(8) .eq. 14)
     & .or. (pid_pdg(8) .eq. 16)
     & .or. (pid_pdg(8) .eq. 11)
     & .or. (pid_pdg(8) .eq. 13)
     & .or. (pid_pdg(8) .eq. 15)
     & )
     & ) .or. (
     & (j7.eq.7 .or. j8.eq.7) .and. (
     &      (-pid_pdg(7).ge.1 .and. -pid_pdg(7).le.5)
     & .or. (pid_pdg(7) .eq. -12)
     & .or. (pid_pdg(7) .eq. -11)
     & .or. (pid_pdg(7) .eq. -13)
     & .or. (pid_pdg(7) .eq. -15)
     & )
     & ) .or. (
     & (j7.eq.8 .or. j8.eq.8) .and. (
     &      (-pid_pdg(8).ge.1 .and. -pid_pdg(8).le.5)
     & .or. (pid_pdg(8) .eq. -12)
     & .or. (pid_pdg(8) .eq. -11)
     & .or. (pid_pdg(8) .eq. -13)
     & .or. (pid_pdg(8) .eq. -15)
     & )
     & ) .or. (! Tests from the opposite direction
     & (j7.eq.1 .or. j8.eq.1) .and. (
     &      (pid_pdg(1).ge.1 .and. pid_pdg(1).le.5)
     & .or. (pid_pdg(1) .eq. 12)
     & .or. (pid_pdg(1) .eq. 14)
     & .or. (pid_pdg(1) .eq. 16)
     & .or. (pid_pdg(1) .eq. 11)
     & .or. (pid_pdg(1) .eq. 13)
     & .or. (pid_pdg(1) .eq. 15)
     & )
     & ) .or. (
     & (j7.eq.2 .or. j8.eq.2) .and. (
     &      (pid_pdg(2).ge.1 .and. pid_pdg(2).le.5)
     & .or. (pid_pdg(2) .eq. 12)
     & .or. (pid_pdg(2) .eq. 14)
     & .or. (pid_pdg(2) .eq. 16)
     & .or. (pid_pdg(2) .eq. 11)
     & .or. (pid_pdg(2) .eq. 13)
     & .or. (pid_pdg(2) .eq. 15)
     & )
     & ) .or. (
     & (j1.eq.1 .or. j2.eq.1) .and. (
     &      (-pid_pdg(1).ge.1 .and. -pid_pdg(1).le.5)
     & .or. (pid_pdg(1) .eq. -12)
     & .or. (pid_pdg(1) .eq. -11)
     & .or. (pid_pdg(1) .eq. -13)
     & .or. (pid_pdg(1) .eq. -15)
     & )
     & ) .or. (
     & (j1.eq.2 .or. j2.eq.2) .and. (
     &      (-pid_pdg(2).ge.1 .and. -pid_pdg(2).le.5)
     & .or. (pid_pdg(2) .eq. -12)
     & .or. (pid_pdg(2) .eq. -11)
     & .or. (pid_pdg(2) .eq. -13)
     & .or. (pid_pdg(2) .eq. -15)
     & )
     & ) ) then
      partSwapOk = .false.
      endif

      return
      end
