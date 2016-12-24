      subroutine couplza_1728(n1,n2,n7,n8)
      implicit none
      include 'constants.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'pid_pdg.f'
      include 'zacouple1728.f'
      integer n1,n2,n7,n8,j

      ! Find 17 and 28 couplings to Z/A
      Q17(:)=0d0;L17(:)=0d0;R17(:)=0d0;
      Q28(:)=0d0;L28(:)=0d0;R28(:)=0d0;
      if (
     & (abs(pid_pdg(n1)).ge.0 .and. abs(pid_pdg(n1)).le.5) .and.
     & (abs(pid_pdg(n7)).ge.0 .and. abs(pid_pdg(n7)).le.5)
     & ) then
         do j=1,nf
            Q17(j) = Q(j)
            L17(j) = L(j)
            R17(j) = R(j)
         enddo
      else if (
     & (abs(pid_pdg(n1)).ge.11 .and. abs(pid_pdg(n1)).le.16) .and.
     & (abs(pid_pdg(n7)).ge.11 .and. abs(pid_pdg(n7)).le.16)
     & ) then
         Q17(1) = qe
         L17(1) = le
         R17(1) = re
         Q17(2) = 0d0
         L17(2) = ln
         R17(2) = rn
      endif
      if (
     & (abs(pid_pdg(n2)).ge.0 .and. abs(pid_pdg(n2)).le.5) .and.
     & (abs(pid_pdg(n8)).ge.0 .and. abs(pid_pdg(n8)).le.5)
     & ) then
         do j=1,nf
            Q28(j) = Q(j)
            L28(j) = L(j)
            R28(j) = R(j)
         enddo
      else if (
     & (abs(pid_pdg(n2)).ge.11 .and. abs(pid_pdg(n2)).le.16) .and.
     & (abs(pid_pdg(n8)).ge.11 .and. abs(pid_pdg(n8)).le.16)
     & ) then
         Q28(1) = qe
         L28(1) = le
         R28(1) = re
         Q28(2) = 0d0
         L28(2) = ln
         R28(2) = rn
      endif

      return
      end
