      subroutine couplzajk()
      implicit none
      include 'constants.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'pid_pdg.f'
      include 'zacouplejk.f'
      integer n1,n2,j

      ! Find 17 and 28 couplings to Z/A
      Q_jk(:,:,:)=0d0
      L_jk(:,:,:)=0d0
      R_jk(:,:,:)=0d0

      do n1=2,mxpart
      do n2=1,n1-1

      if (
     & (
     & (abs(pid_pdg(n1)).ge.0 .and. abs(pid_pdg(n1)).le.5) .and.
     & (abs(pid_pdg(n2)).ge.0 .and. abs(pid_pdg(n2)).le.5) .and.
     & (abs(pid_pdg(n1)).eq.abs(pid_pdg(n2)))
     & )
     & .or.
     & (
     & (abs(pid_pdg(n1)).eq.0) .and.
     & (abs(pid_pdg(n2)).gt.0 .and. abs(pid_pdg(n2)).le.5)
     & )
     & .or.
     & (
     & (abs(pid_pdg(n2)).eq.0) .and.
     & (abs(pid_pdg(n1)).gt.0 .and. abs(pid_pdg(n1)).le.5)
     & )
     & ) then
         do j=1,nf
            Q_jk(n1,n2,j) = Q(j)
            L_jk(n1,n2,j) = L(j)
            R_jk(n1,n2,j) = R(j)
         enddo
      else if (
     & (abs(pid_pdg(n1)).ge.11 .and. abs(pid_pdg(n1)).le.16) .and.
     & (abs(pid_pdg(n2)).ge.11 .and. abs(pid_pdg(n2)).le.16) .and.
     & (abs(pid_pdg(n1)).eq.abs(pid_pdg(n2)))
     & ) then
         Q_jk(n1,n2,1) = qe
         L_jk(n1,n2,1) = le
         R_jk(n1,n2,1) = re
         Q_jk(n1,n2,2) = 0d0
         L_jk(n1,n2,2) = ln
         R_jk(n1,n2,2) = rn
      endif

      do j=1,nf
         Q_jk(n2,n1,j) = Q_jk(n1,n2,j)
         L_jk(n2,n1,j) = L_jk(n1,n2,j)
         R_jk(n2,n1,j) = R_jk(n1,n2,j)
      enddo

      enddo
      enddo

      return
      end
