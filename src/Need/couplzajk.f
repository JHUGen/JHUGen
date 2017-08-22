      subroutine couplzajk()
      implicit none
      include 'constants.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'pid_pdg.f'
      include 'zacouplejk.f'
      integer n1,n2,j
      integer aid(mxpart)

      do j=1,mxpart
         aid(j) = abs(pid_pdg(j))
      enddo

      ! Find 17 and 28 couplings to Z/A
      ! The first two indices are for particle pairs,
      ! the third is for T3 down and up
      Q_jk(:,:,:)=0d0
      L_jk(:,:,:)=0d0
      R_jk(:,:,:)=0d0

      do n1=2,mxpart
      do n2=1,n1-1

      ! id1 unknown and id2 quark
      ! id2 unknown and id1 quark
      if (
     & (aid(n1).le.5 .and. aid(n2).le.5)
     & ) then
         do j=1,nf
            Q_jk(n1,n2,j) = Q(j)
            L_jk(n1,n2,j) = L(j)
            R_jk(n1,n2,j) = R(j)
         enddo
      else if (
     & (aid(n1).eq.11 .and. aid(n2).eq.11) .or.
     & (aid(n1).eq.12 .and. aid(n2).eq.12) .or.
     & (aid(n1).eq.13 .and. aid(n2).eq.13) .or.
     & (aid(n1).eq.14 .and. aid(n2).eq.14) .or.
     & (aid(n1).eq.15 .and. aid(n2).eq.15) .or.
     & (aid(n1).eq.16 .and. aid(n2).eq.16) .or.

     & (aid(n1).eq.11 .and. aid(n2).eq.12) .or.
     & (aid(n1).eq.13 .and. aid(n2).eq.14) .or.
     & (aid(n1).eq.15 .and. aid(n2).eq.16) .or.

     & (aid(n2).eq.11 .and. aid(n1).eq.12) .or.
     & (aid(n2).eq.13 .and. aid(n1).eq.14) .or.
     & (aid(n2).eq.15 .and. aid(n1).eq.16)

     & ) then
         Q_jk(n1,n2,1) = qe
         L_jk(n1,n2,1) = le
         R_jk(n1,n2,1) = re
         Q_jk(n1,n2,2) = 0d0
         L_jk(n1,n2,2) = ln
         R_jk(n1,n2,2) = rn
      endif

      ! Opposite pairing also has the same couplings
      do j=1,nf
         Q_jk(n2,n1,j) = Q_jk(n1,n2,j)
         L_jk(n2,n1,j) = L_jk(n1,n2,j)
         R_jk(n2,n1,j) = R_jk(n1,n2,j)
      enddo

      enddo
      enddo

      return
      end
