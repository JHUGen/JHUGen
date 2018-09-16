      subroutine qqb_ZH1jet(p,msq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  H  + Z +g(p7)
c                           |    |
c                           |    --> fermion(p3)+antifermion(p4)
c                           |
c                           ---> b(p5)+b(p6)
c   for the moment --- radiation only from initial line

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'hbbparams.f'
      include 'hdecaymode.f'
      integer:: j,k
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: radiLL,radiLL_ww
      real(dp):: qqbZHgL,qbqZHgL,qgZHqL,gqZHqL,gqbZHqbL,qbgZHqbL
      real(dp):: qqbZHgR,qbqZHgR,qgZHqR,gqZHqR,gqbZHqbR,qbgZHqbR


      msq(:,:)=zip

      if (hdecaymode=='wpwm') then
         call dotem(9,p,s)
         qqbZHgL=aveqq*radiLL_ww(1,2,9,5,6,7,8,3,4)
         qqbZHgR=aveqq*radiLL_ww(1,2,9,5,6,7,8,4,3)
         qbqZHgL=aveqq*radiLL_ww(2,1,9,5,6,7,8,3,4)
         qbqZHgR=aveqq*radiLL_ww(2,1,9,5,6,7,8,4,3)

         qgZHqL=-radiLL_ww(1,9,2,5,6,7,8,3,4)*aveqg
         qgZHqR=-radiLL_ww(1,9,2,5,6,7,8,4,3)*aveqg
         gqZHqL=-radiLL_ww(2,9,1,5,6,7,8,3,4)*aveqg
         gqZHqR=-radiLL_ww(2,9,1,5,6,7,8,4,3)*aveqg

         gqbZHqbL=-radiLL_ww(9,2,1,5,6,7,8,3,4)*aveqg
         gqbZHqbR=-radiLL_ww(9,2,1,5,6,7,8,4,3)*aveqg

         qbgZHqbL=-radiLL_ww(9,1,2,5,6,7,8,3,4)*aveqg
         qbgZHqbR=-radiLL_ww(9,1,2,5,6,7,8,4,3)*aveqg
      else
         call dotem(7,p,s)

         qqbZHgL=aveqq*radiLL(1,2,7,5,6,3,4)
         qqbZHgR=aveqq*radiLL(1,2,7,5,6,4,3)
         qbqZHgL=aveqq*radiLL(2,1,7,5,6,3,4)
         qbqZHgR=aveqq*radiLL(2,1,7,5,6,4,3)

         qgZHqL=-radiLL(1,7,2,5,6,3,4)*aveqg
         qgZHqR=-radiLL(1,7,2,5,6,4,3)*aveqg
         gqZHqL=-radiLL(2,7,1,5,6,3,4)*aveqg
         gqZHqR=-radiLL(2,7,1,5,6,4,3)*aveqg

         gqbZHqbL=-radiLL(7,2,1,5,6,3,4)*aveqg
         gqbZHqbR=-radiLL(7,2,1,5,6,4,3)*aveqg

         qbgZHqbL=-radiLL(7,1,2,5,6,3,4)*aveqg
         qbgZHqbR=-radiLL(7,1,2,5,6,4,3)*aveqg
      endif

      do j=-nf,nf
      do k=-nf,nf

      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) cycle
      if ((j == 0) .and. (k == 0)) cycle

      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=((L(j)*l1)**2+(R(j)*r1)**2)*qqbZHgL
     &            +((R(j)*l1)**2+(L(j)*r1)**2)*qqbZHgR
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=((L(k)*l1)**2+(R(k)*r1)**2)*qbqZHgL
     &            +((R(k)*l1)**2+(L(k)*r1)**2)*qbqZHgR
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=((L(j)*l1)**2+(R(j)*r1)**2)*qgZHqL
     &            +((R(j)*l1)**2+(L(j)*r1)**2)*qgZHqR
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=((L(-j)*l1)**2+(R(-j)*r1)**2)*qbgZHqbL
     &            +((R(-j)*l1)**2+(L(-j)*r1)**2)*qbgZHqbR
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=((L(k)*l1)**2+(R(k)*r1)**2)*gqZHqL
     &            +((R(k)*l1)**2+(L(k)*r1)**2)*gqZHqR
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=((L(-k)*l1)**2+(R(-k)*r1)**2)*gqbZHqbL
     &            +((R(-k)*l1)**2+(L(-k)*r1)**2)*gqbZHqbR
      endif

      enddo
      enddo
c---  adjust for fixed H->bb BR if necessary
      if ((hdecaymode == 'bqba') .and. (FixBrHbb)) then
         msq(:,:)=msq(:,:)*GamHbb/GamHbb0
      endif


      return
      end


