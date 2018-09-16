      subroutine qqb_zh_ww_g(P,msq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  H  + Z +g(p7)
c                           |    |
c                           |    --> fermion(p3)+antifermion(p4)
c                           |
c                           ---> W^+(nu(p5),e^+(p6))W^-(e^-(p7),nub(p8))
c   for the moment --- radiation only from initial line

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zcouple.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: radiLL_ww
      real(dp):: qqbZHgL,qbqZHgL,qgZHqL,gqZHqL,gqbZHqbL,qbgZHqbL
      real(dp):: qqbZHgR,qbqZHgR,qgZHqR,gqZHqR,gqbZHqbR,qbgZHqbR


      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

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

      do j=-nf,nf
      do k=-nf,nf

      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 40
      if ((j == 0) .and. (k == 0)) goto 40

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
 40   continue
      enddo
      enddo
      return
      end


      function radiLL_ww(j1,j2,j3,j4,j5,j6,j7,j8,j9)
      implicit none
      include 'types.f'
      real(dp):: radiLL_ww

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      integer:: j1,j2,j3,j4,j5,j6,j7,j8,j9
      real(dp):: s4567,s12,s13,s23,s123,prop
      real(dp):: fac,hdecay

      s4567=s(j4,j5)+s(j4,j6)+s(j4,j7)+s(j5,j6)+s(j5,j7)+s(j6,j7)
      s12=s(j1,j2)
      s13=s(j1,j3)
      s23=s(j2,j3)
      s123=s12+s13+s23
c---calculate the 2 Z propagators
      prop=       ((s123-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(j8,j9)-zmass**2)**2+(zmass*zwidth)**2)

      fac=8._dp*cf*xn*(xw/(1._dp-xw))**2*gsq*gwsq**3*wmass**2/prop
      hdecay=gwsq**3*wmass**2*s(j4,j6)*s(j7,j5)
      hdecay=hdecay/(((s4567-hmass**2)**2+(hmass*hwidth)**2)
     &   *((s(j4,j5)-wmass**2)**2+(wmass*wwidth)**2)
     &   *((s(j6,j7)-wmass**2)**2+(wmass*wwidth)**2))
      fac=fac*hdecay

      radiLL_ww=s12/s13/s23
     & *(2*s(j1,j9)*s(j2,j8)+s(j1,j9)*s(j3,j8)+s(j2,j8)*s(j3,j9))
     & +(s(j1,j9)*s(j2,j8)+s(j2,j8)*s(j3,j9)-s(j1,j8)*s(j1,j9))/s13
     & +(s(j1,j9)*s(j2,j8)+s(j1,j9)*s(j3,j8)-s(j2,j8)*s(j2,j9))/s23

      radiLL_ww=fac*radiLL_ww
      return
      end
