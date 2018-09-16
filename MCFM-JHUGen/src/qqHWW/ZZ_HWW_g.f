      subroutine ZZ_HWW_g(p,msq)
      implicit none
      include 'types.f'

c--- Weak Boson Fusion by Z-Z exchange only
c---Matrix element squared averaged over initial colors and spins
c
c     q(-p1)+q(-p2) -->  H(p3,p4)+q(p5)+q(p6)+g(p7)
c                           |
c                           |
c                           |
c                           ---> W+(nu(p3)+e+(p4))+W-(e-(p5)+nub(p6))
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      integer:: j,k,nup,ndo
      real(dp):: p(mxpart,4),facqq,facqg,s3456
      real(dp):: msq(-nf:nf,-nf:nf),hdecay,
     & ud_udg_LL,udb_udbg_LL,ud_udg_LR,udb_udbg_LR,
     & ug_uddb_LL,ug_uddb_LR,gu_ddbu_LL,gu_ddbu_LR

      parameter (nup=2,ndo=3)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call dotem(9,p,s)

      s3456=s(3,4)+s(3,5)+s(3,6)+s(4,5)+s(4,6)+s(5,6)
      hdecay=gwsq**3*wmass**2*s(3,5)*s(6,4)
      hdecay=hdecay/(((s3456-hmass**2)**2+(hmass*hwidth)**2)
     &   *((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
     &   *((s(5,6)-wmass**2)**2+(wmass*wwidth)**2))
      facqq=0.25_dp*gsq*Cf*gwsq**3*hdecay
      facqg=-facqq*(aveqg/aveqq)
C Extra factor of gsq*Cf compared to lowest order

C q-q and qbar-qbar
c--- u(1)+d(2) -> u(7)+d(8)+g(p9)
c--- ub(1)+db(2) -> ub(7)+db(8)+g(p9)
      call msq_gpieces_zz(1,2,7,8,9,ud_udg_LL,ud_udg_LR)

C q-qbar and qbar-q
c--- u(1)+db(2) -> u(7)+db(8)+g(p9)
c--- ub(1)+d(2) -> ub(7)+d(8)+g(p9)
      call msq_gpieces_zz(1,8,7,2,9,udb_udbg_LL,udb_udbg_LR)

C q-g and qbar-g
c--- u(1)+g(2) -> u(7)+[q](8)+[qb](9)
c--- ub(1)+g(2) -> ub(7)+[qb](8)+[q](9)
      call msq_gpieces_zz(1,9,7,8,2,ug_uddb_LL,ug_uddb_LR)

C g-q and g-qbar
c--- g(1)+u(2) -> [q](7)+[qb](8)+u(9)
c--- g(1)+ub(2) -> [qb](7)+[q](8)+ub(9)
      call msq_gpieces_zz(8,2,7,9,1,gu_ddbu_LL,gu_ddbu_LR)

      do j=-nf,nf
      do k=-nf,nf
        if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=facqq*(
     &     +udb_udbg_LL*((L(+j)*L(-k))**2+(R(+j)*R(-k))**2)
     &     +udb_udbg_LR*((L(+j)*R(-k))**2+(R(+j)*L(-k))**2))
        elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=facqq*(
     &     +udb_udbg_LL*((L(-j)*L(k))**2+(R(-j)*R(k))**2)
     &     +udb_udbg_LR*((L(-j)*R(k))**2+(R(-j)*L(k))**2))
        elseif ((j > 0) .and. (k > 0)) then
          msq(j,k)=facqq*(
     &     +ud_udg_LL*((L(+j)*L(+k))**2+(R(+j)*R(+k))**2)
     &     +ud_udg_LR*((L(+j)*R(+k))**2+(R(+j)*L(+k))**2))
        elseif ((j < 0) .and. (k < 0)) then
          msq(j,k)=facqq*(
     &     +ud_udg_LL*((L(-j)*L(-k))**2+(R(-j)*R(-k))**2)
     &     +ud_udg_LR*((L(-j)*R(-k))**2+(R(-j)*L(-k))**2))
        elseif ((j .ne. 0) .and. (k == 0)) then
          msq(j,k)=facqg*real(nup,dp)*(
     &     +ug_uddb_LL*((L(abs(j))*L(+2))**2+(R(abs(j))*R(+2))**2)
     &     +ug_uddb_LR*((L(abs(j))*R(+2))**2+(R(abs(j))*L(+2))**2))
     &            +facqg*real(ndo,dp)*(
     &     +ug_uddb_LL*((L(abs(j))*L(+1))**2+(R(abs(j))*R(+1))**2)
     &     +ug_uddb_LR*((L(abs(j))*R(+1))**2+(R(abs(j))*L(+1))**2))
        elseif ((j == 0) .and. (k .ne. 0)) then
          msq(j,k)=facqg*real(nup,dp)*(
     &     +gu_ddbu_LL*((L(abs(k))*L(+2))**2+(R(abs(k))*R(+2))**2)
     &     +gu_ddbu_LR*((L(abs(k))*R(+2))**2+(R(abs(k))*L(+2))**2))
     &            +facqg*real(ndo,dp)*(
     &     +gu_ddbu_LL*((L(abs(k))*L(+1))**2+(R(abs(k))*R(+1))**2)
     &     +gu_ddbu_LR*((L(abs(k))*R(+1))**2+(R(abs(k))*L(+1))**2))
        endif
      enddo
      enddo

      return
      end


