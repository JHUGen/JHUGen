      subroutine ZZ_HZZ(p,msq)
      implicit none
      include 'types.f'

c--- Weak Bosion Fusion by Z-Z exchange only
c---Matrix element squared averaged over initial colors and spins
c
c     q(-p1)+q(-p2) -->  H(p3,p4)+q(p7)+q(p8)
c                           |
c                           |
c                           |
c                           ---> Z(e-(p3)+e+(p4))+Z(mu-(p5)+mu+(p6))
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: p(mxpart,4),fac,s3456
      real(dp):: msq(-nf:nf,-nf:nf),hdecay,
     & ud_ud_LL,udb_udb_LL,ud_ud_LR,udb_udb_LR

      msq(:,:)=zip

      call dotem(8,p,s)

      s3456=s(3,4)+s(3,5)+s(3,6)+s(4,5)+s(4,6)+s(5,6)
      hdecay=gwsq**3*zmass**2*4._dp*xw**2/(one-xw)*
     &  (((l1*l2)**2+(r1*r2)**2)*s(3,5)*s(4,6)
     & + ((r1*l2)**2+(r2*l1)**2)*s(3,6)*s(4,5))
      hdecay=hdecay/((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)
      hdecay=hdecay/((s(5,6)-zmass**2)**2+(zmass*zwidth)**2)
      hdecay=hdecay/((s3456-hmass**2)**2+(hmass*hwidth)**2)
      fac=0.25_dp*gwsq**3*hdecay
C Color cancels, 0.25_dp is spin average

C q-q and qbar-qbar
c--- u(1)+d(2) -> u(7)+d(8)
c--- ub(1)+db(2) -> ub(7)+db(8)
      call msqpieces_zz(1,2,7,8,ud_ud_LL,ud_ud_LR)

C q-qbar and qbar-q
c--- u(1)+db(2) -> u(7)+db(8)
c--- ub(1)+d(2) -> ub(7)+d(8)
      call msqpieces_zz(1,8,7,2,udb_udb_LL,udb_udb_LR)

      do j=-nf,nf
      do k=-nf,nf
        if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=fac*(
     &     +udb_udb_LL*((L(+j)*L(-k))**2+(R(+j)*R(-k))**2)
     &     +udb_udb_LR*((L(+j)*R(-k))**2+(R(+j)*L(-k))**2))
        elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=fac*(
     &     +udb_udb_LL*((L(-j)*L(k))**2+(R(-j)*R(k))**2)
     &     +udb_udb_LR*((L(-j)*R(k))**2+(R(-j)*L(k))**2))
        elseif ((j > 0) .and. (k > 0)) then
          msq(j,k)=fac*(
     &     +ud_ud_LL*((L(+j)*L(+k))**2+(R(+j)*R(+k))**2)
     &     +ud_ud_LR*((L(+j)*R(+k))**2+(R(+j)*L(+k))**2))
        elseif ((j < 0) .and. (k < 0)) then
          msq(j,k)=fac*(
     &     +ud_ud_LL*((L(-j)*L(-k))**2+(R(-j)*R(-k))**2)
     &     +ud_ud_LR*((L(-j)*R(-k))**2+(R(-j)*L(-k))**2))
        endif
      enddo
      enddo

      return
      end

