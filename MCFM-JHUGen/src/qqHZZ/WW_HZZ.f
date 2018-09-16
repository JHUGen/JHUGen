      subroutine WW_HZZ(p,msq)
      implicit none
      include 'types.f'

c--- Weak Bosion Fusion by W-W exchange only
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
     & ud_du,uub_ddb

      integer,parameter::pn(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

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
c--- u(1)+d(2) -> d(7)+u(8)
c--- ub(1)+db(2) -> db(7)+ub(8)
      call msqpieces_ww(1,2,8,7,ud_du)

C q-qbar and qbar-q
c--- u(1)+ub(2) -> d(7)+db(8)
c--- ub(1)+d(2) -> db(7)+u(8)
      call msqpieces_ww(1,8,2,7,uub_ddb)

c--- Only loop up to (nf-1) to avoid b->t transitions
      do j=-(nf-1),nf-1
      do k=-(nf-1),nf-1
      msq(j,k)=0._dp
        if     ((j > 0) .and. (k < 0)) then
          if (pn(j) == -pn(k)) msq(j,k)=fac*uub_ddb
        elseif ((j < 0) .and. (k > 0)) then
          if (pn(j) == -pn(k)) msq(j,k)=fac*uub_ddb
        elseif ((j > 0) .and. (k > 0)) then
          if (pn(j)+pn(k) == +3) msq(j,k)=fac*ud_du
        elseif ((j < 0) .and. (k < 0)) then
          if (pn(j)+pn(k) == -3) msq(j,k)=fac*ud_du
        endif
      enddo
      enddo

      return
      end

