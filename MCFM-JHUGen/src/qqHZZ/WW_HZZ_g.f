      subroutine WW_HZZ_g(p,msq)
      implicit none
      include 'types.f'

c--- Weak Boson Fusion by W-W exchange only
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
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: p(mxpart,4),facqq,facqg,s3456
      real(dp):: msq(-nf:nf,-nf:nf),hdecay,
     & ud_dug,uub_ddbg,ug_dudb,gu_uddb

      integer,parameter::pn(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call dotem(9,p,s)

      s3456=s(3,4)+s(3,5)+s(3,6)+s(4,5)+s(4,6)+s(5,6)
      hdecay=gwsq**3*zmass**2*4._dp*xw**2/(one-xw)*
     &  (((l1*l2)**2+(r1*r2)**2)*s(3,5)*s(4,6)
     & + ((r1*l2)**2+(r2*l1)**2)*s(3,6)*s(4,5))
      hdecay=hdecay/((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)
      hdecay=hdecay/((s(5,6)-zmass**2)**2+(zmass*zwidth)**2)
      hdecay=hdecay/((s3456-hmass**2)**2+(hmass*hwidth)**2)
      facqq=0.25_dp*gsq*Cf*gwsq**3*hdecay
      facqg=-facqq*(aveqg/aveqq)
C Extra factor of gsq*Cf compared to lowest order

C q-q and qbar-qbar
c--- u(1)+d(2) -> d(7)+u(8)+g(p9)
c--- ub(1)+db(2) -> db(7)+ub(8)+g(p9)
      call msq_gpieces_ww(1,2,8,7,9,ud_dug)

C q-qbar and qbar-q
c--- u(1)+ub(2) -> d(7)+db(8)+g(p9)
c--- ub(1)+u(2) -> db(7)+d(8)+g(p9)
      call msq_gpieces_ww(1,8,2,7,9,uub_ddbg)

C q-g and qbar-g
c--- u(1)+g(2) -> d(7)+[u/c](8)+[db/sb](9)
c--- ub(1)+g(2) -> db(7)+[ub/cb](8)+[d/s](9)
      call msq_gpieces_ww(1,9,8,7,2,ug_dudb)

C g-q and g-qbar
c--- g(1)+u(2) -> [u/c](7)+[db/sb](8)+d(9)
c--- g(1)+ub(2) -> [ub/cb](7)+[d/s](8)+db(9)
      call msq_gpieces_ww(8,2,9,7,1,gu_uddb)

c--- Only loop up to (nf-1) to avoid b->t transitions
      do j=-(nf-1),nf-1
      do k=-(nf-1),nf-1
        if     ((j > 0) .and. (k < 0)) then
          if (pn(j) == -pn(k)) msq(j,k)=facqq*uub_ddbg
        elseif ((j < 0) .and. (k > 0)) then
          if (pn(j) == -pn(k)) msq(j,k)=facqq*uub_ddbg
        elseif ((j > 0) .and. (k > 0)) then
          if (pn(j)+pn(k) == +3) msq(j,k)=facqq*ud_dug
        elseif ((j < 0) .and. (k < 0)) then
          if (pn(j)+pn(k) == -3) msq(j,k)=facqq*ud_dug
        elseif ((j .ne. 0) .and. (k == 0)) then
          msq(j,k)=facqg*2._dp*ug_dudb
        elseif ((j == 0) .and. (k .ne. 0)) then
          msq(j,k)=facqg*2._dp*gu_uddb
        endif
      enddo
      enddo

      return
      end


