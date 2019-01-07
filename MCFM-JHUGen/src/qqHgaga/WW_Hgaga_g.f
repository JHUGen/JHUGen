      subroutine WW_Hgaga_g(p,msq)
      implicit none 
c--- Weak Boson Fusion by W-W exchange only
c---Matrix element squared averaged over initial colors and spins
c
c     q(-p1)+q(-p2) -->  H(p3,p4)+q(p5)+q(p6)+g(p7) 
c                           |
c                           |
c                           |
c                           ---> b(p3)+bbar(p4)
c---- Extension to photon decay contributed by Fabian Stoeckli
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      integer j,k
      double precision p(mxpart,4),facqq,facqg
      double precision msq(-nf:nf,-nf:nf),hdecay,
     . ud_dug,uub_ddbg,ug_dudb,gu_uddb
      double precision msqgamgam

      integer,parameter::pn(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      
      call dotem(7,p,s)

      hdecay=msqgamgam(hmass)/((s(3,4)-hmass**2)**2+(hmass*hwidth)**2)

      facqq=0.25d0*gsq*Cf*gwsq**3*hdecay
      facqg=-facqq*(aveqg/aveqq)
C Extra factor of gsq*Cf compared to lowest order

C q-q and qbar-qbar
c--- u(1)+d(2) -> d(5)+u(6)+g(p7)
c--- ub(1)+db(2) -> db(5)+ub(6)+g(p7)
      call msq_gpieces_gaga_ww(1,2,6,5,7,ud_dug)

C q-qbar and qbar-q
c--- u(1)+ub(2) -> d(5)+db(6)+g(p7)
c--- ub(1)+u(2) -> db(5)+d(6)+g(p7)
      call msq_gpieces_gaga_ww(1,6,2,5,7,uub_ddbg)

C q-g and qbar-g
c--- u(1)+g(2) -> d(5)+[u/c](6)+[db/sb](7)
c--- ub(1)+g(2) -> db(5)+[ub/cb](6)+[d/s](7)
      call msq_gpieces_gaga_ww(1,7,6,5,2,ug_dudb)

C g-q and g-qbar
c--- g(1)+u(2) -> [u/c](5)+[db/sb](6)+d(7)
c--- g(1)+ub(2) -> [ub/cb](5)+[d/s](6)+db(7)
      call msq_gpieces_gaga_ww(6,2,7,5,1,gu_uddb)

c--- Only loop up to (nf-1) to avoid b->t transitions
      do j=-(nf-1),nf-1
      do k=-(nf-1),nf-1
        if     ((j .gt. 0) .and. (k .lt. 0)) then
          if (pn(j) .eq. -pn(k)) msq(j,k)=facqq*uub_ddbg
        elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (pn(j) .eq. -pn(k)) msq(j,k)=facqq*uub_ddbg
        elseif ((j .gt. 0) .and. (k .gt. 0)) then
          if (pn(j)+pn(k) .eq. +3) msq(j,k)=facqq*ud_dug
        elseif ((j .lt. 0) .and. (k .lt. 0)) then
          if (pn(j)+pn(k) .eq. -3) msq(j,k)=facqq*ud_dug
        elseif ((j .ne. 0) .and. (k .eq. 0)) then
          msq(j,k)=facqg*2d0*ug_dudb
        elseif ((j .eq. 0) .and. (k .ne. 0)) then
          msq(j,k)=facqg*2d0*gu_uddb
        endif
      enddo
      enddo

      return
      end


      subroutine msq_gpieces_gaga_ww(i1,i2,i5,i6,i7,wll)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      double precision wll,htheta
      double precision msqll_c,msqll_d,msq_gsamehel
      double precision propw,x
      double precision s16,s25,s167,s257
      integer i1,i2,i5,i6,i7
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(x)=half+sign(half,x)
      propw(s16)=sign(one,(s16-wmass**2))*dsqrt(
     .((s16-wmass**2)**2+htheta(s16)*(wmass*wwidth)**2)/wmass)

c--- Calculate some invariants
      s16=s(i1,i6)
      s25=s(i2,i5)
      s167=s16+s(i1,i7)+s(i6,i7)
      s257=s25+s(i2,i7)+s(i5,i7)

      msqll_c=msq_gsamehel(i1,i2,i6,i5,i7)
      msqll_d=msq_gsamehel(i2,i1,i5,i6,i7)

c--- catch the unwanted diagrams for the gluon-quark processes     
      if (i7 .eq. 1) then
        msqll_d=0d0
      endif
      if (i7 .eq. 2) then
        msqll_c=0d0
      endif

      wll=msqll_c/(propw(s167)*propw(s25))**2
     .   +msqll_d/(propw(s257)*propw(s16))**2
      return
      end
