      subroutine ZZ_Hqq_g(p,msq)
      implicit none
      include 'types.f'

c--- Weak Boson Fusion by Z-Z exchange only
c---Matrix element squared averaged over initial colors and spins
c
c     q(-p1)+q(-p2) -->  H(p3,p4)+q(p5)+q(p6)+g(p7)
c                           |
c                           |
c                           |
c                           ---> b(p3)+bbar(p4)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'hdecaymode.f'
      integer:: j,k,nup,ndo
      real(dp):: p(mxpart,4),facqq,facqg,s34
      real(dp):: msq(-nf:nf,-nf:nf),hdecay,
     & ud_udg_LL,udb_udbg_LL,ud_udg_LR,udb_udbg_LR,
     & ug_uddb_LL,ug_uddb_LR,gu_ddbu_LL,gu_ddbu_LR,
     & msqgamgam

      parameter (nup=2,ndo=3)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call dotem(7,p,s)

      s34=(p(3,4)+p(4,4))**2
     & -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(hmass)
      else
      write(6,*) 'Unimplemented process in ZZ_hqq_g'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
      facqq=0.25_dp*gsq*Cf*gwsq**3*hdecay
      facqg=-facqq*(aveqg/aveqq)
C Extra factor of gsq*Cf compared to lowest order

C q-q and qbar-qbar
c--- u(1)+d(2) -> u(5)+d(6)+g(p7)
c--- ub(1)+db(2) -> ub(5)+db(6)+g(p7)
      call msq_gpieces_zz(1,2,5,6,7,ud_udg_LL,ud_udg_LR)

C q-qbar and qbar-q
c--- u(1)+db(2) -> u(5)+db(6)+g(p7)
c--- ub(1)+d(2) -> ub(5)+d(6)+g(p7)
      call msq_gpieces_zz(1,6,5,2,7,udb_udbg_LL,udb_udbg_LR)

C q-g and qbar-g
c--- u(1)+g(2) -> u(5)+[q](6)+[qb](7)
c--- ub(1)+g(2) -> ub(5)+[qb](6)+[q](7)
      call msq_gpieces_zz(1,7,5,6,2,ug_uddb_LL,ug_uddb_LR)

C g-q and g-qbar
c--- g(1)+u(2) -> [q](5)+[qb](6)+u(7)
c--- g(1)+ub(2) -> [qb](5)+[q](6)+ub(7)
      call msq_gpieces_zz(6,2,5,7,1,gu_ddbu_LL,gu_ddbu_LR)

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


      subroutine msq_gpieces_zz(i1,i2,i5,i6,i7,zll,zlr)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      real(dp):: zll,zlr,htheta
      real(dp):: msqll_a,msqll_b,msqlr_a,msqlr_b,
     & msq_gsamehel,msq_gopphel
      real(dp):: propz,x
      real(dp):: s15,s26,s157,s267
      integer:: i1,i2,i5,i6,i7
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(x)=half+sign(half,x)
      propz(s15)=sign(one,(s15-zmass**2))
     & *sqrt(sqrt(1._dp-xw)/xw/2._dp/zmass
     & *((s15-zmass**2)**2+htheta(s15)*(zmass*zwidth)**2))

c--- Calculate some invariants
      s15=s(i1,i5)
      s26=s(i2,i6)
      s157=s15+s(i1,i7)+s(i5,i7)
      s267=s26+s(i2,i7)+s(i6,i7)

      msqll_a=msq_gsamehel(i1,i2,i5,i6,i7)
      msqll_b=msq_gsamehel(i2,i1,i6,i5,i7)
      msqlr_a=msq_gopphel(i1,i2,i5,i6,i7)
      msqlr_b=msq_gopphel(i2,i1,i6,i5,i7)

c--- catch the unwanted diagrams for the gluon-quark processes
      if (i7 == 1) then
        msqll_b=0._dp
        msqlr_b=0._dp
      endif
      if (i7 == 2) then
        msqll_a=0._dp
        msqlr_a=0._dp
      endif

      zll=msqll_a/(propz(s157)*propz(s26))**2
     &   +msqll_b/(propz(s267)*propz(s15))**2

      zlr=msqlr_a/(propz(s157)*propz(s26))**2
     &   +msqlr_b/(propz(s267)*propz(s15))**2

      return
      end
