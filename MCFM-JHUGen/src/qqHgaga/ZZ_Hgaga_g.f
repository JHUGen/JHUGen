      subroutine ZZ_Hgaga_g(p,msq)
      implicit none 
c--- Weak Boson Fusion by Z-Z exchange only
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
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      integer j,k,nup,ndo
      double precision p(mxpart,4),facqq,facqg
      double precision msq(-nf:nf,-nf:nf),hdecay,
     . ud_udg_LL,udb_udbg_LL,ud_udg_LR,udb_udbg_LR,
     . ug_uddb_LL,ug_uddb_LR,gu_ddbu_LL,gu_ddbu_LR
      double precision msqgamgam

      parameter (nup=2,ndo=3)

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
c--- u(1)+d(2) -> u(5)+d(6)+g(p7)
c--- ub(1)+db(2) -> ub(5)+db(6)+g(p7)
      call msq_gpieces_gaga_zz(1,2,5,6,7,ud_udg_LL,ud_udg_LR)

C q-qbar and qbar-q
c--- u(1)+db(2) -> u(5)+db(6)+g(p7)
c--- ub(1)+d(2) -> ub(5)+d(6)+g(p7)
      call msq_gpieces_gaga_zz(1,6,5,2,7,udb_udbg_LL,udb_udbg_LR)

C q-g and qbar-g
c--- u(1)+g(2) -> u(5)+[q](6)+[qb](7)
c--- ub(1)+g(2) -> ub(5)+[qb](6)+[q](7)
      call msq_gpieces_gaga_zz(1,7,5,6,2,ug_uddb_LL,ug_uddb_LR)

C g-q and g-qbar
c--- g(1)+u(2) -> [q](5)+[qb](6)+u(7)
c--- g(1)+ub(2) -> [qb](5)+[q](6)+ub(7)
      call msq_gpieces_gaga_zz(6,2,5,7,1,gu_ddbu_LL,gu_ddbu_LR)

      do j=-nf,nf
      do k=-nf,nf
        if     ((j .gt. 0) .and. (k .lt. 0)) then
          msq(j,k)=facqq*(
     .     +udb_udbg_LL*((L(+j)*L(-k))**2+(R(+j)*R(-k))**2)
     .     +udb_udbg_LR*((L(+j)*R(-k))**2+(R(+j)*L(-k))**2))
        elseif ((j .lt. 0) .and. (k .gt. 0)) then
          msq(j,k)=facqq*(
     .     +udb_udbg_LL*((L(-j)*L(k))**2+(R(-j)*R(k))**2)
     .     +udb_udbg_LR*((L(-j)*R(k))**2+(R(-j)*L(k))**2))
        elseif ((j .gt. 0) .and. (k .gt. 0)) then
          msq(j,k)=facqq*(
     .     +ud_udg_LL*((L(+j)*L(+k))**2+(R(+j)*R(+k))**2)
     .     +ud_udg_LR*((L(+j)*R(+k))**2+(R(+j)*L(+k))**2))
        elseif ((j .lt. 0) .and. (k .lt. 0)) then
          msq(j,k)=facqq*(
     .     +ud_udg_LL*((L(-j)*L(-k))**2+(R(-j)*R(-k))**2)
     .     +ud_udg_LR*((L(-j)*R(-k))**2+(R(-j)*L(-k))**2))
        elseif ((j .ne. 0) .and. (k .eq. 0)) then
          msq(j,k)=facqg*dfloat(nup)*(
     .     +ug_uddb_LL*((L(abs(j))*L(+2))**2+(R(abs(j))*R(+2))**2)
     .     +ug_uddb_LR*((L(abs(j))*R(+2))**2+(R(abs(j))*L(+2))**2))
     .            +facqg*dfloat(ndo)*(
     .     +ug_uddb_LL*((L(abs(j))*L(+1))**2+(R(abs(j))*R(+1))**2)
     .     +ug_uddb_LR*((L(abs(j))*R(+1))**2+(R(abs(j))*L(+1))**2))
        elseif ((j .eq. 0) .and. (k .ne. 0)) then
          msq(j,k)=facqg*dfloat(nup)*(
     .     +gu_ddbu_LL*((L(abs(k))*L(+2))**2+(R(abs(k))*R(+2))**2)
     .     +gu_ddbu_LR*((L(abs(k))*R(+2))**2+(R(abs(k))*L(+2))**2))
     .            +facqg*dfloat(ndo)*(
     .     +gu_ddbu_LL*((L(abs(k))*L(+1))**2+(R(abs(k))*R(+1))**2)
     .     +gu_ddbu_LR*((L(abs(k))*R(+1))**2+(R(abs(k))*L(+1))**2))
        endif
      enddo
      enddo

      return
      end


      subroutine msq_gpieces_gaga_zz(i1,i2,i5,i6,i7,zll,zlr)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      double precision zll,zlr,htheta
      double precision msqll_a,msqll_b,msqlr_a,msqlr_b,
     . msq_gsamehel,msq_gopphel
      double precision propz,x
      double precision s15,s26,s157,s267
      integer i1,i2,i5,i6,i7
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(x)=half+sign(half,x)
      propz(s15)=sign(one,(s15-zmass**2))
     . *dsqrt(dsqrt(1d0-xw)/xw/2d0/zmass
     . *((s15-zmass**2)**2+htheta(s15)*(zmass*zwidth)**2))

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
      if (i7 .eq. 1) then
        msqll_b=0d0
        msqlr_b=0d0
      endif
      if (i7 .eq. 2) then
        msqll_a=0d0
        msqlr_a=0d0
      endif

      zll=msqll_a/(propz(s157)*propz(s26))**2
     .   +msqll_b/(propz(s267)*propz(s15))**2

      zlr=msqlr_a/(propz(s157)*propz(s26))**2
     .   +msqlr_b/(propz(s267)*propz(s15))**2

      return
      end
