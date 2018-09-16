      subroutine qqb_ZH1jet_v(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Authors: R.K. Ellis and John Campbell                            *
*     May, 2001.                                                       *
*     Matrix element for Z + jet production                            *
*     in order alpha_s^2                                               *
*     averaged over initial colours and spins                          *
*     q(-p1)+qbar(-p2)-->Z^+(l(p3)+a(p4))+H(p5,p6)+g(p7)               *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'hdecaymode.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'nflav.f'
      include 'cutoff.f'
      include 'hbbparams.f'
      include 'noglue.f'
      include 'qlfirst.f'
      integer:: j,k,ig
      real(dp):: msq(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),
     & p(mxpart,4),fac,s34,s56,s127,virt5_VH,subuv,
     & msqhtautau,msqhbb,msqhgamgam,hdecay
      real(dp):: virt5_ZHtop
      real(dp):: qqbZgLL,qqbZgRR,qqbZgLR,qqbZgRL
      real(dp):: gqZqLL,gqZqRR,gqZqLR,gqZqRL
      real(dp):: qgZqLL,qgZqRR,qgZqLR,qgZqRL
      real(dp):: qbqZgLL,qbqZgRR,qbqZgLR,qbqZgRL
      real(dp):: gqbZqbLL,gqbZqbRR,gqbZqbLR,gqbZqbRL
      real(dp):: qbgZqbLL,qbgZqbRR,qbgZqbLR,qbgZqbRL
      complex(dp):: prop
      real(dp):: fac_top,cutoff_orig,sH
      integer,parameter::
     & iqqbgLL(5)=(/1,2,3,4,7/),iqqbgRR(5)=(/2,1,4,3,7/),
     & iqqbgRL(5)=(/2,1,3,4,7/),iqqbgLR(5)=(/1,2,4,3,7/),
     & iqgqLL(5)=(/1,7,3,4,2/),iqgqRR(5)=(/7,1,4,3,2/),
     & iqgqRL(5)=(/7,1,3,4,2/),iqgqLR(5)=(/1,7,4,3,2/),
     & igqqLL(5)=(/2,7,3,4,1/),igqqRR(5)=(/7,2,4,3,1/),
     & igqqRL(5)=(/7,2,3,4,1/),igqqLR(5)=(/2,7,4,3,1/)
      integer,parameter::
     & iqqbgLL9(5)=(/1,2,3,4,9/),iqqbgRR9(5)=(/2,1,4,3,9/),
     & iqqbgRL9(5)=(/2,1,3,4,9/),iqqbgLR9(5)=(/1,2,4,3,9/),
     & iqgqLL9(5)=(/1,9,3,4,2/),iqgqRR9(5)=(/9,1,4,3,2/),
     & iqgqRL9(5)=(/9,1,3,4,2/),iqgqLR9(5)=(/1,9,4,3,2/),
     & igqqLL9(5)=(/2,9,3,4,1/),igqqRR9(5)=(/9,2,4,3,1/),
     & igqqRL9(5)=(/9,2,3,4,1/),igqqLR9(5)=(/2,9,4,3,1/)
      include 'cplx.h'

      scheme='dred'
c--set msq=0 to initialize
      msq(:,:)=zip

      cutoff_orig=cutoff

c--- initialize QCDLoop
      if (qlfirst) then
        qlfirst=.false.
        call qlinit
      endif

      if(hdecaymode.ne.'wpwm') then
         ig=7
      else
         ig=9
      endif

c--calculate spinor and dot-products (using BDK type notation)
      call spinoru(ig,p,za,zb)

c--- calculate lowest order
      msq0=zip
      if (toponly .eqv. .false.) then
        call qqb_ZH1jet(p,msq0)
      endif

c----UV counterterm contains the finite renormalization to arrive
c----at MS bar scheme.
      subuv=ason2pi*xn*(epinv*(11._dp-two*real(nflav,dp)/xn)-one)/6._dp

c--   calculate propagator
      s34=s(3,4)
      s127=s(1,2)+s(1,ig)+s(2,ig)
      prop=cone/cplx2((s34-zmass**2),zmass*zwidth)
      prop=prop*s127/cplx2((s127-zmass**2),zmass*zwidth)

C     Deal with Higgs decay
      if (hdecaymode == 'tlta') then
         sH=s(5,6)+two*mtau**2
         hdecay=msqhtautau(sH)
      elseif (hdecaymode == 'bqba') then
         sH=s(5,6)+two*mb**2
         hdecay=msqhbb(sH)
c---  adjust for fixed H->bb BR if necessary
         if (FixBrHbb) then
            hdecay=hdecay*GamHbb/GamHbb0
         endif
      elseif (hdecaymode == 'gaga') then
         sH=s(5,6)
         hdecay=msqhgamgam(sH)
      elseif (hdecaymode == 'wpwm') then
         sH=s(5,6)+s(5,7)+s(5,8)+s(6,7)+s(6,8)+s(7,8)
         call hwwdecay(p,5,6,7,8,hdecay)
      else
         write(6,*) 'Unimplemented process in qqb_higgs'
         stop
      endif
      hdecay=hdecay/((sH-hmass**2)**2+(hmass*hwidth)**2)

      fac=8._dp*cf*xnsq*esq**2*gsq*hdecay
      fac=fac*gwsq*wmass**2/(one-xw)**2
      fac_top=two*cf*xn*hdecay/abs(prop**2)

!---- smalls cut for top loop piece stability (quick but dirty)
      cutoff=1.E-3_dp
      call smalls(s,ig-2,*101)
      goto 102
 101  fac_top=zip
 102  continue
      cutoff=cutoff_orig

      if(ig==7) then
      qqbZgLL=aveqq*fac_top*virt5_ZHtop(iqqbgLL,za,zb)
      qqbZgLR=aveqq*fac_top*virt5_ZHtop(iqqbgLR,za,zb)
      qqbZgRL=aveqq*fac_top*virt5_ZHtop(iqqbgRL,za,zb)
      qqbZgRR=aveqq*fac_top*virt5_ZHtop(iqqbgRR,za,zb)

      if (toponly .eqv. .false.) then
        qqbZgLL=qqbZgLL+aveqq*fac*virt5_VH(iqqbgLL,za,zb)
        qqbZgLR=qqbZgLR+aveqq*fac*virt5_VH(iqqbgLR,za,zb)
        qqbZgRL=qqbZgRL+aveqq*fac*virt5_VH(iqqbgRL,za,zb)
        qqbZgRR=qqbZgRR+aveqq*fac*virt5_VH(iqqbgRR,za,zb)
      endif

      qbqZgLL=qqbZgRL
      qbqZgLR=qqbZgRR
      qbqZgRL=qqbZgLL
      qbqZgRR=qqbZgLR

      gqZqLL=aveqg*fac_top*virt5_ZHtop(igqqLL,za,zb)
      gqZqLR=aveqg*fac_top*virt5_ZHtop(igqqLR,za,zb)
      gqZqRL=aveqg*fac_top*virt5_ZHtop(igqqRL,za,zb)
      gqZqRR=aveqg*fac_top*virt5_ZHtop(igqqRR,za,zb)

      if (toponly .eqv. .false.) then
        gqZqLL=gqZqLL+aveqg*fac*virt5_VH(igqqLL,za,zb)
        gqZqLR=gqZqLR+aveqg*fac*virt5_VH(igqqLR,za,zb)
        gqZqRL=gqZqRL+aveqg*fac*virt5_VH(igqqRL,za,zb)
        gqZqRR=gqZqRR+aveqg*fac*virt5_VH(igqqRR,za,zb)
      endif

      gqbZqbRL=gqZqLL
      gqbZqbRR=gqZqLR
      gqbZqbLL=gqZqRL
      gqbZqbLR=gqZqRR

      qgZqLL=aveqg*fac_top*virt5_ZHtop(iqgqLL,za,zb)
      qgZqLR=aveqg*fac_top*virt5_ZHtop(iqgqLR,za,zb)
      qgZqRL=aveqg*fac_top*virt5_ZHtop(iqgqRL,za,zb)
      qgZqRR=aveqg*fac_top*virt5_ZHtop(iqgqRR,za,zb)

      if (toponly .eqv. .false.) then
        qgZqLL=qgZqLL+aveqg*fac*virt5_VH(iqgqLL,za,zb)
        qgZqLR=qgZqLR+aveqg*fac*virt5_VH(iqgqLR,za,zb)
        qgZqRL=qgZqRL+aveqg*fac*virt5_VH(iqgqRL,za,zb)
        qgZqRR=qgZqRR+aveqg*fac*virt5_VH(iqgqRR,za,zb)
      endif

      qbgZqbRL=qgZqLL
      qbgZqbRR=qgZqLR
      qbgZqbLL=qgZqRL
      qbgZqbLR=qgZqRR

      elseif(ig==9) then
         qqbZgLL=aveqq*fac_top*virt5_ZHtop(iqqbgLL9,za,zb)
         qqbZgLR=aveqq*fac_top*virt5_ZHtop(iqqbgLR9,za,zb)
         qqbZgRL=aveqq*fac_top*virt5_ZHtop(iqqbgRL9,za,zb)
         qqbZgRR=aveqq*fac_top*virt5_ZHtop(iqqbgRR9,za,zb)

         if (toponly .eqv. .false.) then
            qqbZgLL=qqbZgLL+aveqq*fac*virt5_VH(iqqbgLL9,za,zb)
            qqbZgLR=qqbZgLR+aveqq*fac*virt5_VH(iqqbgLR9,za,zb)
            qqbZgRL=qqbZgRL+aveqq*fac*virt5_VH(iqqbgRL9,za,zb)
            qqbZgRR=qqbZgRR+aveqq*fac*virt5_VH(iqqbgRR9,za,zb)
         endif

         qbqZgLL=qqbZgRL
         qbqZgLR=qqbZgRR
         qbqZgRL=qqbZgLL
         qbqZgRR=qqbZgLR

         gqZqLL=aveqg*fac_top*virt5_ZHtop(igqqLL9,za,zb)
         gqZqLR=aveqg*fac_top*virt5_ZHtop(igqqLR9,za,zb)
         gqZqRL=aveqg*fac_top*virt5_ZHtop(igqqRL9,za,zb)
         gqZqRR=aveqg*fac_top*virt5_ZHtop(igqqRR9,za,zb)

         if (toponly .eqv. .false.) then
            gqZqLL=gqZqLL+aveqg*fac*virt5_VH(igqqLL9,za,zb)
            gqZqLR=gqZqLR+aveqg*fac*virt5_VH(igqqLR9,za,zb)
            gqZqRL=gqZqRL+aveqg*fac*virt5_VH(igqqRL9,za,zb)
            gqZqRR=gqZqRR+aveqg*fac*virt5_VH(igqqRR9,za,zb)
         endif

         gqbZqbRL=gqZqLL
         gqbZqbRR=gqZqLR
         gqbZqbLL=gqZqRL
         gqbZqbLR=gqZqRR

         qgZqLL=aveqg*fac_top*virt5_ZHtop(iqgqLL9,za,zb)
         qgZqLR=aveqg*fac_top*virt5_ZHtop(iqgqLR9,za,zb)
         qgZqRL=aveqg*fac_top*virt5_ZHtop(iqgqRL9,za,zb)
         qgZqRR=aveqg*fac_top*virt5_ZHtop(iqgqRR9,za,zb)

         if (toponly .eqv. .false.) then
            qgZqLL=qgZqLL+aveqg*fac*virt5_VH(iqgqLL9,za,zb)
            qgZqLR=qgZqLR+aveqg*fac*virt5_VH(iqgqLR9,za,zb)
            qgZqRL=qgZqRL+aveqg*fac*virt5_VH(iqgqRL9,za,zb)
            qgZqRR=qgZqRR+aveqg*fac*virt5_VH(iqgqRR9,za,zb)
         endif

         qbgZqbRL=qgZqLL
         qbgZqbRR=qgZqLR
         qbgZqbLL=qgZqRL
         qbgZqbLR=qgZqRR
      endif
      do j=-nflav,nflav
      do k=-nflav,nflav
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j == 0) .and. (k == 0)) then
          msq(j,k)=0._dp
      elseif ((j > 0) .and. (k < 0)) then
          msq(j,k)=+abs(Q(j)*q1+L(j)*l1*prop)**2*qqbZgLL
     &             +abs(Q(j)*q1+R(j)*r1*prop)**2*qqbZgRR
     &             +abs(Q(j)*q1+L(j)*r1*prop)**2*qqbZgLR
     &             +abs(Q(j)*q1+R(j)*l1*prop)**2*qqbZgRL
     &             -subuv*msq0(j,k)
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=+abs(Q(k)*q1+L(k)*l1*prop)**2*qbqZgLL
     &             +abs(Q(k)*q1+R(k)*r1*prop)**2*qbqZgRR
     &             +abs(Q(k)*q1+L(k)*r1*prop)**2*qbqZgLR
     &             +abs(Q(k)*q1+R(k)*l1*prop)**2*qbqZgRL
     &             -subuv*msq0(j,k)
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=+abs(Q(j)*q1+L(j)*l1*prop)**2*qgZqLL
     &             +abs(Q(j)*q1+R(j)*r1*prop)**2*qgZqRR
     &             +abs(Q(j)*q1+L(j)*r1*prop)**2*qgZqLR
     &             +abs(Q(j)*q1+R(j)*l1*prop)**2*qgZqRL
     &             -subuv*msq0(j,k)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=+abs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqbLL
     &             +abs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqbRR
     &             +abs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqbLR
     &             +abs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqbRL
     &             -subuv*msq0(j,k)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=+abs(Q(k)*q1+L(k)*l1*prop)**2*gqZqLL
     &             +abs(Q(k)*q1+R(k)*r1*prop)**2*gqZqRR
     &             +abs(Q(k)*q1+L(k)*r1*prop)**2*gqZqLR
     &             +abs(Q(k)*q1+R(k)*l1*prop)**2*gqZqRL
     &             -subuv*msq0(j,k)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=+abs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqbLL
     &             +abs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqbRR
     &             +abs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqbLR
     &             +abs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqbRL
     &             -subuv*msq0(j,k)
      endif

   19 continue
      enddo
      enddo

      return
      end
