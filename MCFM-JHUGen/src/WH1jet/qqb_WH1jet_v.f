      subroutine qqb_wh1jet_v(p,msq)
      implicit none
      include 'types.f'
!----Matrix element for WH + jet production
!----in order alpha_s^2
!---Matrix element squared averaged over initial colors and spins
!---for nwz=1
!     q(-p1)+qbar(-p2) -->  H  + W +g(p7)
!                           |    |
!                           |    --> nu(p3)+e^+(p4)
!                           |
!                           ---> b(p5)+b(p6)
!---for nwz=-1
!     q(-p1)+qbar(-p2) -->  H  + W +g(p7)
!                           |    |
!                           |    --> e^-(p3)+nubar(p4)
!                           |
!                           ---> b(p5)+b(p6)
!   for the moment --- radiation only from initial line
!---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'hdecaymode.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'nflav.f'
      include 'hbbparams.f'
      include 'cutoff.f'
      include 'noglue.f'
      include 'qlfirst.f'
      integer:: j,k,ig
      real(dp):: msq(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),
     & p(mxpart,4),fac,prop,subuv,hdecay,
     & qqbWg,qbqWg,qgWq,gqWq,qbgWqb,gqbWqb,s127,s34,sH,
     & msqhtautau,msqhbb,msqhgamgam
      integer,parameter::
     & iqqbg(5)=(/1,2,3,4,7/),iqgq(5)=(/1,7,3,4,2/),
     & igqq(5)=(/2,7,3,4,1/),iqbqg(5)=(/2,1,3,4,7/),
     & iqbgqb(5)=(/7,1,3,4,2/),igqbqb(5)=(/7,2,3,4,1/)
      integer,parameter::
     & iqqbg9(5)=(/1,2,3,4,9/),iqgq9(5)=(/1,9,3,4,2/),
     & igqq9(5)=(/2,9,3,4,1/),iqbqg9(5)=(/2,1,3,4,9/),
     & iqbgqb9(5)=(/9,1,3,4,2/),igqbqb9(5)=(/9,2,3,4,1/)
      real(dp):: fac_top,virt5_VHtop,virt5_VH
      real(dp):: cutoff_orig

      cutoff_orig=cutoff

c--set msq=0 to initialize
      msq(:,:)=0._dp
      scheme='dred'

c--- initialize QCDLoop
      if (qlfirst) then
        qlfirst=.false.
        call qlinit
      endif

c-- if Gflag=.false. then only the endpoint contributions from the
c-- 4-quark diagrams are included, ie. no pole subtraction for this
c-- piece. Therefore return 0.
c      if (Gflag .eqv. .false.) return

c---  calculate lowest order
      msq0=zip
      if(toponly.eqv. .false.) then
         call qqb_WH1jet(p,msq0)
      endif

c--- UV counterterm contains the finite renormalization to arrive
c--- at MS bar scheme.
      subuv=ason2pi*xn*(epinv*(11._dp-two*real(nflav,dp)/xn)-1._dp)/six

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          sH=s(5,6)+2._dp*mtau**2
          hdecay=msqhtautau(sH)
          ig=7
      elseif (hdecaymode == 'bqba') then
          sH=s(5,6)+2._dp*mb**2
          hdecay=msqhbb(sH)
          ig=7
c--- adjust for fixed H->bb BR if necessary
          if (FixBrHbb) then
            hdecay=hdecay*GamHbb/GamHbb0
          endif
      elseif (hdecaymode == 'gaga') then
          sH=s(5,6)
          hdecay=msqhgamgam(sH)
          ig=7
      elseif (hdecaymode == 'wpwm') then
          sH=s(5,6)+s(5,7)+s(5,8)+s(6,7)+s(6,8)+s(7,8)
          call hwwdecay(p,5,6,7,8,hdecay)
          ig=9
      else
          write(6,*) 'Unimplemented decay mode in qqb_WH1jet_v'
      stop
      endif

c--- calculate propagators
      s127=s(1,2)+s(1,ig)+s(2,ig)
      s34=s(3,4)
      prop=s127**2/((s127-wmass**2)**2+(wmass*wwidth)**2)
     & /((s34-wmass**2)**2+(wmass*wwidth)**2)

      hdecay=hdecay/((sH-hmass**2)**2+(hmass*hwidth)**2)

      fac=2._dp*cf*xnsq*gsq*gwsq**3*wmass**2*prop*hdecay
!===== factors for top mediated diagram (overall color and Higgs decay)
      fac_top=2._dp*cf*xn*hdecay
!---- smalls cut for top loop piece stability (quick but inefficient)
      cutoff=1.E-3_dp
      call smalls(s,ig-2,*101)
      goto 102
 101  fac_top=zip
 102  continue
      cutoff=cutoff_orig

c--calculate spinor and dot-products (using BDK type notation)
      call spinoru(ig,p,za,zb)

      if (ig == 7) then
        qqbWg=aveqq*(
     &       +fac_top*virt5_VHtop(iqqbg,za,zb))
        qbqWg=aveqq*(
     &       +fac_top*virt5_VHtop(iqbqg,za,zb))
        gqWq=aveqg*(
     &       +fac_top*virt5_VHtop(igqq,za,zb))
        qgWq=aveqg*(
     &       +fac_top*virt5_VHtop(iqgq,za,zb))
        gqbWqb=aveqg*(
     &       +fac_top*virt5_VHtop(igqbqb,za,zb))
        qbgWqb=aveqg*(
     &       +fac_top*virt5_VHtop(iqbgqb,za,zb))

        if(toponly.eqv..false.) then
           qqbWg=qqbWg  +aveqq*fac*virt5_VH(iqqbg,za,zb)
           qbqWg=qbqWg  +aveqq*fac*virt5_VH(iqbqg,za,zb)
           gqWq= gqWq   +aveqg*fac*virt5_VH(igqq,za,zb)
           qgWq= qgWq   +aveqg*fac*virt5_VH(iqgq,za,zb)
           gqbWqb=gqbWqb+aveqg*fac*virt5_VH(igqbqb,za,zb)
           qbgWqb=qbgWqb+aveqg*fac*virt5_VH(iqbgqb,za,zb)
        endif

      else
        qqbWg=aveqq*(
     &       +fac_top*virt5_VHtop(iqqbg9,za,zb))
        qbqWg=aveqq*(
     &       +fac_top*virt5_VHtop(iqbqg9,za,zb))
        gqWq=aveqg*(
     &       +fac_top*virt5_VHtop(igqq9,za,zb))
        qgWq=aveqg*(
     &       +fac_top*virt5_VHtop(iqgq9,za,zb))
        gqbWqb=aveqg*(
     &       +fac_top*virt5_VHtop(igqbqb9,za,zb))
        qbgWqb=aveqg*(
     &       +fac_top*virt5_VHtop(iqbgqb9,za,zb))

        if(toponly.eqv..false.) then
           qqbWg=qqbWg+aveqq*fac*virt5_VH(iqqbg9,za,zb)
           qbqWg=qbqWg+aveqq*fac*virt5_VH(iqbqg9,za,zb)
           gqWq= gqWq +aveqg*fac*virt5_VH(igqq9,za,zb)
           qgWq= qgWq +aveqg*fac*virt5_VH(iqgq9,za,zb)
           gqbWqb=gqbWqb+aveqg*fac*virt5_VH(igqbqb9,za,zb)
           qbgWqb=qbgWqb+aveqg*fac*virt5_VH(iqbgqb9,za,zb)
        endif

      endif

      do j=-nf,nf
      do k=-nf,nf

      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=Vsq(j,k)*qqbWg-subuv*msq0(j,k)
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=Vsq(j,k)*qbqWg-subuv*msq0(j,k)
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=
     &   (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWq
     &     -subuv*msq0(j,k)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=
     &    (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWqb
     &     -subuv*msq0(j,k)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=
     &    (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWq
     &     -subuv*msq0(j,k)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=
     &    (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWqb
     &     -subuv*msq0(j,k)
      endif

      enddo
      enddo

      return
      end
