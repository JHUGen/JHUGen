      subroutine qqb_wh(p,msq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c for nwz=1
c     q(-p1)+qbar(-p2) -->  H  + W
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+b(p6)
c for nwz=-1
c     q(-p1)+qbar(-p2) -->  H  + W
c                           |    |
c                           |    --> e^-(p3)+nubar(p4)
c                           |
c                           ---> b(p5)+b(p6)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'zprods_decl.f'
      include 'hbbparams.f'
      include 'hdecaymode.f'
      integer:: j,k
      real(dp):: p(mxpart,4)
      real(dp):: s,prop,fac,qqbWH,qbqWH,s56
      real(dp):: msq(-nf:nf,-nf:nf),hdecay,msqhtautau,msqhbb
      complex(dp):: test,qqb_WH_HtopEFT

      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

      msq(:,:)=zip

c---calculate the 2 W propagators
      prop=     ((s(1,2)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

      fac=xn*gwsq**3*wmass**2/prop
c--- Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          s56=s(5,6)+2._dp*mtau**2
          hdecay=msqhtautau(s56)
      elseif (hdecaymode == 'bqba') then
          s56=s(5,6)+2._dp*mb**2
          hdecay=msqhbb(s56)
      else
        write(6,*) 'Unimplemented process in qqb_wh'
        stop
      endif
      hdecay=hdecay/((s56-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay

c--- adjust for fixed H->bb BR if necessary
      if ((FixBrHbb) .and. (hdecaymode == 'bqba')) then
        fac=fac*GamHbb/GamHbb0
      endif

c-- Old form of this matrix element (modified to facilitate extension
c--- to H->WW decay)
c      spinave only (color factors cancel)
c       fac=spinave*gw**8*0.5_dp*mbsq*(s56-4._dp*mb**2)/prop

      qqbWH=aveqq*fac*s(1,4)*s(2,3)
      qbqWH=aveqq*fac*s(2,4)*s(1,3)

!==== testing code (debug)
!      qqbWH=aveqq*qqb_wh_HtopEFT(1,2,3,4,p)*hdecay
!      qbqWH=aveqq*qqb_wh_HtopEFT(2,1,3,4,p)*hdecay
!      if(FixBrHbb) then
!         qqbWH=qqbWH*GamHbb/GamHbb0
!         qbqWH=qbqWH*GamHbb/GamHbb0
!      endif

      do j=-nf,nf
      do k=-nf,nf
        if ((j > 0) .and. (k < 0)) msq(j,k)=Vsq(j,k)*qqbWH
        if ((j < 0) .and. (k > 0)) msq(j,k)=Vsq(j,k)*qbqWH
      enddo
      enddo

c      if(1==2) then
c         p(1,4)=  -6.0000000000000000_dp
c         p(1,1)=   3.9291644036717099_dp
c         p(1,2)=   2.2672487230627749_dp
c         p(1,3)=   3.9269899817403862_dp
c         p(2,4)=   4.0000000000000000_dp
c         p(2,1)=   2.6000000000000001_dp
c         p(2,2)=   0.0000000000000000_dp
c         p(2,3)=   3.0397368307141326_dp
c         p(3,4)=   1.7142857142857142_dp
c         p(3,1)= -0.63157894736842102_dp
c         p(3,2)=   1.5937012089614160_dp
c         p(3,3)=   0.0000000000000000_dp
c         p(4,4)=   2.0000000000000000_dp
c         p(4,1)= -0.42857142857142866_dp
c         p(4,2)=  0.90350790290525151_dp
c         p(4,3)=   1.7320508075688772_dp
c         p(5,4)=  -1.7142857142857142_dp
c         p(5,1)=  -5.4690140277318600_dp
c         p(5,2)=  -4.7644578349294422_dp
c         p(5,3)=  -8.6987776200233959_dp
c      endif
c      call spinorz(4,p,za,zb)
c      test=qqb_WH_HtopEFT_loop(1,2,3,4,za,zb)


      return
      end

