      subroutine qqb_WH1jet_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
c----Matrix element for WH production
C----averaged over initial colours and spins
c    contracted with the vector n(mu)
C For nwz=+1
c     u(-p1)+dbar(-p2)--> g(p7)+ W^+(n(p3)+e^+(p4))+H(p5,p6)
C For nwz=-1
c     d(-p1)+ubar(-p2)--> g(p7)+ W^-(e^-(p3)+nbar(p4))+H(p5,p6)
c---
c---ip emitter
c---kp spectator
c---in label of gluon which is contracted with n
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'hdecaymode.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'ckm.f'
      include 'masses.f'
      include 'hbbparams.f'
      integer:: j,k,in,ig
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),w1jetn,
     & p1p2(-1:1,-1:1),n(4),s127,hdecay,msqhgamgam,msqhtautau,
     & msqhbb,fac,prop,sH

      msq(:,:)=zip
      p1p2(:,:)=zip

      if (hdecaymode == 'wpwm') then
        ig=9
      else
        ig=7
      endif
      call dotem(ig,p,s)

      s127=s(1,2)+s(1,ig)+s(2,ig)
      prop=     ((s127-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          sH=s(5,6)+2._dp*mtau**2
          hdecay=msqhtautau(sH)
      elseif (hdecaymode == 'bqba') then
          sH=s(5,6)+2._dp*mb**2
          hdecay=msqhbb(sH)
c--- adjust for fixed H->bb BR if necessary
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
          write(6,*) 'Unimplemented decay mode in qqb_WH1jet_gvec'
          stop
      endif
      hdecay=hdecay/((sH-hmass**2)**2+(hmass*hwidth)**2)

c---calculate the propagator
      fac=two*gsq*V*gwsq**3*wmass**2/prop*hdecay

      if (in == 1) then
      p1p2(0,-1)=-aveqg*fac*w1jetn(ig,2,3,4,1,p,n)
      p1p2(0,+1)=-aveqg*fac*w1jetn(2,ig,3,4,1,p,n)
      elseif (in == 2) then
      p1p2(+1,0)=-aveqg*fac*w1jetn(1,ig,3,4,2,p,n)
      p1p2(-1,0)=-aveqg*fac*w1jetn(ig,1,3,4,2,p,n)
      elseif (in == ig) then
      p1p2(1,-1)=+aveqq*fac*w1jetn(1,2,3,4,ig,p,n)
      p1p2(-1,1)=+aveqq*fac*w1jetn(2,1,3,4,ig,p,n)
      endif

      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=Vsq(j,k)*p1p2(1,-1)
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=Vsq(j,k)*p1p2(-1,1)
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=
     &    (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*p1p2(+1,0)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=
     &    (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*p1p2(-1,0)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=
     &    (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*p1p2(0,+1)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=
     &    (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*p1p2(0,-1)
      endif

      enddo
      enddo

      return
      end



