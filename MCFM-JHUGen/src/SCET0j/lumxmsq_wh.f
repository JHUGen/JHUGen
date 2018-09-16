      subroutine lumxmsq_wh(p,xx,z1,z2,QB,order,xmsq)
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
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'ckm.f'
      include 'hdecaymode.f'
      include 'scet_const.f'
      include 'hbbparams.f'
      include 'noglue.f'
      integer:: j,k,ih1,ih2,m,n,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),s,fac,qqbWH,qbqWH,hdecay,prop,sH,
     & msqhtautau,msqhbb,msqhgamgam,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & z1,z2,QB(2),lum0,lum1(-1:1),lum2(-1:3),bit,
     & msq(-nf:nf,-nf:nf),assemble
      common/density/ih1,ih2
      real(dp):: qqb_wh_HtopEFT,qqbWHtoploop,qbqWHtoploop

      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &           -p(j,2)*p(k,2)-p(j,3)*p(k,3))

c---calculate the two W propagators
      prop=     ((s(1,2)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

      fac=xn*gwsq**3*wmass**2/prop
C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          sH=s(5,6)+2._dp*mtau**2
          hdecay=msqhtautau(sH)
      elseif (hdecaymode == 'bqba') then
          sH=s(5,6)+2._dp*mb**2
          hdecay=msqhbb(sH)
c--- adjust for fixed H->bb BR if necessary
          if (FixBrHbb) then
            fac=fac*GamHbb/GamHbb0
          endif
      elseif (hdecaymode == 'gaga') then
          sH=s(5,6)
          hdecay=msqhgamgam(sH)
      elseif (hdecaymode == 'wpwm') then
          sH=s(5,6)+s(5,7)+s(5,8)+s(6,7)+s(6,8)+s(7,8)
          call hwwdecay(p,5,6,7,8,hdecay)
      else
          write(6,*) 'Unimplemented decay mode in lumxmsq_wh.f'
          stop
      endif
      hdecay=hdecay/((sH-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay

      qqbWH=aveqq*fac*s(1,4)*s(2,3)
      qbqWH=aveqq*fac*s(2,4)*s(1,3)

      call softqqbis(order,soft1,soft2)
      call hardqq(s(1,2),musq,hard)

      if (order >= 0) then
      call fdist(ih1,xx(1),facscale,beama0)
      call fdist(ih2,xx(2),facscale,beamb0)
      endif
      if (order >= 1) then
      call xbeam1bis(ih1,z1,xx(1),QB(1),beama1)
      call xbeam1bis(ih2,z2,xx(2),QB(2),beamb1)
      endif
      if (order >= 2) then
      call xbeam2bis(ih1,z1,xx(1),QB(1),beama2)
      call xbeam2bis(ih2,z2,xx(2),QB(2),beamb2)
      qqbWHtoploop=aveqq*qqb_wh_HtopEFT(1,2,3,4,p)*hdecay
      qbqWHtoploop=aveqq*qqb_wh_HtopEFT(2,1,3,4,p)*hdecay
      if ((FixBrHbb).and.(hdecaymode=='bqba')) then
         qqbWHtoploop=qqbWHtoploop*GamHbb/GamHbb0
         qbqWHtoploop=qbqWHtoploop*GamHbb/GamHbb0
      endif
      endif


      if (toponly) then
        qqbWH=zip
        qbqWH=zip
      endif

      xmsq=zip
      do j=-nf,nf
      do k=-nf,nf
      if (j*k >= 0) cycle ! skip gluons, qq, aa

      bit=assemble(order,
     & beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     & beama2(j,:),beamb2(k,:),soft1,soft2,hard)


      if ((j > 0) .and. (k < 0)) then
        bit=bit*Vsq(j,k)*qqbWH
        if(order>=2) then
           bit=bit+Vsq(j,k)*beama0(j)*beamb0(k)*qqbWHtoploop
        endif
      elseif ((j < 0) .and. (k > 0)) then
        bit=bit*Vsq(j,k)*qbqWH
        if(order>=2) then
           bit=bit+Vsq(j,k)*beama0(j)*beamb0(k)*qbqWHtoploop
        endif
      else
        bit=zip
      endif

      xmsq=xmsq+bit

      enddo
      enddo

      return
      end
