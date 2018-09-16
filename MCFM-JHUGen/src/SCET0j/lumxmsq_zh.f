      subroutine lumxmsq_zh(p,xx,z1,z2,QB,order,xmsq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  H  + Z
c                           |    |
c                           |     ->fermion(p3)+antifermion(p4)
c                           |
c                            ---> b(p5)+b(p6)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
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
      real(dp):: p(mxpart,4),s,fac,qqbWH,qbqWH,hdecay,prop,s56,
     & msqhtautau,msqhbb,msqgamgam,v2(2),q1423,q2413,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & z1,z2,QB(2),lum0,lum1(-1:1),lum2(-1:3),bit,
     &     msq(-nf:nf,-nf:nf),assemble
      real(dp) :: sH
      real(dp):: ggZH,qqbZHtop(1:5),qbqZHtop(1:5),qqb_ZH_VItop,qqb_ZH_VIItop
      common/density/ih1,ih2
      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &     -p(j,2)*p(k,2)-p(j,3)*p(k,3))


      v2(1)=l1
      v2(2)=r1

c---calculate the 2 Z propagators
      prop=     ((s(1,2)-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)

      fac=xn*4._dp*(xw/(1._dp-xw))**2*gwsq**3*wmass**2/prop
C     Deal with Higgs decay
      if (hdecaymode == 'tlta') then
         sH=s(5,6)+2._dp*mtau**2
         call htautaudecay(p,5,6,hdecay)
      elseif (hdecaymode == 'bqba') then
         sH=s(5,6)+2*mb**2
         call hbbdecay(p,5,6,hdecay)
c---  adjust for fixed H->bb BR if necessary
         if (FixBrHbb) then
            fac=fac*GamHbb/GamHbb0
         endif
      elseif (hdecaymode == 'gaga') then
         sH=s(5,6)
         hdecay=msqgamgam(hmass)
      elseif (hdecaymode == 'wpwm') then
         sH=s(5,6)+s(5,7)+s(5,8)+s(6,7)+s(6,8)+s(7,8)
         call hwwdecay(p,5,6,7,8,hdecay)
      else
         write(6,*) 'Unimplemented process in gg_hgg_v',hdecaymode
         stop
      endif
      hdecay=hdecay/((sH-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay

      q1423=aveqq*fac*s(1,4)*s(2,3)
      q2413=aveqq*fac*s(2,4)*s(1,3)

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
      qqbZHtop(:)=zip
      qbqZHtop(:)=zip
      ggZH=zip
!====== three additional pieces, VI, VII and gg initiated loops
!====== VI and VII can be added together
      call gg_zh(p,ggZH)
      endif

      if (FixBrHbb .and. (hdecaymode=='bqba')) then
         ggZH=ggZH*GamHbb/GamHbb0
      endif

      if (toponly) then
        q1423=zip
        q2413=zip
!        ggZH=zip
      endif

      xmsq=zip
      do j=-nf,nf
      k=-j
      bit=zip

      if (j == 0.and. order < 2) cycle ! skip gluons
      if (j==0) then
         bit=zip
      else
         bit=assemble(order,
     &        beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &        beama2(j,:),beamb2(k,:),soft1,soft2,hard)
      endif

      if (j > 0) then
        bit=bit*(((L(j)*v2(1))**2+(R(j)*v2(2))**2)*q1423
     &          +((L(j)*v2(2))**2+(R(j)*v2(1))**2)*q2413)
        if(order >=2) then
           qqbZHtop(j)=aveqq*hdecay*(qqb_ZH_VItop(1,2,3,4,p,abs(j))
     &     +qqb_ZH_VIItop(1,2,3,4,p,abs(j)))
           if (FixBrHbb .and. (hdecaymode=='bqba'))
     &       qqbZHtop(j)=qqbZHtop(j)*GamHbb/GamHbb0
           bit=bit+beama0(j)*beamb0(k)*qqbZHtop(j)
        endif
      elseif (j < 0) then
        bit=bit*(((L(k)*v2(1))**2+(R(k)*v2(2))**2)*q2413
     &          +((L(k)*v2(2))**2+(R(k)*v2(1))**2)*q1423)
       if(order >=2) then
          qbqZHtop(k)=aveqq*hdecay*(qqb_ZH_VItop(2,1,3,4,p,abs(k))
     &         +qqb_ZH_VIItop(2,1,3,4,p,abs(k)))
           if (FixBrHbb .and. (hdecaymode=='bqba'))
     &       qbqZHtop(k)=qbqZHtop(k)*GamHbb/GamHbb0
           bit=bit+beama0(j)*beamb0(k)*qbqZHtop(k)
        endif
      endif
      if (j==0 .and. order >= 2)  then
         bit=bit+beama0(j)*beamb0(k)*ggZH
      endif


! 44   continue
      xmsq=xmsq+bit

      enddo

      return
      end
