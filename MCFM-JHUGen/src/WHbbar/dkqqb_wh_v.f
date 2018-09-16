      subroutine dkqqb_wh_v(p,msq)
      implicit none
      include 'types.f'
c-----Virtual corrections for Higgs decay to b-bbar
c-----Matrix element squared averaged over initial colors and spins
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
c----Formula taken from Braaten and Leveille, PR D22, 715, 1980

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'hbbparams.f'
      integer:: j,k
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),
     & s,prop,fac,qqbWH,qbqWH,s56,hdecay,decayv
C----statement function
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))
C----end statement function

      call qqb_wh(p,msq)
      call hbbdecay_v(p,5,6,decayv)

      msq(:,:)=msq(:,:)*decayv

c--- adjust for fixed H->bb BR if necessary
      if (FixBrHbb) then
        msq(:,:)=msq(:,:)*GamHbb0/GamHbb1
      endif

      return




      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo
      s56=s(5,6)+2*mb**2

c---calculate the 2 W propagators
      prop=     ((s(1,2)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

      fac=xn*gwsq**3*wmass**2/prop
      call hbbdecay_v(p,5,6,hdecay)
      hdecay=hdecay/((s56-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay

      qqbWH=aveqq*fac*s(1,4)*s(2,3)
      qbqWH=aveqq*fac*s(2,4)*s(1,3)

      do j=-nf,nf
      do k=-nf,nf
      if ((j > 0) .and. (k < 0)) msq(j,k)=Vsq(j,k)*qqbWH
      if ((j < 0) .and. (k > 0)) msq(j,k)=Vsq(j,k)*qbqWH
      enddo
      enddo

      return
      end

