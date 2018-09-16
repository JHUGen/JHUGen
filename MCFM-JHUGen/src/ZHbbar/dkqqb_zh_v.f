      subroutine dkqqb_zh_v(p,msq)
      implicit none
      include 'types.f'
c-----Virtual corrections for Higgs decay to b-bbar
c-----Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  H  + Z
c                           |    |
c                           |    --> e^-(p3)+e^+(p4)
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

      call qqb_zh(p,msq)
      call hbbdecay_v(p,5,6,decayv)

      msq(:,:)=msq(:,:)*decayv

c--- adjust for fixed H->bb BR if necessary
      if (FixBrHbb) then
         msq(:,:)=msq(:,:)*GamHbb0/GamHbb1
      endif

      return

      end
