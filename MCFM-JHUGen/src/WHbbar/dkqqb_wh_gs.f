      subroutine dkqqb_wh_gs(p,msqc)
      implicit none
      include 'types.f'

c---Matrix element squared averaged over initial colors and spins
c for nwz=1
c     q(-p1)+qbar(-p2) -->  H  + W
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+b(p6)+g(p7)
c for nwz=-1
c     q(-p1)+qbar(-p2) -->  H  + W
c                           |    |
c                           |    --> e^-(p3)+nubar(p4)
c                           |
c                           ---> b(p5)+b(p6)+g(p7)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'incldip.f'
      include 'masses.f'
      include 'breit.f'
      include 'hbbparams.f'
      real(dp):: msqc(maxd,-nf:nf,-nf:nf),p(mxpart,4)
      real(dp)::
     & msq57_6(-nf:nf,-nf:nf),msq67_5(-nf:nf,-nf:nf),
     & dummyv(-nf:nf,-nf:nf),sub57_6(4),sub67_5(4),dsubv,
     & oldmass2
      integer:: j,k,nd
      real(dp), parameter:: tiny=1e-6_dp
      external qqb_wh,donothing_gvec

      ndmax=2

      msqc(1:ndmax,:,:)=0._dp
      incldip(1:ndmax)=.true.

      if (mb > tiny) then
c--- subtractions for massive quarks
        qqproc=.true.
        qgproc=.false.
        gqproc=.false.
        ggproc=.false.
        oldmass2=mass2
        mass2=mb
        call dips_mass(1,p,5,7,6,sub57_6,dsubv,msq57_6,dummyv,
     &   qqb_wh,donothing_gvec)
        call dips_mass(2,p,6,7,5,sub67_5,dsubv,msq67_5,dummyv,
     &   qqb_wh,donothing_gvec)
        mass2=oldmass2
      else
c--- subtractions for massless quarks
        call dips(1,p,5,7,6,sub57_6,dsubv,msq57_6,dummyv,
     &   qqb_wh,donothing_gvec)
        call dips(2,p,6,7,5,sub67_5,dsubv,msq67_5,dummyv,
     &   qqb_wh,donothing_gvec)
      endif

      do j=-nf,nf
      do k=-nf,nf
        msqc(1,j,k)=sub57_6(qq)*msq57_6(j,k)*two*Cf
        msqc(2,j,k)=sub67_5(qq)*msq67_5(j,k)*two*Cf
      enddo
      enddo

c--- adjust for fixed H->bb BR if necessary
      if (FixBrHbb) then
        msqc(1:2,:,:)=msqc(1:2,:,:)*GamHbb0/GamHbb1
      endif

      return
      end

