      subroutine lumxmsq_h(p,xx,z1,z2,QB,order,xmsq)
      implicit none
      include 'types.f'
c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H -->  b(p3)+b(p4))
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'scet_const.f'
      include 'taucut.f'
      integer:: j,k,ih1,ih2,m,n,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),s,fac,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & z1,z2,QB(2),lum0,lum1(-1:1),lum2(-1:3),full(-1:3),bit,
     & msq(-nf:nf,-nf:nf),assemble
      common/density/ih1,ih2
      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &           -p(j,2)*p(k,2)-p(j,3)*p(k,3))

      call gg_h(p,msq)

      call softggbis(order,soft1,soft2)
      call hardgg(s(1,2),musq,hard)

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
      endif

      xmsq=msq(0,0)*assemble(order,
     & beama0(0),beamb0(0),beama1(0,:),beamb1(0,:),
     & beama2(0,:),beamb2(0,:),soft1,soft2,hard)

      return
      end
