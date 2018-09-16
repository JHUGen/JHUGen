      subroutine lumxmsq_z(p,xx,z1,z2,QB,order,xmsq)
      implicit none
      include 'types.f'
c----Matrix element for W production
C----averaged over initial colours and spins
C For nwz=+1
c     u(-p1)+dbar(-p2)-->W^+(n(p3)+e^+(p4))
C For nwz=-1
c     d(-p1)+ubar(-p2)-->W^-(e^-(p3)+nbar(p4))
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'scet_const.f'
      integer:: j,k,ih1,ih2,m,n,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),fac,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & z1,z2,QB(2),lum0,lum1(-1:1),lum2(-1:3),bit,
     & msq(-nf:nf,-nf:nf),assemble
      complex(dp):: prop,qqb,qbq
      common/density/ih1,ih2
      include 'cplx.h'

      fac=4._dp*esq**2*xn

      call spinoru(4,p,za,zb)

c--   calculate propagator
      fac=aveqq*fac/s(3,4)**2
      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)
c---case dbar-u or ubar-d
      qqb=za(2,3)*zb(4,1)
      qbq=za(1,3)*zb(4,2)

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
      endif

      xmsq=zip
      do j=-nf,nf
      k=-j
      if (j == 0) cycle

      bit=assemble(order,
     & beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     & beama2(j,:),beamb2(k,:),soft1,soft2,hard)

      if (j > 0) then
        bit=bit*fac*(
     &     abs((Q(j)*q1+L(j)*l1*prop)*qqb)**2
     &    +abs((Q(j)*q1+R(j)*r1*prop)*qqb)**2
     &    +abs((Q(j)*q1+L(j)*r1*prop)*qbq)**2
     &    +abs((Q(j)*q1+R(j)*l1*prop)*qbq)**2)
      elseif (j < 0) then
        bit=bit*fac*(
     &     abs((Q(k)*q1+L(k)*l1*prop)*qbq)**2
     &    +abs((Q(k)*q1+R(k)*r1*prop)*qbq)**2
     &    +abs((Q(k)*q1+L(k)*r1*prop)*qqb)**2
     &    +abs((Q(k)*q1+R(k)*l1*prop)*qqb)**2)
      else
        bit=zip
      endif

      xmsq=xmsq+bit

      enddo

      return
      end
