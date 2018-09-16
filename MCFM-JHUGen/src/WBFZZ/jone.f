      subroutine jone(n2,n3,n4,n1,za,zb,zab,zba,j1,jw,j1l)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cmplxmass.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: h12,h34,i1,i2,i3,i4,n1,n2,n3,n4,jdu
      real(dp):: s12,s34,s123,s124,s134,s234,bit,
     & xl,xr,xq
      complex(dp):: zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & j1(4,2,2,2),propz34,propz12,propw12,gmZ(2,2,2),gmZ12(2,2,2),
     & after(4,2,2),before(4,2,2),jw(4,2,2),zab2,WWgmZ(2,2,2),rxw,
     & j1l(4,2,2,2)
C      after(mu,h17,h34),before(mu,h17,h34),jw(mu,jdu,h34)
C---The one Z-current multiplied by i
C---order of indices Lorentz,jdu up or down,
C---quark-line helicity,lepton-line helicity
C---order of indices gmZ(jdu,h2,h34)
C---process 2 
c              2-----<-----------1          2--------<---------1
c                 o         (                    (        o
c                 o         )                    )        o
c                 o         (                    (        o
c                \mu    3--<----4            3--<----4   \mu
c
c                  before                      after
c
c--- j1l represents current corresponding to non-resonant diagrams
c--- such as the one below:
c              2-----<-----------1
c                        (
c                        ) 
c                        (
c                    3--<----4
c                      o 
c                      o 
c                      o 
c                     \mu
c

      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
c--end statement functions

c      if (runstring(1:5) == 'check') then
c      bit=0d0
c      else
      bit=1d0
c      endif
      
C---setting up couplings dependent on whether we are doing 34-line or 56-line
      if (n3+n4 == 7) then
      xl=l1
      xr=r1
      xq=q1
      elseif (n3+n4 == 11) then
      xl=l2
      xr=r2
      xq=q2
      else
      write(6,*) 'Unexpected case jone.f'
      stop
      endif

      s34=s(n3,n4)
      s12=s(n1,n2)
      propz34=cplx1(s34)-czmass2
      propz12=cplx1(s12)-czmass2
      propw12=cplx1(s12)-cwmass2


      do jdu=1,2
      gmZ(jdu,1,1)=(cplx1(Q(jdu)*xq*bit/s34)+cplx1(L(jdu)*xl)/propz34)
      gmZ(jdu,1,2)=(cplx1(Q(jdu)*xq*bit/s34)+cplx1(L(jdu)*xr)/propz34)
      gmZ(jdu,2,1)=(cplx1(Q(jdu)*xq*bit/s34)+cplx1(R(jdu)*xl)/propz34)
      gmZ(jdu,2,2)=(cplx1(Q(jdu)*xq*bit/s34)+cplx1(R(jdu)*xr)/propz34)

      gmZ12(jdu,1,1)=(cplx1(Q(jdu)*xq/s12)+cplx1(L(jdu)*xl)/propz12)
      gmZ12(jdu,1,2)=(cplx1(Q(jdu)*xq/s12)+cplx1(L(jdu)*xr)/propz12)
      gmZ12(jdu,2,1)=(cplx1(Q(jdu)*xq/s12)+cplx1(R(jdu)*xl)/propz12)
      gmZ12(jdu,2,2)=(cplx1(Q(jdu)*xq/s12)+cplx1(R(jdu)*xr)/propz12)

      rxw=sqrt((cone-cxw)/cxw)
      if (jdu == 1) then
      WWgmZ(jdu,1,1)=-(cplx1(bit*xq/s34)+cplx1(xl)*rxw/propz34)
      WWgmZ(jdu,1,2)=-(cplx1(bit*xq/s34)+cplx1(xr)*rxw/propz34)
      WWgmZ(jdu,2,1)=-(cplx1(bit*xq/s34)+cplx1(xl)*rxw/propz34) 
      WWgmZ(jdu,2,2)=-(cplx1(bit*xq/s34)+cplx1(xr)*rxw/propz34)
      endif
      if (jdu == 2) then
      WWgmZ(jdu,1,1)=(cplx1(bit*xq/s34)+cplx1(xl)*rxw/propz34)
      WWgmZ(jdu,1,2)=(cplx1(bit*xq/s34)+cplx1(xr)*rxw/propz34)
      WWgmZ(jdu,2,1)=(cplx1(bit*xq/s34)+cplx1(xl)*rxw/propz34) 
      WWgmZ(jdu,2,2)=(cplx1(bit*xq/s34)+cplx1(xr)*rxw/propz34)
      endif
      enddo

      i1=n1
      i2=n2
      do h34=1,2
         if (h34==1) then 
            i3=n3
            i4=n4
         elseif (h34==2) then
            i3=n4
            i4=n3
         endif

      s134=s(i1,i3)+s(i3,i4)+s(i4,i1)
      s234=s(i2,i3)+s(i3,i4)+s(i4,i2)
      s123=s(i1,i2)+s(i2,i3)+s(i3,i1)
      s124=s(i1,i2)+s(i2,i4)+s(i4,i1)
C---indices before(mu,h12,h34),after(mu,h12,h34)
      before(:,1,h34)=
     & +(zab(i2,:,i1)*za(i1,i3)+zab(i2,:,i4)*za(i4,i3))*zb(i4,i1)/s134
      after(:,1,h34)=
     & -za(i2,i3)*(zb(i4,i2)*zab(i2,:,i1)+zb(i4,i3)*zab(i3,:,i1))/s234
      before(:,2,h34)=
     & +(zba(i2,:,i1)*zb(i1,i4)+zba(i2,:,i3)*zb(i3,i4))*za(i3,i1)/s134
      after(:,2,h34)=
     & -zb(i2,i4)*(za(i3,i2)*zba(i2,:,i1)+za(i3,i4)*zba(i4,:,i1))/s234

      do jdu=1,2
C----Gamma/Z attachment to quark line
      jw(:,jdu,h34)=
     &     gmZ(jdu,1,h34)*before(:,1,h34)
     &    +gmZ(3-jdu,1,h34)*after(:,1,h34)
C----Gamma/Z attachment to exchanged W
     & +(WWgmZ(jdu,1,h34))/propw12
     & *(zab2(i3,i1,i2,i4)*zab(i2,:,i1)-zab2(i2,i3,i4,i1)*zab(i3,:,i4)
     &  -0.5d0*za(i2,i3)*zb(i4,i1)
     & *(zab(i1,:,i1)+zab(i2,:,i2)-zab(i3,:,i3)-zab(i4,:,i4)))
C----non-resonant e-e+ production by exchanged W-line.
      if (h34 == 1) then
      if (((jdu == 1) .and. (xq < 0)) .or.
     &    ((jdu == 2) .and. (xq >-1))) then
      jw(:,jdu,h34)=jw(:,jdu,h34)
     & -(zb(i1,i2)*zab(i2,:,i4)+zb(i1,i3)*zab(i3,:,i4))*za(i3,i2)*bit
     & /(2d0*cxw*propw12*s123)
      elseif (((jdu == 2) .and. (xq < 0)) .or.
     &        ((jdu == 1) .and. (xq >-1))) then
      jw(:,jdu,h34)=jw(:,jdu,h34)
     & +(zab(i3,:,i1)*za(i1,i2)+zab(i3,:,i4)*za(i4,i2))*zb(i1,i4)*bit
     & /(2d0*cxw*propw12*s124)
      endif
      endif

C---add normalization factor
      jw(:,jdu,h34)=rt2/(sqrt(cxw))*jw(:,jdu,h34)

C--Apply couplings for jdu=1,2 to j1
      do h12=1,2
      j1(:,jdu,h12,h34)=
     & 2d0*gmZ(jdu,h12,h34)*(before(:,h12,h34)+after(:,h12,h34))
      enddo

C----non-resonant e-e+ production by exchanged Z/gamma-line.
      j1l(:,jdu,1,h34)=
     & -(zb(i1,i2)*zab(i2,:,i4)+zb(i1,i3)*zab(i3,:,i4))*za(i3,i2)*bit
     & *gmZ12(jdu,1,h34)/s123*2d0
     & +(zab(i3,:,i1)*za(i1,i2)+zab(i3,:,i4)*za(i4,i2))*zb(i1,i4)*bit
     & *gmZ12(jdu,1,h34)/s124*2d0
      j1l(:,jdu,2,h34)=
     & -(zb(i2,i1)*zab(i1,:,i4)+zb(i2,i3)*zab(i3,:,i4))*za(i3,i1)*bit
     & *gmZ12(jdu,2,h34)/s123*2d0
     & +(zab(i3,:,i2)*za(i2,i1)+zab(i3,:,i4)*za(i4,i1))*zb(i2,i4)*bit
     & *gmZ12(jdu,2,h34)/s124*2d0
     
c--- to get opposite helicity, 3<->4 swap needs additional minus sign
      if (h34 == 2) then
        j1l(:,jdu,:,h34)=-j1l(:,jdu,:,h34)
      endif
      
      enddo
      enddo
      
      return
      end
