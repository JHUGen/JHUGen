      subroutine jonezstrong(n2,n3,n4,n1,za,zb,zab,zba,jglue)
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
      real(dp):: s34,s134,s234,xl,xr,xq
      complex(dp):: zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & jglue(4,2,2,2),propz34,gmZ(2,2,2),
     & after(4,2,2),before(4,2,2)
c--- returns gluon current for  mu -> e^- e^+ qq

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
      write(6,*) 'Unexpected case jonez.f'
      stop
      endif

      s34=s(n3,n4)
      propz34=cplx1(s34)-czmass2

      do jdu=1,2
      gmZ(jdu,1,1)=(cplx1(Q(jdu)*xq/s34)+cplx1(L(jdu)*xl)/propz34)
      gmZ(jdu,1,2)=(cplx1(Q(jdu)*xq/s34)+cplx1(L(jdu)*xr)/propz34)
      gmZ(jdu,2,1)=(cplx1(Q(jdu)*xq/s34)+cplx1(R(jdu)*xl)/propz34)
      gmZ(jdu,2,2)=(cplx1(Q(jdu)*xq/s34)+cplx1(R(jdu)*xr)/propz34)
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
C--Apply couplings for jdu=1,2 to j1
      do h12=1,2
      jglue(:,jdu,h12,h34)=
     & 2d0*gmZ(jdu,h12,h34)*(before(:,h12,h34)+after(:,h12,h34))
      enddo
      enddo
      
      enddo
      
      return
      end
