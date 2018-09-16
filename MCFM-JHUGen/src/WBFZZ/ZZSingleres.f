      subroutine ZZSingleres(n1,n2,n3,n4,n5,n6,n7,n8,za,zb,ZZ,WWp,WWm)
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

      real(dp):: s17,s28,s56,t3,
     & s356,s456,s137,s147,s238,s248,xl1,xr1,xq1,xl2,xr2,xq2
      complex(dp):: zab2,ZZ17(2,2,2),ZZ28(2,2,2),
     & prop28,prop17,prop56,propw28,propw17,
     & ZZ(2,2,2,2,2,2),ZZ56(2,2),lZ56(0:1,2,2),
     & srma(2,2,2,2),srmb(2,2,2,2),srmm(2,2,2,2),
     & srpa(2,2,2,2),srpb(2,2,2,2),srpm(2,2,2,2),
     & WWm(2),WWp(2)
      integer:: h17,h28,h34,h56,i1,i2,i3,i4,i5,i6,i7,i8,
     & n1,n2,n3,n4,n5,n6,n7,n8,jdu1,jdu2
C---begin statement functions
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t3(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)
C---end statement functions

C---setting up couplings dependent on whether we are doing 34-line or 56-line
      if (n3+n4 == 7) then
      xl1=l1
      xr1=r1
      xq1=q1
      xl2=l2
      xr2=r2
      xq2=q2
      elseif (n3+n4 == 11) then
      xl1=l2
      xr1=r2
      xq1=q2
      xl2=l1
      xr2=r1
      xq2=q1
      else
      write(6,*) 'Unexpected case ZZSingleres.f'
      stop
      endif

      s17=s(n1,n7)
      s28=s(n2,n8)
      s56=s(n5,n6)

      s356=t3(n3,n5,n6)
      s456=t3(n4,n5,n6)
      s137=t3(n1,n3,n7)
      s147=t3(n1,n4,n7)
      s238=t3(n2,n3,n8)
      s248=t3(n2,n4,n8)

      prop17=cplx1(s17)-czmass2
      prop28=cplx1(s28)-czmass2
      prop56=cplx1(s56)-czmass2

      propw17=cplx1(s17)-cwmass2
      propw28=cplx1(s28)-cwmass2

      do jdu1=1,2
      ZZ17(jdu1,1,1)=cplx1(Q(jdu1)*xq1/s17)+cplx1(L(jdu1)*xl1)/prop17
      ZZ17(jdu1,1,2)=cplx1(Q(jdu1)*xq1/s17)+cplx1(L(jdu1)*xr1)/prop17
      ZZ17(jdu1,2,1)=cplx1(Q(jdu1)*xq1/s17)+cplx1(R(jdu1)*xl1)/prop17
      ZZ17(jdu1,2,2)=cplx1(Q(jdu1)*xq1/s17)+cplx1(R(jdu1)*xr1)/prop17
      ZZ28(jdu1,1,1)=cplx1(Q(jdu1)*xq1/s28)+cplx1(L(jdu1)*xl1)/prop28
      ZZ28(jdu1,1,2)=cplx1(Q(jdu1)*xq1/s28)+cplx1(L(jdu1)*xr1)/prop28
      ZZ28(jdu1,2,1)=cplx1(Q(jdu1)*xq1/s28)+cplx1(R(jdu1)*xl1)/prop28
      ZZ28(jdu1,2,2)=cplx1(Q(jdu1)*xq1/s28)+cplx1(R(jdu1)*xr1)/prop28
      enddo

      ZZ56(1,1)=cplx1(xq1*xq2/s56)+cplx1(xl1*xl2)/prop56
      ZZ56(1,2)=cplx1(xq1*xq2/s56)+cplx1(xl1*xr2)/prop56
      ZZ56(2,1)=cplx1(xq1*xq2/s56)+cplx1(xr1*xl2)/prop56
      ZZ56(2,2)=cplx1(xq1*xq2/s56)+cplx1(xr1*xr2)/prop56

      lZ56(0,1,1)=cplx1(ln*xl2)/prop56
      lZ56(0,1,2)=cplx1(ln*xr2)/prop56
      lZ56(0,2,1)=cplx1(rn*xl2)/prop56
      lZ56(0,2,2)=cplx1(rn*xr2)/prop56
      lZ56(1,1,1)=cplx1(qe*xq2/s56)+cplx1(le*xl2)/prop56
      lZ56(1,1,2)=cplx1(qe*xq2/s56)+cplx1(le*xr2)/prop56
      lZ56(1,2,1)=cplx1(qe*xq2/s56)+cplx1(re*xl2)/prop56
      lZ56(1,2,2)=cplx1(qe*xq2/s56)+cplx1(re*xr2)/prop56

      i3=n3
      i4=n4
      do h17=1,2
         if (h17==1) then
            i7=n7
            i1=n1
         elseif (h17==2) then
            i7=n1
            i1=n7
         endif
      do h28=1,2
         if (h28==1) then
            i8=n8
            i2=n2
         elseif (h28==2) then
            i8=n2
            i2=n8
         endif
      do h56=1,2
         if (h56==1) then
            i5=n5
            i6=n6
         elseif (h56==2) then
            i5=n6
            i6=n5
         endif


C---  Id,srmbl=-8*e^6/s356/s248
C---   *zab2(i8,i2,i4,i1)*za(i3,i5)*zab2(i7,i3,i5,i6)*zb(i2,i4);
      srmb(h17,h28,1,h56)=-8d0/s356/s248
     & *zab2(i8,i2,i4,i1)*za(i3,i5)*zab2(i7,i3,i5,i6)*zb(i2,i4)
C---  Id,srmml=-8*e^6/s137/s248
C---   *zab2(i5,i3,i7,i1)*za(i3,i7)*zab2(i8,i2,i4,i6)*zb(i2,i4);
      srmm(h17,h28,1,h56)=-8d0/s137/s248
     & *zab2(i5,i3,i7,i1)*za(i3,i7)*zab2(i8,i2,i4,i6)*zb(i2,i4)
C---  Id,srmal=-8*e^6/s456/s137
C---   *zab2(i8,i3,i7,i1)*za(i3,i7)*zab2(i5,i4,i6,i2)*zb(i6,i4);
      srma(h17,h28,1,h56)=-8d0/s456/s137
     & *zab2(i8,i3,i7,i1)*za(i3,i7)*zab2(i5,i4,i6,i2)*zb(i6,i4)

C---  Id,srpbl=-8*e^6/s356/s147
C---   *zab2(i8,i3,i5,i6)*za(i3,i5)*zab2(i7,i1,i4,i2)*zb(i1,i4);
      srpb(h17,h28,1,h56)=-8d0/s356/s147
     & *zab2(i8,i3,i5,i6)*za(i3,i5)*zab2(i7,i1,i4,i2)*zb(i1,i4)

C---  Id,srpml=-8*e^6/s147/s238
C---   *zab2(i7,i1,i4,i6)*za(i3,i8)*zab2(i5,i3,i8,i2)*zb(i1,i4);
      srpm(h17,h28,1,h56)=-8d0/s147/s238
     & *zab2(i7,i1,i4,i6)*za(i3,i8)*zab2(i5,i3,i8,i2)*zb(i1,i4)

C---  Id,srpal=-8*e^6/s456/s238
C---   *zab2(i7,i3,i8,i2)*za(i8,i3)*zab2(i5,i4,i6,i1)*zb(i4,i6);
      srpa(h17,h28,1,h56)=-8d0/s456/s238
     & *zab2(i7,i3,i8,i2)*za(i8,i3)*zab2(i5,i4,i6,i1)*zb(i4,i6)


C---  Id,srmbr=-8*e^6/s356/s248
C---   *zab2(i7,i4,i8,i2)*za(i8,i4)*zab2(i5,i3,i6,i1)*zb(i3,i6);
      srmb(h17,h28,2,h56)=-8d0/s356/s248
     & *zab2(i7,i4,i8,i2)*za(i8,i4)*zab2(i5,i3,i6,i1)*zb(i3,i6)

C---  Id,srmmr=-8*e^6/s137/s248
C---   *za(i4,i8)*zab2(i7,i1,i3,i6)*zab2(i5,i4,i8,i2)*zb(i1,i3);
      srmm(h17,h28,2,h56)=-8d0/s137/s248
     & *za(i4,i8)*zab2(i7,i1,i3,i6)*zab2(i5,i4,i8,i2)*zb(i1,i3)

C---  Id,srmar=-8*e^6/s137/s456
C---   *za(i4,i5)*zab2(i7,i1,i3,i2)*zab2(i8,i4,i5,i6)*zb(i1,i3);
      srma(h17,h28,2,h56)=-8d0/s137/s456
     & *za(i4,i5)*zab2(i7,i1,i3,i2)*zab2(i8,i4,i5,i6)*zb(i1,i3)

C---  Id,srpbr=-8*e^6/s356/s147
C---   *zab2(i8,i4,i7,i1)*za(i7,i4)*zab2(i5,i3,i6,i2)*zb(i3,i6);
      srpb(h17,h28,2,h56)=-8d0/s356/s147
     & *zab2(i8,i4,i7,i1)*za(i7,i4)*zab2(i5,i3,i6,i2)*zb(i3,i6)

C---  Id,srpmr=-8*e^6/s147/s238
C---   *za(i4,i7)*zab2(i8,i2,i3,i6)*zab2(i5,i4,i7,i1)*zb(i2,i3);
      srpm(h17,h28,2,h56)=-8d0/s147/s238
     & *za(i4,i7)*zab2(i8,i2,i3,i6)*zab2(i5,i4,i7,i1)*zb(i2,i3)

C---  Id,srpar=-8*e^6/s238/s456
C---   *zab2(i8,i2,i3,i1)*za(i4,i5)*zab2(i7,i4,i5,i6)*zb(i2,i3);
      srpa(h17,h28,2,h56)=-8d0/s238/s456
     & *zab2(i8,i2,i3,i1)*za(i4,i5)*zab2(i7,i4,i5,i6)*zb(i2,i3)
      enddo
      enddo
      enddo

      do h17=1,2
      do h28=1,2
      do h34=1,2
      do h56=1,2
      do jdu1=1,2
      do jdu2=1,2
      ZZ(jdu1,jdu2,h17,h28,h34,h56)=
     & ZZ17(jdu1,h17,h34)*ZZ28(jdu2,h28,h34)*ZZ56(h34,h56)
     & *(srmb(h17,h28,h34,h56)
     &  +srmm(h17,h28,h34,h56)
     &  +srma(h17,h28,h34,h56)
     &  +srpb(h17,h28,h34,h56)
     &  +srpm(h17,h28,h34,h56)
     &  +srpa(h17,h28,h34,h56))
      enddo
      enddo

      if (xq1 < 0) then
      if ((h17==1).and.(h28==1).and.(h34==1)) then
      WWm(h56)=0.25d0/(cxw**2*propw17*propw28)
     & *(srmb(h17,h28,h34,h56)*lZ56(1,h34,h56)
     &  +srmm(h17,h28,h34,h56)*lZ56(0,h34,h56)
     &  +srma(h17,h28,h34,h56)*lZ56(1,h34,h56))

      WWp(h56)=0.25d0/(cxw**2*propw17*propw28)
     & *(srpb(h17,h28,h34,h56)*lZ56(1,h34,h56)
     &  +srpm(h17,h28,h34,h56)*lZ56(0,h34,h56)
     &  +srpa(h17,h28,h34,h56)*lZ56(1,h34,h56))
      endif

      else

      if ((h17==1).and.(h28==1).and.(h34==1)) then
      WWm(h56)=0.25d0/(cxw**2*propw17*propw28)
     & *(srpb(h17,h28,h34,h56)*lZ56(0,h34,h56)
     &  +srpm(h17,h28,h34,h56)*lZ56(1,h34,h56)
     &  +srpa(h17,h28,h34,h56)*lZ56(0,h34,h56))

      WWp(h56)=0.25d0/(cxw**2*propw17*propw28)
     & *(srmb(h17,h28,h34,h56)*lZ56(0,h34,h56)
     &  +srmm(h17,h28,h34,h56)*lZ56(1,h34,h56)
     &  +srma(h17,h28,h34,h56)*lZ56(0,h34,h56))
      endif

      endif

      enddo
      enddo
      enddo
      enddo

      return
      end
