      subroutine jtwo3456(n2,n3,n4,n5,n6,n1,za,zb,zab,zba,j2,j2w)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: h1,h34,h56,i1,i2,i3,i4,i5,i6,n1,n2,n3,n4,n5,n6,jdu,al,nu

      real(dp):: s34,s56,s234,s156,s356,s456,s3456,
     & s23456,s13456,t3,t4,t5,xl1,xr1,xq1,xl2,xr2,xq2

      complex(dp):: zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & j2(4,2,2,2,2),zab2,zba2,prop34,prop56,prop3456,gmZ1(2,2,2),
     & gmZ2(2,2,2),gmZ3(2,2,2),gml56(2,2),
     & before(4,2,2,2),straddle(4,2,2,2),after(4,2,2,2),j2w(4,2,2,2),
     & before1(4,2,2,2),after1(4,2,2,2),
     & j4l(4,2,2),jqb(4,4,2),jqa(4,4,2)
C---  The one Z-current multiplied by i
C---  order of indices:-
C---- Lorentz,jdu,quark-line helicity,lepton34 helicity,lepton56 helicity
C---  jdu=flavour of quark
C---  j2w(mu,jdu,h34,h56)

C--- statement functions
      t3(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t4(i1,i2,i3,i4)=s(i1,i2)+s(i1,i3)+s(i1,i4)
     &               +s(i2,i3)+s(i2,i4)+s(i3,i4)
      t5(i1,i2,i3,i4,i5)=
     & +s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i1,i5)
     &          +s(i2,i3)+s(i2,i4)+s(i2,i5)
     &                   +s(i3,i4)+s(i3,i5)
     &                            +s(i4,i5)
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      zba2(i1,i2,i3,i4)=zb(i1,i2)*za(i2,i4)+zb(i1,i3)*za(i3,i4)
C--- end statement functions

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
      write(6,*) 'Unexpected case jtwo3456.f'
      stop
      endif

      s34=s(n3,n4)
      prop34=cplx2(s34-zmass**2,zmass*zwidth)
      s56=s(n5,n6)
      prop56=cplx2(s56-zmass**2,zmass*zwidth)
      s3456=t4(n3,n4,n5,n6)
      prop3456=cplx2(s3456-zmass**2,zmass*zwidth)

      do jdu=1,2
      gmZ1(jdu,1,1)=cplx1(Q(jdu)*xq1/s34)+cplx1(L(jdu)*xl1)/prop34
      gmZ1(jdu,1,2)=cplx1(Q(jdu)*xq1/s34)+cplx1(L(jdu)*xr1)/prop34
      gmZ1(jdu,2,1)=cplx1(Q(jdu)*xq1/s34)+cplx1(R(jdu)*xl1)/prop34
      gmZ1(jdu,2,2)=cplx1(Q(jdu)*xq1/s34)+cplx1(R(jdu)*xr1)/prop34

      gmZ2(jdu,1,1)=cplx1(Q(jdu)*xq2/s56)+cplx1(L(jdu)*xl2)/prop56
      gmZ2(jdu,1,2)=cplx1(Q(jdu)*xq2/s56)+cplx1(L(jdu)*xr2)/prop56
      gmZ2(jdu,2,1)=cplx1(Q(jdu)*xq2/s56)+cplx1(R(jdu)*xl2)/prop56
      gmZ2(jdu,2,2)=cplx1(Q(jdu)*xq2/s56)+cplx1(R(jdu)*xr2)/prop56

      gmZ3(jdu,1,1)=cplx1(Q(jdu)*xq1/s3456)+cplx1(L(jdu)*xl1)/prop3456
      gmZ3(jdu,1,2)=cplx1(Q(jdu)*xq1/s3456)+cplx1(L(jdu)*xr1)/prop3456
      gmZ3(jdu,2,1)=cplx1(Q(jdu)*xq1/s3456)+cplx1(R(jdu)*xl1)/prop3456
      gmZ3(jdu,2,2)=cplx1(Q(jdu)*xq1/s3456)+cplx1(R(jdu)*xr1)/prop3456
      enddo

      gml56(1,1)=cplx1(xq1*xq2/s56)+cplx1(xl1*xl2)/prop56
      gml56(1,2)=cplx1(xq1*xq2/s56)+cplx1(xl1*xr2)/prop56
      gml56(2,1)=cplx1(xq1*xq2/s56)+cplx1(xr1*xl2)/prop56
      gml56(2,2)=cplx1(xq1*xq2/s56)+cplx1(xr1*xr2)/prop56

      i1=n1
      i2=n2

      s234=t3(n2,n3,n4)
      s156=t3(n1,n5,n6)
      s356=t3(n3,n5,n6)
      s456=t3(n4,n5,n6)
      s23456=t5(n2,n3,n4,n5,n6)
      s13456=t5(n1,n3,n4,n5,n6)

      do nu=1,4
      do al=1,4
      jqb(nu,al,1)=
     & -(zab(n2,nu,n1)*zab(n1,al,n1)
     &  +zab(n2,nu,n3)*zab(n3,al,n1)
     &  +zab(n2,nu,n4)*zab(n4,al,n1)
     &  +zab(n2,nu,n5)*zab(n5,al,n1)
     &  +zab(n2,nu,n6)*zab(n6,al,n1))/s13456

      jqa(nu,al,1)=
     & +(zab(n2,al,n2)*zab(n2,nu,n1)
     &  +zab(n2,al,n3)*zab(n3,nu,n1)
     &  +zab(n2,al,n4)*zab(n4,nu,n1)
     &  +zab(n2,al,n5)*zab(n5,nu,n1)
     &  +zab(n2,al,n6)*zab(n6,nu,n1))/s23456

      jqb(nu,al,2)=
     & -(zba(n2,nu,n1)*zba(n1,al,n1)
     &  +zba(n2,nu,n3)*zba(n3,al,n1)
     &  +zba(n2,nu,n4)*zba(n4,al,n1)
     &  +zba(n2,nu,n5)*zba(n5,al,n1)
     &  +zba(n2,nu,n6)*zba(n6,al,n1))/s13456

      jqa(nu,al,2)=
     & +(zba(n2,al,n2)*zba(n2,nu,n1)
     &  +zba(n2,al,n3)*zba(n3,nu,n1)
     &  +zba(n2,al,n4)*zba(n4,nu,n1)
     &  +zba(n2,al,n5)*zba(n5,nu,n1)
     &  +zba(n2,al,n6)*zba(n6,nu,n1))/s23456

      enddo
      enddo
      do h56=1,2
         if (h56==1) then 
            i5=n5
            i6=n6
         elseif (h56==2) then
            i5=n6
            i6=n5
         endif

      i3=n3
      i4=n4
      j4l(:,1,h56)=
     & +za(i3,i5)/s356*(zb(i6,i3)*zab(i3,:,i4)+zb(i6,i5)*zab(i5,:,i4))
     & -zb(i6,i4)/s456*(zab(i3,:,i4)*za(i4,i5)+zab(i3,:,i6)*za(i6,i5))
      j4l(:,2,h56)=
     & +zb(i3,i6)/s356*(za(i5,i3)*zba(i3,:,i4)+za(i5,i6)*zba(i6,:,i4))
     & -za(i5,i4)/s456*(zba(i3,:,i4)*zb(i4,i6)+zba(i3,:,i5)*zb(i5,i6))

      do h34=1,2
         if (h34==1) then 
            i3=n3
            i4=n4
        elseif (h34==2) then
            i3=n4
            i4=n3
        endif
      after(:,1,h34,h56)=
     & +za(i2,i3)*zba2(i4,i2,i3,i5)/(s234*s23456)
     & *(zb(i6,i2)*zab(i2,:,i1)
     &  +zb(i6,i3)*zab(i3,:,i1)
     &  +zb(i6,i4)*zab(i4,:,i1)
     &  +zb(i6,i5)*zab(i5,:,i1))

      straddle(:,1,h34,h56)=
     & -za(i2,i3)*zb(i6,i1)/(s234*s156)
     & *(zb(i4,i2)*zab(i2,:,i1)*za(i1,i5)
     &  +zb(i4,i3)*zab(i3,:,i1)*za(i1,i5)
     &  +zb(i4,i2)*zab(i2,:,i6)*za(i6,i5)
     &  +zb(i4,i3)*zab(i3,:,i6)*za(i6,i5))

      before(:,1,h34,h56)=
     & +zba2(i4,i1,i6,i5)*zb(i6,i1)/(s13456*s156)
     & *(zab(i2,:,i1)*za(i1,i3)
     &  +zab(i2,:,i4)*za(i4,i3)
     &  +zab(i2,:,i5)*za(i5,i3)
     &  +zab(i2,:,i6)*za(i6,i3))

      before1(:,1,h34,h56)=
     & +jqb(:,4,1)*j4l(4,h34,h56)
     & -jqb(:,1,1)*j4l(1,h34,h56)
     & -jqb(:,2,1)*j4l(2,h34,h56)
     & -jqb(:,3,1)*j4l(3,h34,h56)

      after1(:,1,h34,h56)=
     & +jqa(:,4,1)*j4l(4,h34,h56)
     & -jqa(:,1,1)*j4l(1,h34,h56)
     & -jqa(:,2,1)*j4l(2,h34,h56)
     & -jqa(:,3,1)*j4l(3,h34,h56)

      after(:,2,h34,h56)=
     & +zb(i2,i4)*zab2(i3,i2,i4,i6)/(s234*s23456)
     & *(za(i5,i2)*zba(i2,:,i1)
     &  +za(i5,i4)*zba(i4,:,i1)
     &  +za(i5,i3)*zba(i3,:,i1)
     &  +za(i5,i6)*zba(i6,:,i1))

      straddle(:,2,h34,h56)=
     & -zb(i2,i4)*za(i5,i1)/(s234*s156)
     & *(za(i3,i2)*zba(i2,:,i1)*zb(i1,i6)
     &  +za(i3,i4)*zba(i4,:,i1)*zb(i1,i6)
     &  +za(i3,i2)*zba(i2,:,i5)*zb(i5,i6)
     &  +za(i3,i4)*zba(i4,:,i5)*zb(i5,i6))

      before(:,2,h34,h56)=
     & +zab2(i3,i1,i5,i6)*za(i5,i1)/(s13456*s156)
     & *(zba(i2,:,i1)*zb(i1,i4)
     &  +zba(i2,:,i3)*zb(i3,i4)
     &  +zba(i2,:,i6)*zb(i6,i4)
     &  +zba(i2,:,i5)*zb(i5,i4))

      before1(:,2,h34,h56)=
     & +jqb(:,4,2)*j4l(4,h34,h56)
     & -jqb(:,1,2)*j4l(1,h34,h56)
     & -jqb(:,2,2)*j4l(2,h34,h56)
     & -jqb(:,3,2)*j4l(3,h34,h56)

      after1(:,2,h34,h56)=
     & +jqa(:,4,2)*j4l(4,h34,h56)
     & -jqa(:,1,2)*j4l(1,h34,h56)
     & -jqa(:,2,2)*j4l(2,h34,h56)
     & -jqa(:,3,2)*j4l(3,h34,h56)

      do jdu=1,2
      do h1=1,2
      j2(:,jdu,h1,h34,h56)=
     & 4d0*gmZ1(jdu,h1,h34)*gmZ2(jdu,h1,h56)*
     & (before(:,h1,h34,h56)+straddle(:,h1,h34,h56)+after(:,h1,h34,h56))
     & +2d0*gmZ3(jdu,h1,h34)*gml56(h34,h56)
     & *(before1(:,h1,h34,h56)+after1(:,h1,h34,h56))
      enddo
      j2w(:,jdu,h34,h56)=
     & +4d0*gmZ1(jdu,1,h34)*gmZ2(jdu,1,h56)*before(:,1,h34,h56)
     & +4d0*gmZ1(3-jdu,1,h34)*gmZ2(jdu,1,h56)*straddle(:,1,h34,h56)
     & +4d0*gmZ1(3-jdu,1,h34)*gmZ2(3-jdu,1,h56)*after(:,1,h34,h56)
     & +2d0*gmZ3(jdu,1,h34)*before1(:,1,h34,h56)*gml56(h34,h56)
     & +2d0*gmZ3(3-jdu,1,h34)*after1(:,1,h34,h56)*gml56(h34,h56)
      enddo

      enddo
      enddo

      return
      end
