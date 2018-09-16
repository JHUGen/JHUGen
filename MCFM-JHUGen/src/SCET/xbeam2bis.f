!########################################################
!#### next-to-next-to-leading order beam functions
!#### normalized to (as/(2 pi))^2 which has been extracted
!########################################################
      subroutine xbeam2bis(ih,zin,xb,QB,btau)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      include 'facscale.f'
      include 'scale.f'
      real(dp), intent(in) ::zin,xb,QB
      real(dp), intent(out) :: btau(-5:5,-1:3)
      real(dp)::
     & p0qiqiz(-1:0),p0qiqiz1(-1:0),p0ggz(-1:0),p0ggz1(-1:0),
     & p0qgz,p0gqz,logB,
     & I1qiqiz(0:2),I1qiqiz1(0:2),I1ggz(0:2),I1ggz1(0:2),I1gqz,
     & p0ijterm,I1ijterm,p1ijterm,p0ikp0kjterm,I1ikp0kjterm,
     & p0qg,p0gq,
     & p0gqp0qqz,p0ggp0gqz,
     & pgqpqq,pggpgq,
     & I1gqp0qqz,I1ggp0gqz,
     & Igqpqq,Iggpgq,
     & I1qig,I1qigz,I1gqi,
     & L0,L1,L2,L3,L01,L11,L21,L31,
     & fxq,fxqb,fxqsum,fx(-5:5),fx0(-5:5),zb,jaco,
     & I2ijterm,I1gqp0qgz,p0gqp0qgz,p1gqz,p1gqi,
     & p1qigz,p1qiqbiz,p1qiqjz,p1qiqbjz,
     & p1qig,p1qiqbi,p1qiqj,p1qiqbj,
     & p1ggz(-1:1),p1ggz1(-1:1),p0ggp0ggz(-1:2),p0ggp0ggz1(-1:2),
     & p1qiqiz(-1:1),p1qiqiz1(-1:1),p0qqp0qqz(-1:2),p0qqp0qqz1(-1:2),
     & pgqpqg,pqgpgq,p0qgp0gqz,p0qgp0ggz,p0qqp0qgz,pqgpgg,pqqpqg,
     & I1ggp0ggz(-1:3),I1ggp0ggz1(-1:3),
     & I1qqp0qqz(-1:3),I1qqp0qqz1(-1:3),
     & I1qgp0gqz,Igqpqg,plus,
     & I1qgp0ggz,I1qqp0qgz,Iqgpgg,Iqqpqg,Iqgpgq,
     & I2qiqiz(-1:4),I2qiqiz1(-1:4),I2ggz(-1:4),I2ggz1(-1:4),I2gqz,
     & I2gqi,I2qiqbiz,I2qiqjz,I2qiqbjz,I2qigz,
     & I2qiqbi,I2qiqj,I2qiqbj,I2qig
      integer :: ih,j,k
      real(dp), parameter:: tiny=1.e-15_dp

c--- changing variables so z integral ranges from xb to 1
      zb = (one - xb)*zin + xb
      jaco = abs(one - xb)

c--- catch numerical problems
      if (zb > one-tiny) then
        btau(:,:)=zip
        return
      endif
      
      L01 = log(one - xb)
      L11 = L01**2/2.0_dp
      L21 = L01**3/3.0_dp
      L31 = L01**4/4.0_dp

      L0 = one/(one-zb)
      L1 = log(one-zb)/(one-zb)
      L2 = log(one-zb)**2/(one-zb)
      L3 = log(one-zb)**3/(one-zb)
      
!---------------------------------------------
!----- Components of gluon beam function -----
!---------------------------------------------
! fill arrays for apbit
      call xp0gg(zb,p0ggz)
      call xp0gg(one,p0ggz1)
      p0gqz = p0gq(zb)

! fill arrays for other bit
      call xI1gg(zb,I1ggz)
      call xI1gg(one,I1ggz1)
      I1gqz = I1gqi(zb)

! fill arrays for p1ij
      call xp1gg(zb,p1ggz)
      call xp1gg(one,p1ggz1)
      p1gqz = p1gqi(zb)

! fill arrays for p0ikp0kj
      call xpggpgg(zb,p0ggp0ggz)
      call xpggpgg(one,p0ggp0ggz1)
      p0gqp0qqz = pgqpqq(zb)
      p0ggp0gqz = pggpgq(zb)
! factor of 2*nf for sum over intermediate quarks and antiquarks 
      p0gqp0qgz = 2*nf*pgqpqg(zb)

! fill arrays for I1ikp0kj
      call xIggpgg(zb,I1ggp0ggz)
      call xIggpgg(one,I1ggp0ggz1)
      I1gqp0qqz = Igqpqq(zb)
      I1ggp0gqz = Iggpgq(zb)
! factor of 2*nf for sum over intermediate quarks and antiquarks 
      I1gqp0qgz = 2*nf*Igqpqg(zb)

! fill arrays for I2ij
      call xI2gg(zb,I2ggz)
      call xI2gg(one,I2ggz1)
      I2gqz = I2gqi(zb)

!---------------------------------------------
!----- Components of quark beam function -----
!---------------------------------------------
! fill arrays for apbit
      call xp0qiqi(zb,p0qiqiz)
      call xp0qiqi(one,p0qiqiz1)
      p0qgz = p0qg(zb)

! fill arrays for other bit
      call xI1qiqi(zb,I1qiqiz)
      call xI1qiqi(one,I1qiqiz1)
      I1qigz = I1qig(zb)

! fill arrays for p1ij
      call xp1qiqi(zb,p1qiqiz)
      call xp1qiqi(one,p1qiqiz1)
      p1qiqbiz = p1qiqbi(zb)
      p1qiqjz = p1qiqj(zb)
      p1qiqbjz = p1qiqbj(zb)
      p1qigz = p1qig(zb)

! fill arrays for p0ikp0kj
      call xpqqpqq(zb,p0qqp0qqz)
      call xpqqpqq(one,p0qqp0qqz1)
      p0qgp0gqz = pqgpgq(zb)
      p0qgp0ggz = pqgpgg(zb)
      p0qqp0qgz = pqqpqg(zb)

! fill arrays for I1ikp0kj
      call xIqqpqq(zb,I1qqp0qqz)
      call xIqqpqq(one,I1qqp0qqz1)
      I1qgp0gqz = Iqgpgq(zb)
      I1qgp0ggz = Iqgpgg(zb)
      I1qqp0qgz = Iqqpqg(zb)

! fill arrays for I2ij
      call xI2qiqi(zb,I2qiqiz)
      call xI2qiqi(one,I2qiqiz1)
      I2qiqbiz = I2qiqbi(zb)
      I2qiqjz = I2qiqj(zb)
      I2qiqbjz = I2qiqbj(zb)
      I2qigz = I2qig(zb)

! calculate parton distribution
      call fdist(ih,xb,facscale,fx0)
      call fdist(ih,xb/zb,facscale,fx)

! calculate quark+antiquark sum
      fxq = zip
      fxqb = zip
      do j = 1,nf
         fxq=fxq+fx(j)
         fxqb=fxqb+fx(-j)
      enddo
      fxqsum=fxq+fxqb

      do j=-nf,nf
      if (j == 0) then
!---------------------------------------------
!----- Calculation of gluon beam function ----
!---------------------------------------------
        p0ijterm=p0ggz(-1)*fx0(0)
     &  +plus(p0ggz(0),fx(0),L0,p0ggz1(0),fx0(0),L01,zb,jaco)
     &  +p0gqz*fxqsum/zb*jaco

        I1ijterm = I1ggz(0)*fx0(0)
     &  +plus(I1ggz(1),fx(0),L1,I1ggz1(1),fx0(0),L11,zb,jaco)
     &  +I1ggz(2)*fx(0)/zb*jaco
     &  +I1gqz*fxqsum/zb*jaco

        p1ijterm =p1ggz(-1)*fx0(0)
     &  +plus(p1ggz(0),fx(0),L0,p1ggz1(0),fx0(0),L01,zb,jaco)
     &  +p1ggz(1)*fx(0)/zb*jaco
     &  +p1gqz*fxqsum/zb*jaco

        p0ikp0kjterm =
     &  +p0ggp0ggz(-1)*fx0(0)
     &  +plus(p0ggp0ggz(0),fx(0),L0,p0ggp0ggz1(0),fx0(0),L01,zb,jaco)
     &  +plus(p0ggp0ggz(1),fx(0),L1,p0ggp0ggz1(1),fx0(0),L11,zb,jaco)
     &  +p0ggp0ggz(2)*fx(0)/zb*jaco
     &  +p0gqp0qgz*fx(0)/zb*jaco
     &  +(p0gqp0qqz+p0ggp0gqz)*fxqsum/zb*jaco

        I1ikp0kjterm =
     &  +I1ggp0ggz(-1)*fx0(0)
     &  +plus(I1ggp0ggz(0),fx(0),L0,I1ggp0ggz1(0),fx0(0),L01,zb,jaco)
     &  +plus(I1ggp0ggz(1),fx(0),L1,I1ggp0ggz1(1),fx0(0),L11,zb,jaco)
     &  +plus(I1ggp0ggz(2),fx(0),L2,I1ggp0ggz1(2),fx0(0),L21,zb,jaco)
     &  +I1ggp0ggz(3)*fx(0)/zb*jaco
     &  +I1gqp0qgz*fx(0)/zb*jaco
     &  +(I1gqp0qqz+I1ggp0gqz)*fxqsum/zb*jaco

        I2ijterm = I2ggz(-1)*fx0(0)
     &  +plus(I2ggz(0),fx(0),L0,I2ggz1(0),fx0(0),L01,zb,jaco) 
     &  +plus(I2ggz(1),fx(0),L1,I2ggz1(1),fx0(0),L11,zb,jaco)
     &  +plus(I2ggz(2),fx(0),L2,I2ggz1(2),fx0(0),L21,zb,jaco) 
     &  +plus(I2ggz(3),fx(0),L3,I2ggz1(3),fx0(0),L31,zb,jaco) 
     &  +I2ggz(4)*fx(0)/zb*jaco
     &  +I2gqz*fxqsum/zb*jaco

!---implementation of 1405.1044v2, Eq. (3.3)
!---additional factorization scale dependence in 1505.04794, Eq. (A.20)
      btau(j,3)=Gag0**2/two*fx0(0)
      btau(j,2)=Gag0*(-(three/four*gBg0+be0/two)*fx0(0)+three*p0ijterm)
      btau(j,1)=(Gag1-Gag0**2*pisqo6+gBg0**2/four+be0/two*gBg0)*fx0(0)
     & +two*Gag0*I1ijterm-two*(gBg0+be0)*p0ijterm+four*p0ikp0kjterm
     & +four*Gag0*p0ijterm*log(scale/facscale)
      btau(j,0)=(Gag0**2*zeta3+Gag0*gBg0*pisq/twelve-gBg1/two)*fx0(0)
     & -Gag0*pisq/three*p0ijterm-(gBg0+two*be0)*I1ijterm
     & +four*(I1ikp0kjterm+p1ijterm)
     & +(-gBg0*p0ijterm+four*p0ikp0kjterm)*two*log(scale/facscale)
      btau(j,-1)=four*I2ijterm
     & +four*(I1ikp0kjterm+p1ijterm)*two*log(scale/facscale)
     & +(be0*p0ijterm+two*p0ikp0kjterm)*four*log(scale/facscale)**2
      
      else
!---------------------------------------------
!----- Calculation of quark beam function ----
!---------------------------------------------
      p0ijterm = 
     & +p0qiqiz(-1)*fx0(j)
     & +plus(p0qiqiz(0),fx(j),L0,p0qiqiz1(0),fx0(j),L01,zb,jaco)
     & +p0qgz*fx(0)/zb*jaco

      I1ijterm = 
     & +I1qiqiz(0)*fx0(j)
     & +plus(I1qiqiz(1),fx(j),L1,I1qiqiz1(1),fx0(j),L11,zb,jaco) 
     & +I1qiqiz(2)*fx(j)/zb*jaco
     & +I1qigz*fx(0)/zb*jaco

      p0ikp0kjterm =
     &  +p0qqp0qqz(-1)*fx0(j)
     &  +plus(p0qqp0qqz(0),fx(j),L0,p0qqp0qqz1(0),fx0(j),L01,zb,jaco)
     &  +plus(p0qqp0qqz(1),fx(j),L1,p0qqp0qqz1(1),fx0(j),L11,zb,jaco)
     &  +p0qqp0qqz(2)*fx(j)/zb*jaco
     &  +p0qgp0gqz*fxqsum/zb*jaco
     &  +(p0qgp0ggz+p0qqp0qgz)*fx(0)/zb*jaco

      I1ikp0kjterm = 
     &  +I1qqp0qqz(-1)*fx0(j)
     &  +plus(I1qqp0qqz(0),fx(j),L0,I1qqp0qqz1(0),fx0(j),L01,zb,jaco)
     &  +plus(I1qqp0qqz(1),fx(j),L1,I1qqp0qqz1(1),fx0(j),L11,zb,jaco)
     &  +plus(I1qqp0qqz(2),fx(j),L2,I1qqp0qqz1(2),fx0(j),L21,zb,jaco)
     &  +I1qqp0qqz(3)*fx(j)/zb*jaco
     &  +I1qgp0gqz*fxqsum/zb*jaco
     &  +(I1qgp0ggz+I1qqp0qgz)*fx(0)/zb*jaco

      p1ijterm =p1qiqiz(-1)*fx0(j)
     &  +plus(p1qiqiz(0),fx(j),L0,p1qiqiz1(0),fx0(j),L01,zb,jaco)
     &  +p1qiqiz(1)*fx(j)/zb*jaco
     &  +p1qiqbiz*fx(-j)/zb*jaco
     &  +p1qigz*fx(0)/zb*jaco
      if (j > 0) then
      p1ijterm = p1ijterm
     &  +p1qiqjz*(fxq-fx(j))/zb*jaco
     &  +p1qiqbjz*(fxqb-fx(-j))/zb*jaco
      else
      p1ijterm = p1ijterm
     &  +p1qiqjz*(fxqb-fx(j))/zb*jaco
     &  +p1qiqbjz*(fxq-fx(-j))/zb*jaco
      endif

      I2ijterm = I2qiqiz(-1)*fx0(j)
     &  +plus(I2qiqiz(0),fx(j),L0,I2qiqiz1(0),fx0(j),L01,zb,jaco) 
     &  +plus(I2qiqiz(1),fx(j),L1,I2qiqiz1(1),fx0(j),L11,zb,jaco)
     &  +plus(I2qiqiz(2),fx(j),L2,I2qiqiz1(2),fx0(j),L21,zb,jaco) 
     &  +plus(I2qiqiz(3),fx(j),L3,I2qiqiz1(3),fx0(j),L31,zb,jaco) 
     &  +I2qiqiz(4)*fx(j)/zb*jaco
     &  +I2qiqbiz*fx(-j)/zb*jaco
     &  +I2qigz*fx(0)/zb*jaco
      if (j > 0) then
      I2ijterm = I2ijterm
     &  +I2qiqjz*(fxq-fx(j))/zb*jaco
     &  +I2qiqbjz*(fxqb-fx(-j))/zb*jaco
      else
      I2ijterm = I2ijterm
     &  +I2qiqjz*(fxqb-fx(j))/zb*jaco
     &  +I2qiqbjz*(fxq-fx(-j))/zb*jaco
      endif

!---implementation of 1401.5478v2, Eq. (2.28)
!---additional factorization scale dependence in 1505.04794, Eq. (A.20)
      btau(j,3)=Gaq0**2/two*fx0(j)
      btau(j,2)=Gaq0*(-(three/four*gBq0+be0/two)*fx0(j)+three*p0ijterm)
      btau(j,1)=(Gaq1-Gaq0**2*pisqo6+gBq0**2/four+be0/two*gBq0)*fx0(j)
     & +two*Gaq0*I1ijterm-two*(gBq0+be0)*p0ijterm+four*p0ikp0kjterm
     & +four*Gaq0*p0ijterm*log(scale/facscale)
      btau(j,0)=(Gaq0**2*zeta3+Gaq0*gBq0*pisq/twelve-gBq1/two)*fx0(j)
     & -Gaq0*pisq/three*p0ijterm-(gBq0+two*be0)*I1ijterm
     & +four*(I1ikp0kjterm+p1ijterm)
     & +(-gBq0*p0ijterm+four*p0ikp0kjterm)*two*log(scale/facscale)
      btau(j,-1)=four*I2ijterm
     & +four*(I1ikp0kjterm+p1ijterm)*two*log(scale/facscale)
     & +(be0*p0ijterm+two*p0ikp0kjterm)*four*log(scale/facscale)**2

      endif

c--- we want result as a distribution in log(taucut/scale)
c--- rather than log(taucut*QB/musq) => extra +log(QB/scale)

      logB=log(QB/scale)

      btau(j,-1)=btau(j,-1)
     & +btau(j,0)*logB
     & +btau(j,1)*logB**2/2._dp
     & +btau(j,2)*logB**3/3._dp
     & +btau(j,3)*logB**4/4._dp
      btau(j,0)=btau(j,0)
     & +btau(j,1)*logB+btau(j,2)*logB**2+btau(j,3)*logB**3
      btau(j,1)=btau(j,1)
     & +2._dp*btau(j,2)*logB+3._dp*btau(j,3)*logB**2
      btau(j,2)=btau(j,2)
     & +3._dp*btau(j,3)*logB

c--- adjust normalization to (as/2/pi)^2 rather than (as/4/pi)^2
      btau(j,:)=btau(j,:)/four

c      Lt1cut = log(taucut*QB/musq)
c      Lt1cut2 = Lt1cut**2/2.0_dp
c      Lt1cut3 = Lt1cut**3/3.0_dp
c      Lt1cut4 = Lt1cut**4/4.0_dp

c      beam2(j)=btau(-1)+btau(0)*Lt1cut+btau(1)*Lt1cut2
c     &        +btau(2)*Lt1cut3+btau(3)*Lt1cut4
      enddo
      
      do j=-nf,nf
      do k=-1,3
      if (btau(j,k) == btau(j,k)) then
        continue
      else
        write(6,*) 'NaN in beam function:'
        write(6,*) 'btau(',j,',-1)',btau(j,-1)
        write(6,*) 'btau(',j,', 0)',btau(j, 0)
        write(6,*) 'btau(',j,', 1)',btau(j, 1)
        write(6,*) 'btau(',j,', 2)',btau(j, 2)
        write(6,*) 'btau(',j,', 3)',btau(j, 3)
        write(6,*) 'ih,zin,xb,QB',ih,zin,xb,QB
        call flush(6)
        stop
      endif
      enddo
      enddo
      
      return
      end


