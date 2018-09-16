!########################################################
!#### next-to-leading order quark beam functions
!#### normalized to as/(2pi) which has been extracted
!########################################################
      subroutine xbeam1bis(ih,zin,xb,QB,btau)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      include 'facscale.f'
      include 'scale.f'
      real(dp), intent(in) ::zin,xb,QB
      real(dp), intent(out) :: btau(-5:5,-1:1)
      real(dp)::
     & p0qiqiz(-1:0),p0qiqiz1(-1:0),p0ggz(-1:0),p0ggz1(-1:0),
     & p0qgz,p0gqz,logB,
     & I1qiqiz(0:2),I1qiqiz1(0:2),I1ggz(0:2),I1ggz1(0:2),I1gqz,
     & p0ijterm,I1ijterm,
     & p0qg,p0gq,I1qig,I1qigz,I1gqi,
     & L0,L1,L01,L11,fxq,fx(-5:5),fx0(-5:5),
     & zb,jaco,plus
      integer :: ih,j
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

!      lB = log(scale/QB)
      L0 = one/(one-zb)
      L1 = log(one-zb)/(one-zb)

! fill arrays for apbit
      call xp0gg(zb,p0ggz)
      call xp0gg(one,p0ggz1)
      p0gqz = p0gq(zb)

      call xp0qiqi(zb,p0qiqiz)
      call xp0qiqi(one,p0qiqiz1)
      p0qgz = p0qg(zb)

! fill arrays for other bit
      call xI1gg(zb,I1ggz)
      call xI1gg(one,I1ggz1)
      I1gqz = I1gqi(zb)

      call xI1qiqi(zb,I1qiqiz)
      call xI1qiqi(one,I1qiqiz1)
      I1qigz = I1qig(zb)

! calculate parton distribution
      call fdist(ih,xb,facscale,fx0)
      call fdist(ih,xb/zb,facscale,fx)
      
! calculate quark+antiquark sum
      fxq = zip  
      do j = 1,nf
         fxq=fxq+fx(j)+fx(-j)
      enddo         

      do j=-nf,nf
      if (j == 0) then
!---------------------------------------------
!----- Calculation of gluon beam function ----
!---------------------------------------------
        p0ijterm =p0ggz(-1)*fx0(0)
     &  +plus(p0ggz(0),fx(0),L0,p0ggz1(0),fx0(0),L01,zb,jaco)
     &  +p0gqz*fxq/zb*jaco

        I1ijterm = I1ggz(0)*fx0(0)
     &  +plus(I1ggz(1),fx(0),L1,I1ggz1(1),fx0(0),L11,zb,jaco)
     &  +I1ggz(2)*fx(0)/zb*jaco
     &  +I1gqz*fxq/zb*jaco

        btau(j,1)= Gag0*fx0(0)
        btau(j,0)= -gBg0/2.0_dp*fx0(0) + 2*p0ijterm
        btau(j,-1)= 2*I1ijterm + 4*p0ijterm*log(scale/facscale)
        
      else
!---------------------------------------------
!----- Calculation of quark beam function ----
!---------------------------------------------
        p0ijterm = p0qiqiz(-1)*fx0(j)
     &   +plus(p0qiqiz(0),fx(j),L0,p0qiqiz1(0),fx0(j),L01,zb,jaco)
     &   +p0qgz*fx(0)/zb*jaco

        I1ijterm = I1qiqiz(0)*fx0(j)
     &   +plus(I1qiqiz(1),fx(j),L1,I1qiqiz1(1),fx0(j),L11,zb,jaco)
     &   +I1qiqiz(2)*fx(j)/zb*jaco
     &   +I1qigz*fx(0)/zb*jaco

        btau(j,1)=Gaq0*fx0(j)
        btau(j,0)=-gBq0/2.0_dp*fx0(j) + 2*p0ijterm
        btau(j,-1)=2*I1ijterm + 4*p0ijterm*log(scale/facscale)

      endif
      
c--- we want result as a distribution in log(taucut/scale)
c--- rather than log(taucut*QB/musq) => extra +log(QB/scale)

      logB=log(QB/scale)
      
      btau(j,-1)=btau(j,-1)+btau(j,0)*logB+btau(j,1)*logB**2/2._dp
      btau(j,0)=btau(j,0)+logB*btau(j,1)

c--- adjust normalization to (as/2/pi) rather than (as/4/pi)
      btau(j,:)=btau(j,:)/two

c      Lt1cut = log(taucut*QB/musq)
c      Lt1cut2 = Lt1cut**2/2.0_dp
c      beam1(j)=half*(btau(-1)+btau(0)*Lt1cut+btau(1)*Lt1cut2)
      enddo
      return
      end


