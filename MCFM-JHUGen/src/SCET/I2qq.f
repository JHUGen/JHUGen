      subroutine xI2qiqi(z,I2qiqi)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp), intent(in) :: z
      real(dp) :: I2qiqi(-1:4)
      real(dp) :: I2qqVz(-1:4),I2qqSz

      call xI2qqV(z,I2qqVz)
      call xI2qqS(z,I2qqSz)

      I2qiqi(-1:3) = I2qqVz(-1:3)*CF
      I2qiqi(4) = (I2qqVz(4) + I2qqSz)*CF

      return
      end

!---
      function I2qiqj(z)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp), intent(in) :: z
      real(dp) :: I2qiqj
      real(dp) :: I2qqSz

      call xI2qqS(z,I2qqSz)
      I2qiqj = (I2qqSz)*CF

      return
      end

!---
      function I2qiqbi(z)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp), intent(in) :: z
      real(dp) :: I2qiqbi
c      real(dp) :: I2qqbV,I2qqS
      real(dp) :: I2qqbVz,I2qqSz

      call xI2qqbV(z,I2qqbVz)
      call xI2qqS(z,I2qqSz)

      I2qiqbi = (I2qqbVz + I2qqSz)*CF

      return
      end

!---
      function I2qiqbj(z)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp), intent(in) :: z
      real(dp) :: I2qiqbj
      real(dp) :: I2qqSz

      call xI2qqS(z,I2qqSz)

      I2qiqbj = (I2qqSz)*CF

      return
      end



!---
      subroutine xI2qqV(z,I2qqV)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'

      real(dp), intent(in) :: z
      real(dp) :: I2qqV(-1:4),Li2_TR,omz,opz
      real(dp) :: T3z,V3z,U3z,splitT3,splitV3,splitU3

      omz=one-z
      opz=one+z
      T3z = splitT3(z)
      V3z = splitV3(z)
      U3z = splitU3(z)


      I2qqV(3) = CF*(1+z**2)

      I2qqV(2) = be0*(1+z**2)*(-1/4.0_dp)

      I2qqV(1) = CF*(1+z**2)*(-5/6.0_dp*pisq)
     & +CA*(1+z**2)*(2/3.0_dp-pisq/6.0_dp)
     & +be0*(1+z**2)*5/6.0_dp

      I2qqV(0) = CF*(1+z**2)*4*zeta3
     & +CA*(1+z**2)*(-8/9.0_dp+7/2.0_dp*zeta3)
     & +be0*(1+z**2)*(-7/9.0_dp+pisq/12.0_dp)


      I2qqV(-1) = CF*7/120.0_dp*pisq**2
     & +CA*(52/27.0_dp-pisq/6.0_dp-pisq**2/36.0_dp)
     & +be0*(41/27.0_dp-5/24.0_dp*pisq-5/6.0_dp*zeta3)

      I2qqV(4) = CF*(-2/(1.0_dp-z)*T3z+(1+z**2)/omz*(V3z-2*U3z)
     & +(1/(1.0_dp-z)+z)*(3*Li2_TR(z)+9/4.0_dp*log(z)**2
     & +4*log(z)-pisq/2.0_dp )
     & -omz*(log(omz)*log(z)-pisq/6.0_dp)
     & +(-1/2.0_dp+3/4.0_dp*z)*log(z)**2
     & +(-6+11/2.0_dp*z)*log(omz)
     & -(5/2.0_dp+18*z)*log(z)-(13-15*z)/2.0_dp )
     & +CA*((1+z**2)/(1.0_dp-z)*(U3z-log(omz)*log(omz/z)*log(z)
     & -5/4.0_dp*log(z) )

     & -opz*Li2_TR(z) -3/2.0_dp*z*log(z)**2
     & +(3-5/2.0_dp*z)*log(omz)
     & -(1-11*z)/2.0_dp*log(z)
     & +(7-11*z)/4.0_dp + (1+3*z)/12.0_dp*pisq)
     & +be0*((1+z**2)/omz*(Li2_TR(omz)/2.0_dp+log(omz)*log(z)
     & -5/8.0_dp*log(z)**2-5/4.0_dp*log(z))
     & +omz/2.0_dp*(log(omz)+1/2.0_dp)+z*log(z)/2.0_dp )
      return
      end



!---
      subroutine xI2qqbV(z,I2qqbV)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp), intent(in) :: z
      real(dp) :: I2qqbV,omz,opz
      real(dp) :: S3z,Li2_TR,splitS3

      S3z = splitS3(z)
      omz=one-z
      opz=one+z
      I2qqbV = (2*CF-CA)*((1+z**2)/opz*S3z
     &   -opz/2.0_dp*(Li2_TR(z**2)+2*log(z)*log(opz)-pisq/6.0_dp)
     &   -z*log(z)**2+2*omz*log(omz)
     &   +(3+19*z)/4.0_dp*log(z)+15/4.0_dp*omz)
      return
      end

!---
      subroutine xI2qqS(z,I2qqS)
      implicit none
      include 'types.f'
      include 'constants.f'
c      include 'nf.f'
c      include 'scet_const.f'

      real(dp), intent(in) :: z
      real(dp) :: I2qqS,omz,opz
      real(dp) :: T3z,Li2_TR,splitT3

      omz=one-z
      opz=one+z
      T3z = splitT3(z)

      I2qqS = TR*(-2*opz*T3z
     & -(3+4/3.0_dp/z+5*z+8/3.0_dp*z**2)
     & *(Li2_TR(z)-pisq/6.0_dp)
     & +(1+4/3.0_dp/z-z-4/3.0_dp*z**2)
     & *(log(omz)**2/2.0_dp
     & -log(omz)*log(z) - pisq/6.0_dp)
     & -(13/4.0_dp*opz+10/3.0_dp*z**2)*log(z)**2
     & +(26/9.0_dp/z-11/3.0_dp+17/3.0_dp*z-44/9.0_dp*z**2)
     & *log(omz)
     & +(23/3.0_dp-5/3.0_dp*z+76/9.0_dp*z**2)*log(z)
     & +104/27.0_dp/z-41/18.0_dp+17/18.0_dp*z-68/27.0_dp*z**2)
      return
      end


!---

      function I2qig(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'

      real(dp), intent(in) :: z
      real(dp) :: I2qig,omz,opz
      real(dp) :: T3z,V3z,U3z,S3z,S2z
      real(dp) :: Li2_TR,splitS2,splitS3,splitT3,splitU3,splitV3
      real(dp) :: pqgz,pqgmz

      omz=one-z
      opz=one+z
      T3z = splitT3(z)
      V3z = splitV3(z)
      U3z = splitU3(z)
      S3z = splitS3(z)
      S2z = splitS2(z)
      pqgz = omz**2+z**2
      pqgmz = opz**2+z**2

      I2qig=CF*(-2*omz**2*T3z
     &           +pqgz*(V3z+5/6.0_dp*log(omz)**3
     &                  -log(omz/z)*log(omz)*log(z)
     &                  -5/6.0_dp*pisq*log(omz)
     &                  -pisq/3.0_dp*log(z)+11*zeta3)
     &           -(7/4.0_dp-7*z+6*z**2)*log(omz/z)**2
     &           +(7+4*z)/8.0_dp*log(z)**2
     &           +Li2_TR(z)/2.0_dp
     &           +(13/2.0_dp-23*z+39/2.0_dp*z**2)*log(omz/z)
     &           +(3/2.0_dp+9/4.0_dp*z)*log(z)
     &           -(34-145*z+121*z**2)/4.0_dp
     &           +(3-12*z+10*z**2)*pisq/6.0_dp)
     &      +CA*(-2*(1+4*z)*T3z-pqgz*(U3z-log(omz)**3/6.0_dp
     &           +pisq/6.0_dp*log(omz)-pisq/3.0_dp*log(z)
     &           +7/2.0_dp*zeta3)+pqgmz*S3z-z*opz*S2z
     &           -2*z*omz*log(omz)*log(z)
     &           +(2/3.0_dp/z+1/2.0_dp+3*z-25*z**2/6.0_dp)
     &            *log(omz/z)**2
     &           -(15/4.0_dp+2/3.0_dp/z+2*z+47/3.0_dp*z**2)
     &            *log(z)**2
     &           -(4/3.0_dp/z+3+8*z+44/3.0_dp*z**2)*Li2_TR(z)
     &           +(26/9.0_dp/z+4+31/2.0_dp*z+50/9.0_dp*z**2)
     &            *log(omz)
     &           -(49/6.0_dp+4/3.0_dp*z+323/18.0_dp*z**2)
     &            *log(omz/z)
     &           + 104/27.0_dp/z - 19/36.0_dp + 40*z/9.0_dp
     &           - 947/108.0_dp*z**2 + (2+z+24*z**2)/6.0_dp*pisq)

       I2qig = I2qig*TR

       return
       end


      function I2gqi(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'

      real(dp), intent(in) :: z
      real(dp) :: I2gqi,omz,opz
      real(dp) :: T3z,V3z,U3z,S3z,S2z
      real(dp) :: Li2_TR,splitS2,splitS3,splitT3,splitU3,splitV3
      real(dp) :: pgqz,pgqmz

      omz=one-z
      opz=one+z
      T3z = splitT3(z)
      V3z = splitV3(z)
      U3z = splitU3(z)
      S3z = splitS3(z)
      S2z = splitS2(z)
      pgqz = (1+omz**2)/z
      pgqmz = (1+opz**2)/(-z)

      I2gqi = CA*(pgqz*(V3z-U3z-T3z
     &            +5/6.0_dp*log(omz)**3-2/3.0_dp*pisq*log(omz)
     &            -pisq/3.0_dp*log(z)+45/6.0_dp*zeta3)
     &      +2*(4+z)*T3z+pgqmz*S3z-z/2.0_dp*(S2z-pisq/2.0_dp)
     &      +(-31/6.0_dp/z+4+2*z+2/3.0_dp*z**2)
     &       *(log(omz/z)**2-pisq/3.0_dp)
     &      +(11/3.0_dp/z+15+11/4.0_dp*z+8/3.0_dp*z**2)*log(z)**2
     &      +(22/3.0_dp/z+14+5*z+8/3.0_dp*z**2)*(Li2_TR(z)-pisq/6.0_dp)
     &      +(-181/18.0_dp/z+37/3.0_dp-6*z+44/9.0_dp*z**2)*log(omz)
     &      -(43/6.0_dp/z+51/2.0_dp+13/6.0_dp*z+88/9.0_dp*z**2)*log(z)
     &      -2351/108.0_dp/z+101/6.0_dp-83/36.0_dp*z+152/27.0_dp*z**2 )
     &  +CF*(pgqz*(log(omz)**3/6.0_dp-log(omz/z)*log(omz)*log(z)
     &             -pisq/3.0_dp*log(omz/z)) - (2-z)*T3z
     &      +(-9/2.0_dp/z+11/2.0_dp-3*z)*(log(omz/z)**2-pisq/3.0_dp)
     &      +(-3/z+3-5/2.0_dp*z)*(log(omz)*log(z)+pisq/6.0_dp)
     &      +(9/2.0_dp/z-7+9/8.0_dp*z)*log(z)**2
     &      -(6+5*z)/2.0_dp*(Li2_TR(z)-pisq/6.0_dp)
     &      +(21/2.0_dp/z-75/6.0_dp)*log(omz)
     &      +(-15/6.0_dp/z+9/4.0_dp+15/4.0_dp*z)*log(z)
     &      -43/4.0_dp/z+16-3*z )
     &  +be0*(pgqz*(log(omz)**2/4.0_dp-log(z)**2/2.0_dp
     &        +5/6.0_dp*log(omz)-5/3.0_dp*log(z)-14/9.0_dp)
     &       +z/2.0_dp*(log(omz)+5/3.0_dp) )

      I2gqi = I2gqi*CF

 !     write(*,*)"I2gqi:",I2gqi,"z:",z
 !     pause

      return

      end function I2gqi

