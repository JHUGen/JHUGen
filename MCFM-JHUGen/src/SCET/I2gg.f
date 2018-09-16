      subroutine xI2gg(z,I2gg)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'

      real(dp), intent(in) :: z
      real(dp) :: I2gg(-1:4)
      real(dp) :: I2ggAz(-1:4),I2ggFz

      call xI2ggA(z,I2ggAz)
      call xI2ggF(z,I2ggFz)

      I2gg(-1:3) = I2ggAz(-1:3)*CA
      I2gg(4) = I2ggAz(4)*CA + I2ggFz*TR*NF

      return
      end



      subroutine xI2ggA(z,I2ggA)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'

      real(dp), intent(in) :: z
      real(dp) :: I2ggA(-1:4)
      real(dp) :: T3z,V3z,U3z,S3z,pggz,pggmz
      real(dp) :: splitT3,splitV3,splitU3,splitS3,Li2

      T3z = splitT3(z)
      V3z = splitV3(z)
      U3z = splitU3(z)
      S3z = splitS3(z)
      pggz = 2*(one-z+z**2)**2/z/(one-z)
      pggmz = 2*(one+z+z**2)**2/(-z)/(one+z)

      I2ggA(3) = CA*2*(one-z+z**2)**2/z

      I2ggA(2) = be0*2*(one-z+z**2)**2/z*(-1/4.0_dp)

      I2ggA(1) = CA*2*(one-z+z**2)**2/z*(2/3.0_dp-pisq)
     &         +be0*2*(one-z+z**2)**2/z*(5/6.0_dp)

      I2ggA(0) = CA*2*(one-z+z**2)**2/z*(-8/9.0_dp+15/2.0_dp*zeta3)
     &         +be0*2*(one-z+z**2)**2/z*(-7/9.0_dp+pisq/12.0_dp)

      I2ggA(-1) = CA*(52/27.0_dp-zeta2+11/360.0_dp*pisq**2)
     &           +be0*(41/27.0_dp-5/24.0_dp*pisq-5/6.0_dp*zeta3)

      I2ggA(4) = CA*(pggz*(V3z-U3z-T3z-log((one-z)/z)*log(one-z)*log(z))
     &     +8*(one+z)*T3z + pggmz*S3z
     &     +(-22/3.0_dp/z+6._dp-6*z+22/3.0_dp*z**2)
     &      *(log((one-z)/z)**2-pisq/3.0_dp)
     &     +(11/3.0_dp/z+15._dp+4*z+44/3.0_dp*z**2)*log(z)**2
     &     +(22/3.0_dp/z+14._dp+8*z+44/3.0_dp*z**2)*(Li2(z)-zeta2)
     &     +(-143/18.0_dp/z+34/3.0_dp-145/12.0_dp*z+143/18.0_dp*z**2)
     &      *log(one-z)
     &     -(4/3.0_dp/(one-z)+149/18.0_dp/z+155/6.0_dp
     &       +101/12.0_dp*z+43/2.0_dp*z**2)*log(z)
     &     -(one-z)*(209/9.0_dp/z+8._dp+403/18.0_dp*z) )
     & +be0*(
     &  pggz*(log(one-z)*log(z)/2.0_dp-log(z)**2/4.0_dp-5/6.0_dp*log(z))
     &  +(13/6.0_dp/z-3/2.0_dp+7/4.0_dp*z-13/6.0_dp*z**2)*log((one-z)/z)
     &   -(one+z)*(Li2(z)+3/4.0_dp*log(z)**2+7/12.0_dp*log(z)-zeta2)
     &   -34/9.0_dp/z+19/6.0_dp-11/3.0_dp*z+77/18.0_dp*z**2 )

      return
      end


      subroutine xI2ggF(z,I2ggF)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'

      real(dp), intent(in) :: z
      real(dp) :: I2ggF
      real(dp) :: T3z,splitT3,Li2

      T3z = splitT3(z)

      I2ggF = CF*(-4*(one+z)*T3z
     &   +(4/3.0_dp/z+one-z-4/3.0_dp*z**2)
     &    *(log((one-z)/z)**2-pisq/3.0_dp)+(5._dp+7*z)/2.0_dp*log(z)**2
     &   +2*(2._dp+3*z)*(Li2(z)+log(z)-zeta2)
     &   +(-14/9.0_dp/z-40/3.0_dp+28/3.0_dp*z+50/9.0_dp*z**2)
     &    *log((one-z)/z)
     &   +23/27.0_dp/z+247/9.0_dp-211/9.0_dp*z-131/27.0_dp*z**2 )

      return
      end


