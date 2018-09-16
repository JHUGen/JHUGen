!--- return p0gg(0) the coefficient of L0(1-z)
!---        p0gg(-1) the coefficient of delta(1-z)
       subroutine xp0gg(z,p0gg)
       implicit none
       include 'types.f'
       include 'constants.f'
       include 'nf.f'
       include 'scet_const.f'
       real(dp), intent(in) :: z
       real(dp) :: p0gg(-1:0)
       p0gg(-1) = be0/2.0_dp  !delta function
       p0gg(0) = 2.0_dp*(1.0_dp-z+z**2)**2/z*CA  !coeff of 1/(1-z)
       return
       end

!---
       function p0gq(z)
       implicit none
       include 'types.f'
       include 'constants.f'
       real(dp), intent(in) :: z
       real(dp) :: p0gq
       p0gq = CF*(1.0_dp+(1.0_dp-z)**2)/z
       return
       end


!---
       function p0qg(z)
       implicit none
       include 'types.f'
       include 'constants.f'
       real(dp), intent(in) :: z
       real(dp) :: p0qg
       p0qg = ((1.0_dp-z)**2+z**2)*TR
       return
       end


!---
       subroutine xp0qiqi(z,p0qiqi)
       implicit none
       include 'types.f'
       include 'constants.f'
       real(dp), intent(in) :: z
       real(dp) :: p0qiqi(-1:0)
       p0qiqi(-1) = 3.0_dp/2.0_dp*CF !delta function
       p0qiqi(0) = (1.0_dp+z**2)*CF  !coeff of 1/(1-z)
       return
       end




c----- used at NNLO
c----- used at NNLO
c----- used at NNLO
c----- used at NNLO
c----- used at NNLO
c----- used at NNLO
c----- used at NNLO
c----- used at NNLO
c----- used at NNLO
c----- used at NNLO
c----- used at NNLO
c----- used at NNLO
c----- used at NNLO
c----- used at NNLO
c----- used at NNLO
c----- used at NNLO

      subroutine xp1gg(z,p1gg)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      real(dp), intent(in) :: z
      real(dp) :: p1gg(-1:1)
      real(dp) :: pggA(-1:1)
      real(dp) :: pggF(-1:1)
      real(dp) :: fpgg,fpggmz,s2z,splits2

      fpgg = 2*(1-z+z**2)**2/z/(1.0_dp-z)
      fpggmz = -2*(1+z+z**2)**2/z/(1.0_dp+z)
      s2z = splits2(z)

      pggA(0) = Ga1/8.0_dp*2*(1-z+z**2)**2/z
      pggA(-1) = CA*(-1+3*zeta3)+be0
      pggA(1) = CA*(fpgg*(-2*log(1-z)+log(z)/2.0_dp)*log(z)
     &             + fpggmz*(s2z+log(z)**2/2.0_dp)
     &             + 4*(1+z)*log(z)**2
     &             - 4/3.0_dp*(9+11*z**2)*log(z)
     &             -277/18.0_dp/z + 19*(1-z) + 277/18.0_dp*z**2 )
     &          + be0*(13/6.0_dp/z-3/2.0_dp*(1-z)
     &                -13/6.0_dp*z**2+(1+z)*log(z))

      pggF(0) = 0.0_dp
      pggF(-1) = -CF
      pggF(1) = CF*(4/3.0_dp/z-16+8*z+20/3.0_dp*z**2
     &              -2*(1+z)*log(z)**2 - 2*(3+5*z)*log(z) )


      p1gg(:) = CA*pggA(:) + TR*NF*pggF(:)

      return
      end

 !---
      function p1gqi(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      real(dp), intent(in) :: z
      real(dp) :: p1gqi
      real(dp) :: fpgq,fpgqmz,s2z,splits2

      fpgq = (1+(1-z)**2)/z
      fpgqmz = (1+(1+z)**2)/(-z)

      s2z = splits2(z)

      p1gqi =
     & CA*(fpgq*(log(1-z)**2-2*log(1-z)*log(z)-101/18.0_dp-pisq/6.0_dp)
     &  +fpgqmz*s2z + 2*z*log(1-z)+(2+z)*log(z)**2
     &  -(36+15*z+8*z**2)/3.0_dp*log(z)
     &  +(56-z+88*z**2)/18.0_dp )
     &  -CF*(fpgq*log(1-z)**2+(3*fpgq+2*z)*log(1-z)
     &   +(2-z)/2.0_dp*log(z)**2
     &   -(4+7*z)/2.0_dp*log(z) + (5+7*z)/2.0_dp )
     &  +be0*(fpgq*(log(1-z)+5/3.0_dp)+z)

      p1gqi = p1gqi*CF

      return
      end


 !---
      function p1qig(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp), intent(in) :: z
      real(dp) :: p1qig
      real(dp) :: fpqg,fpqgmz,s2z,splits2

      fpqg = (1-z)**2+z**2
      fpqgmz = (1+z)**2+z**2

      s2z = splits2(z)

      p1qig = CF*(fpqg*(log((1-z)/z)**2-2*log((1-z)/z)-pisq/3.0_dp+5)
     &             + 2*log(1-z) - (1-2*z)/2.0_dp*log(z)**2
     &             - (1-4*z)/2.0_dp*log(z) + 2 - 9/2.0_dp*z )
     &        +CA*(fpqg*(-log(1-z)**2+2*log(1-z)+22/3.0_dp*log(z)
     &             -109/9.0_dp+pisq/6.0_dp) +fpqgmz*s2z
     &             -2*log(1-z)-(1+2*z)*log(z)**2
     &             +(68*z-19)/3.0_dp*log(z)
     &             +20/9.0_dp/z+91/9.0_dp+7/9.0_dp*z )

      p1qig = p1qig*TR

      return
      end

 !---
      subroutine xp1qiqi(z,p1qiqi)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp),intent(in) :: z
      real(dp) :: p1qiqi(-1:1)
      real(dp) :: fp1qqV(-1:1),p1qqS

      call xp1qqV(z,fp1qqV)

      p1qiqi(-1) = fp1qqV(-1)
      p1qiqi(0) = fp1qqV(0)
      p1qiqi(1) = fp1qqV(1) + p1qqS(z)

      p1qiqi(:) = p1qiqi(:)*CF

      return
      end

 !---
      function p1qiqbi(z)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp),intent(in) :: z
      real(dp) :: p1qiqbi,p1qqbV,p1qqS

      p1qiqbi = p1qqbV(z) + p1qqS(z)

      p1qiqbi = p1qiqbi*CF

      return
      end


 !---
      function p1qiqj(z)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp),intent(in) :: z
      real(dp) :: p1qiqj,p1qqS

      p1qiqj = p1qqS(z)*CF

      return
      end


 !---
      function p1qiqbj(z)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp),intent(in) :: z
      real(dp) :: p1qiqbj,p1qqS

      p1qiqbj = p1qqS(z)*CF

      return
      end function p1qiqbj



 !---
      subroutine xp1qqV(z,p1qqV)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'

      real(dp), intent(in) :: z
      real(dp) :: p1qqV(-1:1)

      p1qqV(0) = Ga1/8.0_dp*(1+z**2)
      p1qqV(-1) = CF*(3/8.0_dp-pisq/2.0_dp+6*zeta3)
     &           +CA*(1/4.0_dp-3*zeta3)
     &           +be0*(1/8.0_dp+pisq/6.0_dp)

      p1qqV(1) = -CF*((1+z**2)/(1-z)*(2*log(1-z)+3/2.0_dp)*log(z)
     &                 +(1+z)/2.0_dp*log(z)**2
     &                 +(3+7*z)/2.0_dp*log(z)+5*(1-z) )
     &            +CA*(1/2.0_dp*(1+z**2)/(1-z)*log(z)**2
     &                 +(1+z)*log(z) + 3*(1-z) )
     &            +be0*(1/2.0_dp*(1+z**2)/(1-z)*log(z)+1-z)
       return
       end

 !---
      function p1qqbV(z)
      implicit none
      include 'types.f'
      include'constants.f'
      real(dp), intent(in) :: z
      real(dp) :: p1qqbV

      real(dp) :: s2z,splits2

      s2z = splits2(z)

      p1qqbV = (2*CF-CA)*((1+z**2)/(1+z)*(s2z+log(z)**2/2.0_dp)
     &                    +(1+z)*log(z) + 2*(1-z)  )

      end

 !---
      function p1qqS(z)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp),intent(in) :: z
      real(dp) :: p1qqS

      p1qqS = TR*(-(1+z)*log(z)**2+(1+5*z+8/3.0_dp*z**2)*log(z)
     &            +20/9.0_dp/z-2+6*z-56/9.0_dp*z**2 )

      return
      end

 !---
      subroutine xpqqpqq(z,pqqpqq)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp), intent(in) :: z
      real(dp) :: pqqpqq(-1:2)

      pqqpqq(2) = -2*(1-z)+(1+z-2*(1+z**2)/(1-z))*log(z)
      pqqpqq(1) = 4*(1+z**2)
      pqqpqq(0) = 3*(1+z**2)
      pqqpqq(-1) = -(9/4.0_dp+2*pisq/3.0_dp) + 3*3/2.0_dp

      pqqpqq(:)=pqqpqq(:)*CF**2

      return
      end

 !---
      function pqgpgq(z)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp),intent(in) :: z
      real(dp) :: pqgpgq

      pqgpgq = 2*(1+z)*log(z) + 4/3.0_dp/z + 1 - z - 4/3.0_dp*z**2

      pqgpgq=pqgpgq*CF*TR

      return
      end

 !---
      function pqgpgg(z)
! This is an implementation of the full convolution of
! P(0)qg x P(0)gg (z), including all color factors and
! endpoint contributions.
! It is obtained by dressing Eq. (A.11) of 1405.1044 with a color
! factor CA*TR and adding the endpoint explicit in Eq. (A.5)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'

      real(dp), intent(in) :: z
      real(dp) :: pqgpgg,p0qg

      pqgpgg = 2*((1-z)**2+z**2)*log(1-z)+2*(1+4*z)*log(z)
     &        +4/3.0_dp/z+1+8*z-31/3.0_dp*z**2
     &        +be0/two*p0qg(z)/(CA*TR)

      pqgpgg=pqgpgg*CA*TR

      return
      end

 !---
      function pqqpqg(z)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp), intent(in) :: z
      real(dp) :: pqqpqg

      pqqpqg =2*((1-z)**2+z**2)*log((1-z)/z)+(1-2*z)*log(z)-1/2.0_dp+2*z

      pqqpqg=pqqpqg*CF*TR

      return
      end


 !---
      subroutine xpggpgg(z,pggpgg)
! This is an implementation of the full convolution of
! P(0)gg x P(0)gg (z), including all color factors and
! endpoint contributions.
! It is obtained by dressing Eq. (A.11) of 1405.1044 with a color
! factor CA^2 and adding the endpoint explicit in Eq. (A.5)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'

      real(dp), intent(in) :: z
      real(dp) :: pggpgg(-1:2)
      real(dp) :: p0ggz(-1:0)

      call xp0gg(z,p0ggz)
      p0ggz(:) = p0ggz(:)*be0/2.0_dp*2.0_dp/CA**2

      pggpgg(1) = 8*(1-z+z**2)**2/z
      pggpgg(0) = 0.0_dp + p0ggz(0)
      pggpgg(-1) = -2*pisq/3.0_dp + be0**2/4.0_dp/CA**2
      pggpgg(2) = - 2*(2*(1-z+z**2)**2/z/(1-z)+4*(1+z))*log(z)
     &            - 44/3.0_dp/z + 12*(1-z)+44/3.0_dp*z**2

      pggpgg(:)=pggpgg(:)*CA**2

      return
      end

 !---
      function pgqpqg(z)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp), intent(in) :: z
      real(dp) :: pgqpqg

      pgqpqg = 2*(1+z)*log(z)+4/3.0_dp/z+1-z-4/3.0_dp*z**2

      pgqpqg=pgqpqg*CF*TR

      return
      end

 !---
      function pgqpqq(z)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp), intent(in) :: z
      real(dp) :: pgqpqq

      pgqpqq = 2*(1+(1-z)**2)/z*log(1-z)+(2-z)*log(z)+2-z/2.0_dp

      pgqpqq=pgqpqq*CF**2

      return
      end


 !---
      function pggpgq(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'

      real(dp), intent(in) :: z
      real(dp) :: pggpgq,p0gq
! This is an implementation of the full convolution of
! P(0)gg x P(0)gq (z), including all color factors and
! endpoint contributions.
! It is obtained by dressing Eq. (A.11) of 1405.1044 with a color
! factor CA*CF and adding the endpoint explicit in Eq. (A.5)

      pggpgq = 2*(1+(1-z)**2)/z*log((1-z)/z)-2*(4+z)*log(z)
     &        -31/3.0_dp/z+8+z+4/3.0_dp*z**2
     &         +be0/two*p0gq(z)/(CA*CF)

      pggpgq=pggpgq*CA*CF

      return
      end


 !---
      subroutine xIqqpqq(z,Iqqpqq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      real(dp), intent(in) :: z
      real(dp) :: Iqqpqq(-1:3)
      real(dp) :: fLI2z,Li2_TR
 !     complex(8) :: DiLog,zcom

 !     zcom = z + ci*0.0_dp
 !     fLI2z = real(dilog(zcom))
      fLI2z = Li2_TR(z)

      Iqqpqq(2) = 3*(1+z**2)
      Iqqpqq(1) = 1.5_dp*(1+z**2)
      Iqqpqq(0) = -pisq/2.0_dp*(1+z**2)
      Iqqpqq(-1) = 4*zeta3 + 1.5_dp*(-pisq/6.0_dp)
      Iqqpqq(3) = (log(z)**2-4*log(1-z)*log(z))*(1+z**2)/(1-z)
     &             - (1+z)*(fLI2z+log(z)**2/2.0_dp-pisq/6.0_dp)
     &             - z*log(z)-1+z
     &             + 1.5_dp*(1-z-(1+z**2)/(1-z)*log(z))

      Iqqpqq(:)=Iqqpqq(:)*CF**2

      return
      end

 !---
      function Iqgpgq(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp), intent(in) :: z
      real(dp) :: Iqgpgq
      real(dp) :: fLI2z,Li2_TR
 !     complex(8) :: DiLog,zcom

 !     zcom = z + ci*0.0_dp
 !     fLI2z = real(dilog(zcom))
      fLI2z = Li2_TR(z)

      Iqgpgq = -2*(1+z)*(fLI2z+log(z)**2/2.0_dp-pisq/6.0_dp)
     &         -(1+z+4/3.0_dp*z**2)*log((1-z)/z)
     &         +(4/3.0_dp/z+2)*log(1-z)
     &         +(2/z-5-z+4*z**2)/3.0_dp

      Iqgpgq=Iqgpgq*CF*TR

      return
      end


 !---
      function Iqgpgg(z)
! This is an implementation of the full convolution of
! I(1)qg x P(0)gg (z), including all color factors and
! endpoint contributions.
! It is obtained by dressing Eq. (A.12) of 1405.1044 with a color
! factor CA*TR and adding the endpoint explicit in Eq. (A.5)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      real(dp), intent(in) :: z
      real(dp) :: Iqgpgg
      real(dp) :: fLI2z,Li2_TR,I1qig
 !     complex(8) :: DiLog,zcom

 !     zcom = z + ci*0.0_dp
 !     fLI2z = real(dilog(zcom))
      fLI2z = Li2_TR(z)

      Iqgpgg = 2*((1-z)**2+z**2)*(log(1-z)*log((1-z)/z)-pisq/6.0_dp)
     &         -2*(1+4*z)*(fLI2z+log(z)**2/2.0_dp-pisq/6.0_dp)
     &         +(4/3.0_dp/z+1+12*z-43/3.0_dp*z**2)*log(1-z)
     &         +(1-8*z+31/3.0_dp*z**2)*log(z)
     &         +2/3.0_dp/z-13/6.0_dp-37/3.0_dp*z+83/6.0_dp*z**2
     &         +be0/two*I1qig(z)/(CA*TR)

      Iqgpgg=Iqgpgg*CA*TR

      return
      end


 !---
      function Iqqpqg(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp), intent(in) :: z
      real(dp) :: Iqqpqg
      real(dp) :: fLI2z,Li2_TR
 !     complex(8) :: DiLog,zcom

 !     zcom = z + ci*0.0_dp
 !     fLI2z = real(dilog(zcom))
      fLI2z = Li2_TR(z)

      Iqqpqg =  ((1-z)**2+z**2)*(log((1-z)/z)**2-pisq/6.0_dp)
     &        - (1-2*z)*(fLI2z+log(z)**2/2.0_dp-pisq/6.0_dp)
     &        + z*(7-3*z)*log((1-z)/z)
     &        - 2*(1+z)*log(1-z)
     &        - 1/2.0_dp-4*z+9/2.0_dp*z**2

      Iqqpqg=Iqqpqg*CF*TR

      return
      end

 !---
      subroutine xIggpgg(z,Iggpgg)
! This is an implementation of the full convolution of
! I(1)gg x P(0)gg (z), including all color factors and
! endpoint contributions.
! It is obtained by dressing Eq. (A.12) of 1405.1044 with a color
! factor CA^2 and adding the endpoint explicit in Eq. (A.5)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      real(dp), intent(in) :: z
      real(dp) :: Iggpgg(-1:3)
      real(dp) :: I1ggz(0:2)
      real(dp) :: fLI2z,fpgg,Li2_TR
 !     complex(8) :: DiLog,zcom

 !     zcom = z + ci*0.0_dp
 !     fLI2z = real(dilog(zcom))
      fLI2z = Li2_TR(z)

      fpgg = 2*(1-z+z**2)**2/z/(1-z)
      call xI1gg(z,I1ggz)
      I1ggz(:) = I1ggz(:)*be0/2.0_dp/CA**2

      Iggpgg(2) = 6*(1-z+z**2)**2/z
      Iggpgg(1) = 0.0_dp + I1ggz(1)
      Iggpgg(0) = -pisq/2.0_dp*2*(1-z+z**2)**2/z
      Iggpgg(-1) = 4*zeta3 + I1ggz(0)
      Iggpgg(3) = fpgg*(log(z)**2-4*log(1-z)*log(z))
     &            + 8*(1+z)*(fLI2z+log(z)**2/2.0_dp-pisq/6.0_dp)
     &            + (-22/3.0_dp/z+14-4*z+44/3.0_dp*z**2)*log((1-z)/z)
     &            - (22/3.0_dp/z+2+8*z)*log(1-z)
     &            + 67/9.0_dp/z-23/3.0_dp*(1-z)-67/9.0_dp*z**2
     &            + I1ggz(2)

      Iggpgg(:)=Iggpgg(:)*CA**2

      return
      end

 !---
      function Igqpqg(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp), intent(in) :: z
      real(dp) :: Igqpqg
      real(dp) :: fLI2z,Li2_TR
 !     complex(8) :: DiLog,zcom


 !     zcom = z + ci*0.0_dp
 !     fLI2z = real(dilog(zcom))
      fLI2z = Li2_TR(z)

      Igqpqg = -2*(1+z)*(fLI2z+log(z)**2/2.0_dp-pisq/6.0_dp)
     &          +(4/3.0_dp/z-3*z-4/3.0_dp*z**2)*log((1-z)/z)
     &          +(1+2*z)*log(1-z)
     &          -13/9.0_dp/z+4/3.0_dp+2/3.0_dp*z-5/9.0_dp*z**2

      Igqpqg=Igqpqg*CF*TR

      return
      end

 !---
      function Igqpqq(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp), intent(in) :: z
      real(dp) :: Igqpqq
      real(dp) :: fLI2z,fpgq,Li2_TR
 !     complex(8) :: DiLog,zcom

 !     zcom = z + ci*0.0_dp
 !     fLI2z = real(dilog(zcom))
      fLI2z = Li2_TR(z)

      fpgq = (1+(1-z)**2)/z

      Igqpqq = 2*fpgq*(log(1-z)*log((1-z)/z)-pisq/6.0_dp+5/8.0_dp)
     &        -(2-z)*(fLI2z+log(z)**2/2.0_dp-pisq/6.0_dp-1/4.0_dp)
     &        +(4+3*z)/2.0_dp*log(1-z) - (2+z)/2.0_dp*log(z)

      Igqpqq=Igqpqq*CF**2

      return
      end


 !---
      function Iggpgq(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp), intent(in) :: z
      real(dp) :: Iggpgq
      real(dp) :: fLI2z,fpgq,Li2_TR
 !     complex(8) :: DiLog, zcom

 !     zcom = z + ci*0.0_dp
 !     fLI2z = real(dilog(zcom))
      fLI2z = Li2_TR(z)

      fpgq = (1.0_dp+(1.0_dp-z)**2)/z

      Iggpgq = fpgq*(log((1-z)/z)**2-pisq/6.0_dp)
     &          +2*(4+z)*(fLI2z+log(z)**2/2.0_dp-pisq/6.0_dp)
     &          +(21-26*z+5*z**2)/6.0_dp/z
     &          +(-3/z+10+3*z+4/3.0_dp*z**2)*log((1-z)/z)
     &          -(22/3.0_dp/z+2+2*z)*log(1.0_dp-z)

      Iggpgq=Iggpgq*CA*CF

      return
      end
