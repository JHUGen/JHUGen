      subroutine box3(k1,k2,k3,k4,k5,k6,za,zb,app,apm,amp,amm)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'Higgsint.f'
      real(dp):: mtsq,s12,s34,s56,s134,s156
      integer:: k1,k2,k3,k4,k5,k6
      complex(dp):: app,apm,amp,amm,boxamp,zab2
c--- statement functions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
c--- end statement functions

      apm=czip
      amp=czip

      mtsq=mt**2
      s12=s(k1,k2)
      s34=s(k3,k4)
      s56=s(k5,k6)
      s134=s(k1,k3)+s(k1,k4)+s(k3,k4)
      s156=s(k1,k5)+s(k1,k6)+s(k5,k6)

      boxamp=-(s134*s156-s56*s34+mtsq*s12)*za(k2,k3)/s12**2/za(k1,k2)**2
     & *(za(k1,k5)*zb(k5,k6)*zb(k1,k2)
     & +zb(k1,k6)*mtsq*s12/zab2(k2,k5,k6,k1))
     & *(za(k1,k5)*za(k2,k3)*zb(k3,k4)*zb(k1,k2)
     &  -za(k2,k5)*zb(k1,k4)*mtsq*s12/zab2(k2,k5,k6,k1))
      app=boxamp*(-half*im)/(s34*s56)
c      write(6,*)
c      write(6,*) 'box3:app',app
      boxamp=-(s134*s156-s56*s34+mtsq*s12)*zb(k2,k6)/s12**2/zb(k1,k2)**2
     &  *(za(k3,k4)*za(k1,k2)*zb(k1,k4)
     &  +za(k1,k3)*mtsq*s12/zab2(k1,k5,k6,k2))
     & *(za(k5,k6)*za(k1,k2)*zb(k2,k6)*zb(k1,k4)
     & -za(k1,k5)*zb(k2,k4)*mtsq*s12/zab2(k1,k5,k6,k2))
      amm=boxamp*(-half*im)/(s34*s56)
c      write(6,*) 'box3:amm',amm

      if (Higgsint) return

      boxamp=-(s134*s156-s56*s34+mtsq*s12)*zb(k2,k6)
     & *mtsq/zab2(k1,k5,k6,k2)**2/s12**2
     & *(za(k1,k2)*za(k4,k3)*zb(k1,k4)
     & -za(k1,k3)*mtsq*s12/zab2(k1,k5,k6,k2))
     & *(za(k1,k2)*zb(k1,k4)*zab2(k5,k3,k4,k2)
     &  -za(k1,k5)*zb(k2,k4)*mtsq*s12/zab2(k1,k5,k6,k2))
      apm=boxamp*(-half*im)/(s34*s56)
c      write(6,*) 'box3:apm',apm
      boxamp=-(s134*s156-s56*s34+mtsq*s12)*za(k2,k3)
     & *mtsq/zab2(k2,k5,k6,k1)**2/s12**2
     & *(za(k1,k5)*zb(k1,k2)*zb(k5,k6)
     &  +zb(k1,k6)*mtsq*s12/zab2(k2,k5,k6,k1))
     & *(za(k1,k5)*zb(k1,k2)*zab2(k2,k5,k6,k4)
     & +za(k2,k5)*zb(k1,k4)*mtsq*s12/zab2(k2,k5,k6,k1))
      amp=boxamp*(-half*im)/(s34*s56)
c      write(6,*) 'box3:amp',amp


      return
      end
