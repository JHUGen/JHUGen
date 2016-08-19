      subroutine Acalc(k1,k2,k3,k4,k5,k6,mt,A)
      implicit none 
      include 'constants.f'
      include 'sprods_com.f'
      include 'ZZclabels.f'
      include 'ZZdlabels.f'
      include 'ggZZintegrals.f'
      double complex A(6)
      double precision mt,mtsq,s12,s134,s234,s134bar,s234bar,p3sq,Y
      integer k1,k2,k3,k4,k5,k6
      mtsq=mt**2
      s12=s(k1,k2)
      s134bar=s(k1,k3)+s(k1,k4)
      s234bar=s(k2,k3)+s(k2,k4)
      p3sq=s(k3,k4)
      s234=s234bar+p3sq
      s134=s134bar+p3sq
      Y=s134bar*s234bar-s12*p3sq

      A(1)=+0.5d0*mtsq/s12*(
     & +2d0*C0(c1_34)*s134bar
     & +2d0*C0(c2_56)*(-s234bar-s12)
     & +2d0*C0(c2_34)*s234bar
     & +2d0*C0(c1_56)*(-s12-s134bar)
     & -2d0*Y*D0(d1_34_2)
     & +s12*(s12-4*mtsq)*(D0(d1_2_34)+D0(d2_1_34)+D0(d1_34_2)))

      A(2)=2d0*mtsq*(D0(d6_1_2_34)+D0(d6_2_1_34)+C0(c12_34)
     & +mtsq*(D0(d1_2_34)+D0(d2_1_34)-D0(d1_34_2)))


c      A(2)=mtsq/Y*(
c     & +C0(c1_34)*s134*s134bar
c     & +C0(c2_56)*s134*(-s234bar-s12)
c     & +C0(c2_34)*s234*s234bar
c     & +C0(c1_56)*s234*(-s134bar-s12)
c     & +C0(c1_2)*s12*(s234+s134)
c     & +C0(c12_34)*(2d0*Y-((s134bar+s234bar)**2-4d0*s12*p3sq))
c     & -D0(d2_1_34)*s12*s134**2-D0(d1_2_34)*s12*s234**2
c     & -2d0*Y*mtsq*(D0(d1_2_34)+D0(d1_34_2)+D0(d2_1_34)))

      A(3)=0.5d0*mtsq*s12*(D0(d1_2_34)-D0(d2_1_34)-D0(d1_34_2))
      A(4)=0.5d0*mtsq*s12*(D0(d2_1_34)-D0(d1_2_34)-D0(d1_34_2))

      A(5)=0.5d0*mtsq*s12/(s134*s234)*(
     & +2d0*s134*D0(d6_1_2_34)+2d0*s234*D0(d6_2_1_34)
     & -s134*s234*D0(d1_34_2)
     & +4d0*mtsq*(s134*D0(d1_2_34)+s234*D0(d2_1_34))
     & +2d0*(s234+s134)*C0(c12_34))

c      A(5)=0.5d0*mtsq*s12/Y*(
c     & +C0(c1_34)*2d0*s134bar
c     & -C0(c1_56)*2d0*(s12+s134bar)
c     & +C0(c1_2)*2d0*s12
c     & -D0(d1_2_34)*s12*s234-D0(d2_1_34)*s12*s134-Y*D0(d1_34_2))


      A(6)=czip
c      A(6)=mtsq*s12/Y*(
c     & -C0(c1_34)*s134bar
c     & +C0(c1_56)*(s12+s134bar)
c     & +C0(c2_34)*s234bar
c     & -C0(c2_56)*(s12+s234bar))
c      write(6,*) 'A(6)',A(6)

      return
      end
