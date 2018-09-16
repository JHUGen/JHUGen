      subroutine Acalc(k1,k2,k3,k4,k5,k6,mt,A)
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'ZZclabels.f'
      include 'ZZdlabels.f'
      include 'ggZZintegrals.f'
      complex(dp):: A(6)
      real(dp):: mt,mtsq,s12,s134,s234,s134bar,s234bar,p3sq,Y
      integer:: k1,k2,k3,k4,k5,k6
      mtsq=mt**2
      s12=s(k1,k2)
      s134bar=s(k1,k3)+s(k1,k4)
      s234bar=s(k2,k3)+s(k2,k4)
      p3sq=s(k3,k4)
      s234=s234bar+p3sq
      s134=s134bar+p3sq
      Y=s134bar*s234bar-s12*p3sq

      A(1)=+half*mtsq/s12*(
     & +two*C0(c1_34)*s134bar
     & +two*C0(c2_56)*(-s234bar-s12)
     & +two*C0(c2_34)*s234bar
     & +two*C0(c1_56)*(-s12-s134bar)
     & -two*Y*D0(d1_34_2)
     & +s12*(s12-4*mtsq)*(D0(d1_2_34)+D0(d2_1_34)+D0(d1_34_2)))

      A(2)=two*mtsq*(D0(d6_1_2_34)+D0(d6_2_1_34)+C0(c12_34)
     & +mtsq*(D0(d1_2_34)+D0(d2_1_34)-D0(d1_34_2)))


c      A(2)=mtsq/Y*(
c     & +C0(c1_34)*s134*s134bar
c     & +C0(c2_56)*s134*(-s234bar-s12)
c     & +C0(c2_34)*s234*s234bar
c     & +C0(c1_56)*s234*(-s134bar-s12)
c     & +C0(c1_2)*s12*(s234+s134)
c     & +C0(c12_34)*(two*Y-((s134bar+s234bar)**2-four*s12*p3sq))
c     & -D0(d2_1_34)*s12*s134**2-D0(d1_2_34)*s12*s234**2
c     & -two*Y*mtsq*(D0(d1_2_34)+D0(d1_34_2)+D0(d2_1_34)))

      A(3)=half*mtsq*s12*(D0(d1_2_34)-D0(d2_1_34)-D0(d1_34_2))
      A(4)=half*mtsq*s12*(D0(d2_1_34)-D0(d1_2_34)-D0(d1_34_2))

      A(5)=half*mtsq*s12/(s134*s234)*(
     & +two*s134*D0(d6_1_2_34)+two*s234*D0(d6_2_1_34)
     & -s134*s234*D0(d1_34_2)
     & +four*mtsq*(s134*D0(d1_2_34)+s234*D0(d2_1_34))
     & +two*(s234+s134)*C0(c12_34))

c      A(5)=half*mtsq*s12/Y*(
c     & +C0(c1_34)*two*s134bar
c     & -C0(c1_56)*two*(s12+s134bar)
c     & +C0(c1_2)*two*s12
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
