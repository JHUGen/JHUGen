      double complex function D2six(k1,k2,k3,k4,k5,k6)
      implicit none
C----- six dimensional box corresponding to 
C----- qlI4(s56,0d0,0d0,s34,s156,s12,mtsq,0d0,0d0,0d0,musq,e)
C----- multiplied by a factor of -C(0)/2
      include 'constants.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'scale.f'
      double complex qlI4,qlI3,IntC(4),IntD
      double precision s12,s34,s56,s134,s156,C(0:4),mtsq
      integer e,k1,k2,k3,k4,k5,k6
      mtsq=mt**2
      s12=s(k1,k2)
      s34=s(k3,k4)
      s56=s(k5,k6)
      s134=s(k1,k3)+s(k1,k4)+s(k3,k4)
      s156=s(k1,k5)+s(k1,k6)+s(k5,k6)

      C(1)=2d0/(mtsq-s156)
      C(2)=2d0*(s156-s34)/(s12*(mtsq-s156))
      C(3)=2d0*((s156+mtsq)*(s56+s34-s12)-2d0*s156*mtsq-2d0*s34*s56)
     & /(s12*(mtsq-s156)**2)
      C(4)=2d0*(s156-s56)/(s12*(mtsq-s156))
      C(0)=4d0*(s134*s156-s34*s56)/(s12*(mtsq-s156)**2)

      e=0
      IntC(1)=qlI3(0d0,0d0,s12,0d0,0d0,0d0,musq,e)   !C1
      IntC(2)=qlI3(0d0,s34,s156,0d0,0d0,mtsq,musq,e) !C5
      IntC(3)=qlI3(s12,s56,s34,0d0,0d0,mtsq,musq,e)  !C6
      IntC(4)=qlI3(0d0,s156,s56,0d0,0d0,mtsq,musq,e) !C4
      IntD=qlI4(s56,0d0,0d0,s34,s156,s12,mtsq,0d0,0d0,0d0,musq,e)
      D2six=0.5d0*(C(1)*IntC(1)+C(2)*IntC(2)+C(3)*IntC(3)+C(4)*IntC(4)
     & +2d0*IntD)

c---debug
C---real sixdim box
c      D2six=-(C(1)*IntC(1)+C(2)*IntC(2)+C(3)*IntC(3)+C(4)*IntC(4)
c     & +2d0*IntD)/C(0)
c---debug
      return
      end
