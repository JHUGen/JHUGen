      function aqqb_wbb(i1,i2,i5,i6,i3,i4)
      implicit none
      include 'types.f'
      complex(dp):: aqqb_wbb

C---This is the amplitude for
C---q_L(p1)+q_L(p6) --> q_L(p2)+q_L(p5)+W(l(p3)+antilepton(p4))
c---with no couplings included
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
C---'prods.f' includes both the declaration and the common for s and za,zb
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: t2
      real(dp):: s234,s256,prop
      s234=s(i2,i3)+s(i2,i4)+s(i3,i4)
      s256=s(i2,i6)+s(i2,i5)+s(i5,i6)
      prop=s(i5,i6)*sqrt((s(i3,i4)-wmass**2)**2+(wmass*wwidth)**2)
      aqqb_wbb=
     & +za(i3,i2)*zb(i6,i1)*t2(i4,i2,i3,i5)/(prop*s234)
     & +za(i5,i2)*zb(i4,i1)*t2(i6,i2,i5,i3)/(prop*s256)
      return
      end
