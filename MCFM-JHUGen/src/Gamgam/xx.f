      include 'types.f'
      AGTYAs=
     & (128._dp*Li4z-128._dp*Li4x+128._dp*Li4y+
     & (-64._dp/3._dp+128._dp*Ly)*Li3x
     & +(64._dp/
     & 3._dp*Lx-64._dp/3._dp*pisq)*Li2x+16._dp/3._dp*Lx**4
     & -64._dp/3._dp*Lx**3._dp*Ly+
     & (-16+32._dp/3._dp*pisq+32._dp*Ly**2)*Lx**2
     & +
     & (-64._dp/3._dp*pisq*Ly+48+160._dp/9._dp*pisq)*Lx
     & +64._dp/3._dp*zeta3+224._dp/45._dp*pisq**2-128._dp*Ly*zeta3
     & )*t/u
     & +(32._dp/3._dp*Li3x-32._dp/3._dp*Li3y+
     & (-32._dp/3._dp*Lx-32._dp/3._dp*Ly)*Li2x
     & +(-32._dp/3._dp*Ly**2-80._dp/9._dp*pisq-64._dp/3._dp
     & )*Lx+(64._dp/3._dp+32._dp/3._dp*pisq)*Ly
     & )*t**2/s**2+24._dp*Lx**2._dp*t**2/u**2
     & +(416._dp/3._dp*Li3x+64._dp*Li3y*Lx-416._dp/
     & 3._dp*Li2x*Lx+(8._dp*Ly**2+16)*Lx**2
     & +(-8._dp/3._dp*Ly+80._dp/3._dp+112._dp/
     & 9._dp*pisq-64._dp*zeta3-64._dp*Ly**2)*Lx-416._dp/
     & 3._dp*zeta3-148._dp/9._dp*pisq+44._dp/45._dp*pisq**2)

!     +t leftrightarrow u

      AGTYBs=(-112._dp*Li4z-88._dp*Li4y+
     & (-128._dp*Ly+48._dp*Lx-64)*Li3x
     & +(-16._dp*Ly-16._dp*Lx+12)*Li3y+
     & (12._dp*Ly-4._dp*Ly**2+8._dp*Lx**2-8._dp*pisq+64._dp*Lx)*Li2x
     & +2._dp/3._dp*Lx**4+56._dp/3._dp*Lx**3._dp*Ly+
     & (44._dp*Ly-4._dp*pisq+2-32._dp*Ly**2)*Lx**2
     & +(-4._dp*Ly**3-8._dp-32._dp*zeta3-80._dp/
     & 3._dp*pisq+6._dp*Ly**2+56._dp/3._dp*pisq*Ly)*Lx
     & +Ly**4+6._dp*Ly**3+(-10._dp/3._dp*pisq-5._dp)*Ly**2
     & +(-39-18._dp*pisq+144._dp*zeta3)*Ly
     & +3._dp*Ls+187._dp/4._dp-4._dp*pisq*Ls+4._dp/
     & 45._dp*pisq**2-5._dp*pisq-20._dp*zeta3+48._dp*zeta3*Ls)*t/u
     & +(-12._dp*Lx**2+(24._dp*Ly+24._dp
     & )*Lx-12._dp*Ly**2-24._dp*Ly-12._dp*pisq
     & )*t**2/s**2+8._dp*Lx**2._dp*t**2/u**2
     & +(-80*Li4y+32._dp*Lx*Li3x+
     & (-128._dp*Lx-152
     & )*Li3y+152._dp*Li2x*Lx
     & +8._dp*Ly**2._dp*Li2y+(-16._dp*Ly**2-24)*Lx**2+
     & (60*Ly**2+(28+32._dp/3._dp*pisq
     & )*Ly-58)*Lx
     & +14._dp/3._dp*Ly**4+44._dp/3._dp*Ly**3+8._dp/
     & 3._dp*Ly**2*pisq+(96._dp*zeta3-32._dp/3._dp*pisq
     & )*Ly+32._dp/45._dp*pisq**2+16._dp*zeta3-86._dp/3._dp*pisq-2._dp)

!     +t leftrightarrow u

      AGTYCs=(-20*Li4z+28._dp*Li4x+
     & (-28._dp*Ly-10._dp*Lx+1._dp/3._dp)*Li3x
     & +(6._dp*Lx**2+(-1._dp/
     & 3-4._dp*Ly)*Lx-2._dp*Ly**2+58._dp/3._dp*Ly+4._dp/
     & 3._dp*pisq)*Li2x
     & +(58._dp/
     & 3-12._dp*Lx-12._dp*Ly)*Li3y-1._dp/
     & 6._dp*Lx**4+(10._dp/3._dp*Ly+13._dp/9._dp)*Lx**3
     & +(-9._dp*Ly**2-1._dp/3._dp*pisq+4._dp/9._dp+11._dp/
     & 2._dp*Ls+13._dp*Ly)*Lx**2
     & +(-8._dp/3._dp*Ly**3+50._dp/3._dp*Ly**2+(17._dp/
     & 3._dp*pisq-28._dp/3._dp+11._dp*Ls)*Ly-563._dp/
     & 27-233._dp/36._dp*pisq+55._dp/3._dp*Ls)*Lx
     & +(-2._dp/3._dp*pisq+80._dp/9._dp)*Ly**2+
     & (-284._dp/ 27-299._dp/36._dp*pisq+26._dp*zeta3+55._dp/
     & 3._dp*Ls)*Ly-209._dp/36._dp*pisq*Ls
     & -2._dp*zeta3*Ls+121._dp/12._dp*Ls**2-13._dp*Ls-1142._dp/
     & 81._dp-197._dp/360*pisq**2+461._dp/36._dp*pisq-55._dp/
     & 18._dp*zeta3)*t/u
     & +(-5._dp/4._dp*Lx**2+(5._dp/ 2+5._dp/
     & 2._dp*Ly)*Lx-5._dp/4._dp*Ly**2-5._dp/ 2._dp*Ly-5._dp/
     & 4._dp*pisq)*t**2/s**2+1._dp/ 2._dp*Lx**2._dp*t**2/u**2
     & +(24._dp*Li4y-20*Lx*Li3x+
     & (-40*Lx-22
     & )*Li3y+22._dp*Li2x*Lx+8._dp*Ly**2._dp*Li2y

     & +4._dp/3._dp*Lx**3._dp*Ly+(-6._dp*Ly**2-575._dp/36
     & )*Lx**2+(46._dp/3._dp*Ly**2+(73._dp/
     & 12+4._dp*pisq)*Ly-637._dp/18._dp)*Lx
     & +1._dp/3._dp*Ly**4+59._dp/9._dp*Ly**3+(2._dp/
     & 3._dp*pisq+11._dp*Ls)*Ly**2+
     & (44._dp*zeta3-4._dp/9._dp*pisq+11._dp*Ls)*Ly
     & -38._dp/45._dp*pisq**2+77._dp/72._dp*pisq+2._dp*zeta3)
!     +t leftrightarrow u

      AGTYD1s=
     & (96._dp*Li4z-48._dp*Li4x+52._dp*Li4y+
     & (124._dp*Ly-8._dp*Lx+46)*Li3x
     & +
     & (-16._dp*Lx**2+(-46+8._dp*Ly
     & )*Lx+6._dp*Ly**2-30*Ly+4._dp/3._dp*pisq
     & )*Li2x
     & +(-30+28._dp*Ly+36._dp*Lx)*Li3y+1._dp/
     & 2._dp*Lx**4+(-56._dp/3._dp*Ly-100._dp/9._dp)*Lx**3

     & +(-125._dp/3._dp*Ly+39._dp*Ly**2+214._dp/
     & 9+3._dp*pisq-22._dp*Ls)*Lx**2
     & +(14._dp/3._dp*Ly**3-73._dp/3._dp*Ly**2+
     & (-24._dp*pisq+4._dp)*Ly+155._dp/9._dp*pisq+148._dp/
     & 3)*Lx
     & -74._dp/9._dp*Ly**3+(10._dp/3._dp*pisq-55._dp/
     & 9-11._dp*Ls)*Ly**2+(-33._dp*Ls+136._dp/9._dp*pisq-140*zeta3
     & +241._dp/3._dp)*Ly
     & -43417._dp/324._dp+23._dp/
     & 6._dp*pisq*Ls-52._dp*zeta3*Ls-173._dp/18._dp*pisq+227._dp/
     & 180*pisq**2+1834._dp/ 27._dp*Ls+515._dp/9._dp*zeta3
     & )*t/u
     & +(14._dp*Lx**2+(-28-28._dp*Ly
     & )*Lx+14._dp*Ly**2+28._dp*Ly+14._dp*pisq
     & )*t**2/s**2-5._dp*Lx**2._dp*t**2/u**2
     & +(-8._dp*Li4y+24._dp*Lx*Li3x+
     & (144._dp*Lx+120
     & )*Li3y-120*Li2x*Lx-20*Ly**2._dp*Li2y
     & -8._dp/3._dp*Lx**3._dp*Ly+(20*Ly**2+472._dp/9._dp
     & )*Lx**2+(-182._dp/3._dp*Ly**2+(-40._dp/
     & 3._dp*pisq-104._dp/3._dp)*Ly+898._dp/9._dp)*Lx
     & -3._dp*Ly**4-184._dp/9._dp*Ly**3+(-22._dp*Ls-8._dp/
     & 3._dp*pisq)*Ly**2+(-22._dp*Ls+56._dp/
     & 9._dp*pisq-136._dp*zeta3)*Ly
     & +4._dp/3._dp*pisq**2+148._dp/9._dp*pisq+2-12._dp*zeta3)
!     +t leftrightarrow u

      AGTYE1s=(22._dp/9._dp*Lx**3+(-76._dp/
     & 9+4._dp*Ls+2._dp/3._dp*Ly)*Lx**2+(1._dp/
     & 3._dp*Ly**2+Ly+16._dp/9._dp*pisq-7._dp/3._dp)*Lx
     & +11._dp/9._dp*Ly**3+(7._dp/9._dp+2._dp*Ls
     & )*Ly**2+(6._dp*Ls+8._dp/9._dp*pisq-37._dp/3._dp
     & )*Ly
     & +19._dp/9._dp*pisq-328._dp/ 27._dp*Ls+3401._dp/162._dp-2._dp/
     & 9._dp*zeta3-1._dp/3._dp*pisq*Ls)*t/u
     & +(-46._dp/9._dp*Lx**2+(2._dp/3._dp*Ly+2._dp/
     & 3._dp*Ly**2-76._dp/9._dp)*Lx
     & +22._dp/9._dp*Ly**3+4._dp*Ly**2._dp*Ls+(16._dp/
     & 9._dp*pisq+4._dp*Ls)*Ly+8._dp/9._dp*pisq)
!     + t leftrightarrow u

      AGTYE2s=(-4._dp/3._dp*Li3x-4._dp/
     & 3._dp*Li3y+(-4._dp/3._dp*Ly+4._dp/3._dp*Lx
     & )*Li2x
     & -11._dp/18._dp*Lx**3+(-1._dp/ 2._dp*Ly+5._dp/
     & 6-Ls)*Lx**2
     & +(-5._dp/3._dp*Ly**2+(10._dp/9._dp-2._dp*Ls
     & )*Ly-1._dp/9._dp*pisq+37._dp/9._dp-31._dp/6._dp*Ls
     & )*Lx
     & -41._dp/18._dp*Ly**2+(5._dp/9._dp*pisq+43._dp/9._dp-31._dp/
     & 6._dp*Ls)*Ly+65._dp/81._dp+19._dp/
     & 18._dp*pisq*Ls-11._dp/3._dp*Ls**2
     & -13._dp/9._dp*zeta3+206._dp/
     & 27._dp*Ls-275._dp/108._dp*pisq)*t/u
     & +(-Lx**2+(2._dp*Ly+2
     & )*Lx-Ly**2-2._dp*Ly-pisq)*t**2/s**2
     & +(14._dp/9._dp*Lx**2+(2._dp/3._dp*Ly-1._dp/
     & 3._dp*Ly**2+38._dp/9._dp)*Lx-11._dp/
     & 9._dp*Ly**3-2._dp*Ly**2._dp*Ls
     & -17._dp/18._dp*pisq+(-8._dp/
     & 9._dp*pisq-2._dp*Ls)*Ly)
!     + t leftrightarrow u

      AGTYF1s=(5._dp/36._dp*Lx**2+(-10._dp/ 27+1._dp/
     & 3._dp*Ls+1._dp/18._dp*Ly)*Lx+5._dp/36._dp*Ly**2+
     & (-10._dp/ 27+1._dp/3._dp*Ls)*Ly
     & +1._dp/3._dp*Ls**2+1._dp/54._dp*pisq-20._dp/ 27._dp*Ls
     & )*t/u
!      + t leftrightarrow u

      AGTYD2s=
     & (48._dp*Li4z-16._dp*Li4x+24._dp*Li4y+
     & (56._dp*Ly-8._dp*Lx+20)*Li3x
     & +(8._dp*Lx-12+16._dp*Ly)*Li3y+
     & (16._dp/3._dp*pisq-20*Lx-12._dp*Ly-8._dp*Lx**2+4._dp*Ly**2
     & )*Li2x
     & +1._dp/3._dp*Lx**4+(-8._dp*Ly-70._dp/9._dp
     & )*Lx**3+(-4._dp*pisq+286._dp/
     & 9-16._dp*Ly+14._dp*Ly**2-44._dp/3._dp*Ls)*Lx**2
     & +(-22._dp/
     & 9._dp*pisq+4._dp*Ly**3-8._dp*pisq*Ly-6._dp*Ly**2
     & )*Lx-44._dp/9._dp*Ly**3+(-4._dp/3._dp*pisq+35._dp/
     & 9-22._dp/3._dp*Ls)*Ly**2
     & +(57-26._dp/9._dp*pisq-72._dp*zeta3-22._dp*Ls
     & )*Ly+479._dp/9._dp*zeta3+19._dp/
     & 60*pisq**2-52._dp*zeta3*Ls+1141._dp/ 27._dp*Ls-215._dp/
     & 18._dp*pisq
     & -43417._dp/324._dp+23._dp/6._dp*pisq*Ls)*t/u
     & +(6._dp*Lx**2+(-12-12._dp*Ly
     & )*Lx+6._dp*Ly**2+12._dp*Ly+6._dp*pisq
     & )*t**2/s**2-6._dp*Lx**2._dp*t**2/u**2
     & +
     & (16._dp*Li4y+48._dp*Li3x*Ly+64._dp*Li3y
     & -8._dp*Ly**2._dp*Li2y-64._dp*Li2x*Lx
     & -4._dp/3._dp*Lx**4+(-20._dp/3._dp*pisq+6._dp*Ly**2
     & )*Lx**2+(-24._dp*Ly**2+(-16._dp/
     & 3._dp*pisq-14)*Ly-148._dp/9._dp*pisq)*Lx
     & -112._dp/9._dp*Ly**3+(-44._dp/3._dp*Ls+298._dp/9._dp
     & )*Ly**2+(538._dp/9._dp-48._dp*zeta3-44._dp/3._dp*Ls
     & )*Ly-8._dp*zeta3-1._dp/3._dp*pisq**2+61._dp/9._dp*pisq)

!     +t leftrightarrow u

      AGTYE3s=(16._dp/9._dp*Lx**3+(-76._dp/9._dp+8._dp/
     & 3._dp*Ls)*Lx**2+16._dp/9._dp*pisq*Lx
     & +8._dp/9._dp*Ly**3+(4._dp/3._dp*Ls-2._dp/9._dp
     & )*Ly**2+(8._dp/9._dp*pisq+4._dp*Ls-10._dp
     & )*Ly
     & -1._dp/3._dp*pisq*Ls-202._dp/ 27._dp*Ls+19._dp/
     & 9._dp*pisq-2._dp/9._dp*zeta3+3401._dp/162._dp)*t/u
     & +(16._dp/9._dp*pisq*Lx+16._dp/9._dp*Ly**3+
     & (8._dp/3._dp*Ls-52._dp/9._dp)*Ly**2
     & +(-76._dp/9._dp+8._dp/3._dp*Ls)*Ly+8._dp/
     & 9._dp*pisq)
!     + t leftrightarrow u

      AGTYE4s=(16._dp/9._dp*Ly**3-4._dp/9._dp*Ly**2+16._dp/
     & 3._dp*Ly**2._dp*Ls-20*Ly+64._dp/
     & 9._dp*pisq*Ly+16._dp*Ly*Ls
     & +32._dp/9._dp*Lx**3-152._dp/
     & 9._dp*Lx**2+32._dp/3._dp*Lx**2._dp*Ls+128._dp/
     & 9._dp*pisq*Lx
     & -2._dp/3._dp*pisq*Ls+110._dp/
     & 9._dp*pisq+3401._dp/81._dp-908._dp/ 27._dp*Ls-4._dp/9._dp*zeta3
     & )*t/u
     & +(32._dp/9._dp*Lx**3-104._dp/9._dp*Lx**2+32._dp/
     & 3._dp*Lx**2._dp*Ls-152._dp/9._dp*Lx+32._dp/3._dp*Lx*Ls
     & +128._dp/9._dp*pisq*Lx+64._dp/9._dp*pisq)
!       + t leftrightarrow u

      AGTYF2s=(-92._dp/ 27._dp*pisq+32._dp/9._dp*Ls**2
     & -160._dp/ 27._dp*Ls)*t/u
!      + t leftrightarrow u




      AGTYAu=
     & (24._dp*pisq-48._dp*Lx*Ly+24._dp*Ly**2+24._dp*Lx**2
     & )*t**2/s**2 +24._dp*Ly**2._dp*s**2/t**2
     & +((64._dp*Ly+32._dp/3._dp)
     & *Li3x+64._dp*Li3y*Lx-32._dp/
     & 3._dp*Li2x*Lx+(-8._dp+16._dp*Ly**2)
     & *Lx**2
     & +((-64._dp/3._dp*pisq+16
     & )*Ly-16._dp/9._dp*pisq+24-64._dp*zeta3+32._dp/
     & 3._dp*Ly**3-32._dp/3._dp*Ly**2)*Lx+64._dp/9._dp*Ly**3
     & +(-48+64._dp*zeta3+32._dp/9._dp*pisq
     & )*Ly+(-16._dp/3._dp*pisq-16)*Ly**2-16._dp/
     & 3._dp*Ly**4-8._dp*pisq+32._dp/3._dp*zeta3
     & +88._dp/45._dp*pisq**2)*(t**2+s**2)/(s*t)
     & +
     & (-128._dp*Li4x-128._dp*Li4y+128._dp*Li4z+
     & (64._dp*Ly+32._dp/3._dp)*Li3x
     & +(128._dp*Ly+64._dp/3._dp-64._dp*Lx
     & )*Li3y+(-32._dp/3._dp*Lx+64._dp/3._dp*Ly+64._dp/
     & 3._dp*pisq)*Li2x
     & +16._dp/3._dp*Lx**4-64._dp/3._dp*Lx**3._dp*Ly+
     & (16._dp*Ly**2-8._dp+32._dp/3._dp*pisq)*Lx**2
     & +(32._dp/3._dp*Ly**2+32._dp/
     & 3._dp*Ly**3+64._dp*zeta3-16._dp/9._dp*pisq+24+16._dp*Ly
     & )*Lx+88._dp/15._dp*pisq**2-32._dp/3._dp*zeta3
     & +(-32._dp/9._dp*pisq-64._dp*zeta3
     & )*Ly+16._dp*Ly**2._dp*pisq-8._dp*pisq)*(t**2-s**2)/(s*t)
     & +(-32._dp/3._dp*Li3x-64._dp/3._dp*Li3y+
     & (32._dp/3._dp*Lx-64._dp/3._dp*Ly)*Li2x
     & +(-32._dp/3._dp*Ly**2-64._dp/3._dp+16._dp/9._dp*pisq
     & )*Lx+32._dp/9._dp*pisq*Ly+32._dp/3._dp*zeta3
     & )*(t**2-s**2)/u**2
     & +((64._dp*Ly-416._dp/3._dp
     & )*Li3x+64._dp*Li3y*Lx+416._dp/
     & 3._dp*Li2x*Lx+(64._dp*Ly+16+16._dp*Ly**2
     & )*Lx**2
     & +(-160._dp/3._dp*Ly**2+32._dp/3._dp*Ly**3+
     & (-80._dp/3._dp-64._dp/3._dp*pisq)*Ly+208._dp/
     & 9._dp*pisq-64._dp*zeta3+80._dp/3._dp)*Lx+320._dp/
     & 9._dp*Ly**3
     & +(-160._dp/3._dp+160._dp/9._dp*pisq+64._dp*zeta3
     & )*Ly+(-16._dp/3._dp*pisq+80._dp/3._dp
     & )*Ly**2-16._dp/3._dp*Ly**4+88._dp/45._dp*pisq**2-200._dp/9._dp*pisq

     & -416._dp/3._dp*zeta3)


      AGTYBu=
     & -12._dp*Lx**2._dp*(t**2+s**2)/u**2+24._dp*Lx*(t**2-s**2)/u**2
     & +8._dp*Ly**2._dp*s**2/t**2
     & +(-16._dp*Lx*Ly+8._dp*Ly**2+8._dp*pisq+8._dp*Lx**2)*t**2/s**2
     & +(44._dp*Li4z+44._dp*Li4y+
     & (-16._dp*Lx-56._dp*Ly+26
     & )*Li3x-88._dp*Li3y*Lx
     & +(-6._dp*Lx**2+(12._dp*Ly-26
     & )*Lx-6._dp*pisq)*Li2x+5._dp*Lx**4+
     & (-20*Ly+3._dp)*Lx**3
     & +(-3._dp/ 2+10._dp/3._dp*pisq-28._dp*Ly+Ly**2
     & )*Lx**2+(-28._dp/3._dp*Ly**3+40*Ly**2+
     & (20._dp/3._dp*pisq+3._dp)*Ly
     & -1._dp/3._dp*pisq-47._dp/ 2+72._dp*zeta3)*Lx+14._dp/
     & 3._dp*Ly**4-80._dp/3._dp*Ly**3+(-3._dp-10._dp/
     & 3._dp*pisq)*Ly**2
     & +(47-56._dp*zeta3-28._dp/3._dp*pisq
     & )*Ly-13._dp/ 2._dp*pisq-46._dp*zeta3-4._dp*pisq*Lu-28._dp/
     & 9._dp*pisq**2
     & +3._dp*Lu+48._dp*zeta3*Lu+187._dp/4._dp)*(t**2+s**2)/(s*t)
     & +(-44._dp*Li4z+44._dp*Li4y+112._dp*Li4x
     & +(-32._dp*Lx+38-24._dp*Ly)*Li3x+(24._dp*Lx-48._dp*Ly+76)*Li3y
     & +(-2._dp*Lx**2+(-38+4._dp*Ly
     & )*Lx+76._dp*Ly-4._dp*Ly**2+6._dp*pisq)*Li2x
     & +1._dp/3._dp*Lx**4+(-4._dp/3._dp*Ly-3._dp
     & )*Lx**3+(-6._dp*pisq-Ly**2+7._dp/ 2-16._dp*Ly
     & )*Lx**2
     & +(-8._dp*Ly**3+54._dp*Ly**2+(12._dp*pisq-7._dp
     & )*Ly-7._dp/3._dp*pisq-56._dp*zeta3+31._dp/ 2)*Lx
     & -4._dp/3._dp*Ly**2._dp*pisq+(24._dp*zeta3-86._dp/
     & 3._dp*pisq)*Ly-206._dp/45._dp*pisq**2-38._dp*zeta3+7._dp/
     & 2._dp*pisq)*(t**2-s**2)/(s*t)
     & +(80*Li4z+80*Li4y+
     & (-96._dp*Ly-32._dp*Lx+152)*Li3x-160*Li3y*Lx
     & +(-8._dp*Lx**2+(-152+16._dp*Ly
     & )*Lx-8._dp*pisq)*Li2x+8._dp*Lx**4+
     & (44._dp/3._dp-32._dp*Ly)*Lx**3
     & +(16._dp/3._dp*pisq-4._dp*Ly**2-104._dp*Ly-24
     & )*Lx**2+(-16._dp*Ly**3+72._dp*Ly**2+
     & (16._dp*pisq-8._dp)*Ly
     & +4._dp/3._dp*pisq+128._dp*zeta3-58
     & )*Lx+8._dp*Ly**4-48._dp*Ly**3+(8._dp-8._dp/3._dp*pisq
     & )*Ly**2
     & +(116-56._dp/3._dp*pisq-96._dp*zeta3
     & )*Ly-4._dp-76._dp/3._dp*pisq-120*zeta3-24._dp/5._dp*pisq**2)


      AGTYCu=(-Lx*Ly+1._dp/ 2._dp*pisq+1._dp/
     & 2._dp*Ly**2+1._dp/ 2._dp*Lx**2)*t**2/s**2 +1._dp/
     & 2._dp*Ly**2._dp*s**2/t**2
     & -5._dp/4._dp*Lx**2._dp*(t**2+s**2)/u**2
     & +5._dp/ 2._dp*Lx*(t**2-s**2)/u**2
     & +(-14._dp*Li4y-14._dp*Li4z+(-59._dp/
     & 6+11._dp*Lx-31._dp*Ly)*Li3x-9._dp*Li3y*Lx
     & +(-4._dp*Lx**2+(59._dp/6._dp+8._dp*Ly
     & )*Lx-4._dp*pisq)*Li2x-1._dp/4._dp*Lx**4+
     & (4._dp/3._dp*Ly+13._dp/18._dp)*Lx**3
     & +(-7._dp/ 2._dp*Ly**2-22._dp/3._dp*Ly-5._dp/
     & 2._dp*pisq+11._dp/4._dp*Lu+14._dp/3._dp)*Lx**2
     & +(-4._dp/3._dp*Ly**3+55._dp/ 2._dp*Ly**2+(35._dp/
     & 6._dp*pisq-33._dp/ 2._dp*Lu)*Ly+2._dp*zeta3+55._dp/
     & 3._dp*Lu-847._dp/54._dp+1._dp/18._dp*pisq)*Lx
     & +2._dp/3._dp*Ly**4-55._dp/3._dp*Ly**3+(-1._dp/
     & 2._dp*pisq+33._dp/ 2._dp*Lu)*Ly**2+(847._dp/
     & 27-28._dp/9._dp*pisq+5._dp*zeta3-110._dp/3._dp*Lu
     & )*Ly
     & -1142._dp/81._dp+61._dp/9._dp*zeta3-13._dp*Lu-166._dp/
     & 9._dp*pisq+47._dp/360*pisq**2-2._dp*zeta3*Lu+143._dp/
     & 18._dp*pisq*Lu+121._dp/12._dp*Lu**2)*(t**2+s**2)/(s*t)
     & +
     & (20*Li4x+14._dp*Li4y-14._dp*Li4z+
     & (-7._dp*Ly+19._dp/ 2-Lx)*Li3x
     & +(-14._dp*Ly+7._dp*Lx+19)*Li3y+
     & (-2._dp*Lx**2-19._dp/ 2._dp*Lx+2._dp/3._dp*pisq+19._dp*Ly
     & )*Li2x
     & -1._dp/4._dp*Lx**4+(2._dp/3._dp*Ly+13._dp/18._dp
     & )*Lx**3+(-3._dp/ 2._dp*Ly**2-38._dp/9._dp-3._dp/
     & 2._dp*pisq-10._dp*Ly+11._dp/4._dp*Lu)*Lx**2
     & +(-4._dp/3._dp*Ly**3+39._dp/ 2._dp*Ly**2+(76._dp/
     & 9+13._dp/6._dp*pisq-11._dp/ 2._dp*Lu)*Ly+77._dp/
     & 36._dp*pisq-31._dp/6._dp-12._dp*zeta3)*Lx
     & -3._dp/ 2._dp*Ly**2._dp*pisq+(7._dp*zeta3-79._dp/
     & 6._dp*pisq)*Ly-19._dp/ 2._dp*zeta3-38._dp/
     & 9._dp*pisq-5._dp/6._dp*pisq**2+11._dp/4._dp*pisq*Lu
     & )*(t**2-s**2)/(s*t)
     & +(-24._dp*Li4y-24._dp*Li4z+
     & (22-60*Ly+20*Lx
     & )*Li3x-20*Li3y*Lx
     & +(-8._dp*Lx**2+(-22+16._dp*Ly
     & )*Lx-8._dp*pisq)*Li2x-2._dp/3._dp*Lx**4+
     & (4._dp/3._dp*Ly+59._dp/9._dp)*Lx**3
     & +(-35._dp*Ly-4._dp/3._dp*pisq-575._dp/
     & 36-6._dp*Ly**2+11._dp*Lu)*Lx**2
     & +(-8._dp/3._dp*Ly**3+131._dp/3._dp*Ly**2+(26._dp/
     & 3._dp*pisq-22._dp*Lu+178._dp/9._dp)*Ly+11._dp*Lu-637._dp/
     & 18+24._dp*zeta3+53._dp/9._dp*pisq)*Lx
     & +4._dp/3._dp*Ly**4-262._dp/9._dp*Ly**3+(-178._dp/
     & 9+2._dp/3._dp*pisq+22._dp*Lu)*Ly**2+(637._dp/
     & 9-28._dp*zeta3-67._dp/9._dp*pisq-22._dp*Lu)*Ly
     & -18._dp*zeta3+2._dp/45._dp*pisq**2+11._dp*pisq*Lu-71._dp/
     & 3._dp*pisq)

      AGTYD1u=
     & (10._dp*Lx*Ly-5._dp*pisq-5._dp*Ly**2-5._dp*Lx**2
     & )*t**2/s**2 -5._dp*Ly**2._dp*s**2/t**2
     & +14._dp*Lx**2._dp*(t**2+s**2)/u**2-28._dp*Lx*(t**2-s**2)/u**2
     & +(-2._dp*Li4y-2._dp*Li4z+
     & (-8._dp-10._dp*Lx+90*Ly
     & )*Li3x+70*Li3y*Lx
     & +(11._dp*Lx**2+(8._dp-22._dp*Ly
     & )*Lx+11._dp*pisq)*Li2x-11._dp/
     & 6._dp*Lx**4+(28._dp/3._dp*Ly-29._dp/3._dp)*Lx**3
     & +(-33._dp/ 2._dp*Lu+53._dp/6._dp+47._dp*Ly+13._dp/
     & 2._dp*Ly**2)*Lx**2+(22._dp/
     & 3._dp*Ly**3-75._dp*Ly**2+(-65._dp/
     & 3-15._dp*pisq+33._dp*Lu)*Ly
     & +389._dp/6._dp-60*zeta3-31._dp/3._dp*pisq-33._dp/ 2._dp*Lu
     & )*Lx-11._dp/3._dp*Ly**4+50*Ly**3+(8._dp/
     & 3._dp*pisq+65._dp/3._dp-33._dp*Lu)*Ly**2
     & +(50*zeta3+29._dp/3._dp*pisq-389._dp/
     & 3+33._dp*Lu)*Ly+587._dp/9._dp*zeta3+326._dp/
     & 9._dp*pisq-52._dp*zeta3*Lu+61._dp/36._dp*pisq**2
     & +1834._dp/ 27._dp*Lu-43417._dp/324._dp-38._dp/3._dp*pisq*Lu
     & )*(t**2+s**2)/(s*t)
     & +
     & (-96._dp*Li4x-50*Li4y+50*Li4z+
     & (18._dp*Lx+26._dp*Ly-38)*Li3x
     & +(52._dp*Ly-26._dp*Lx-76)*Li3y+
     & (5._dp*Lx**2+(-2._dp*Ly+38
     & )*Lx-76._dp*Ly+2._dp*Ly**2-13._dp/3._dp*pisq
     & )*Li2x
     & +1._dp/3._dp*Lx**4+(-2._dp/3._dp*Ly-13._dp/9._dp
     & )*Lx**3+(269._dp/18._dp+6._dp*pisq-11._dp/
     & 2._dp*Lu+28._dp*Ly+7._dp/ 2._dp*Ly**2)*Lx**2
     & +(20._dp/3._dp*Ly**3-66._dp*Ly**2+(-269._dp/
     & 9-31._dp/3._dp*pisq+11._dp*Lu)*Ly+8._dp/
     & 9._dp*pisq+52._dp*zeta3+33._dp/ 2._dp*Lu-31._dp/ 2)*Lx

     & +11._dp/3._dp*Ly**2._dp*pisq+(-26._dp*zeta3+122._dp/
     & 3._dp*pisq)*Ly+38._dp*zeta3+178._dp/
     & 45._dp*pisq**2+269._dp/18._dp*pisq-11._dp/ 2._dp*pisq*Lu
     & )*(t**2-s**2)/(s*t)
     & +(8._dp*Li4y+8._dp*Li4z+
     & (-120+168._dp*Ly-24._dp*Lx
     & )*Li3x+120*Li3y*Lx
     & +(20*Lx**2+(120-40*Ly
     & )*Lx+20*pisq)*Li2x-8._dp/3._dp*Lx**4+
     & (40._dp/3._dp*Ly-184._dp/9._dp)*Lx**3
     & +(472._dp/9._dp+14._dp*Ly**2+122._dp*Ly-22._dp*Lu
     & )*Lx**2+(40._dp/3._dp*Ly**3-370._dp/3._dp*Ly**2+
     & (-76._dp/3._dp*pisq+44._dp*Lu-320._dp/9._dp)*Ly
     & -112._dp/9._dp*pisq+898._dp/9._dp-112._dp*zeta3-22._dp*Lu
     & )*Lx-20._dp/3._dp*Ly**4+740._dp/9._dp*Ly**3+(320._dp/
     & 9-44._dp*Lu)*Ly**2
     & +(44._dp*Lu+218._dp/9._dp*pisq+104._dp*zeta3-1796._dp/
     & 9)*Ly+96._dp*zeta3+104._dp/
     & 45._dp*pisq**2+4._dp-22._dp*pisq*Lu+60*pisq)

      AGTYE1u=(11._dp/6._dp*Lx**3+(3._dp*Lu-23._dp/
     & 6-6._dp*Ly)*Lx**2+(7._dp*Ly**2+
     & (-6._dp*Lu+20._dp/3._dp)*Ly+11._dp/6._dp*pisq-22._dp/
     & 3+3._dp*Lu)*Lx
     & -14._dp/3._dp*Ly**3+(-20._dp/3._dp+6._dp*Lu
     & )*Ly**2+(-6._dp*Lu+44._dp/3._dp+4._dp/3._dp*pisq
     & )*Ly
     & -2._dp/9._dp*zeta3+8._dp/3._dp*pisq*Lu-121._dp/
     & 18._dp*pisq +3401._dp/162._dp-328._dp/ 27._dp*Lu)
     & *(t**2+s**2)/(s*t)
     & +(11._dp/18._dp*Lx**3+(-83._dp/
     & 18+LU-2._dp*Ly)*Lx**2+(2._dp*Ly**2+
     & (83._dp/9._dp-2._dp*Lu)*Ly-3._dp*Lu+5._dp+11._dp/
     & 18._dp*pisq)*Lx
     & -2._dp*pisq*Ly+1._dp/18._dp*pisq*
     & (-83+18._dp*Lu))*(t**2-s**2)/(s*t)
     & +(22._dp/9._dp*Lx**3+(-46._dp/
     & 9+4._dp*Lu-8._dp*Ly)*Lx**2+(28._dp/
     & 3._dp*Ly**2+(80._dp/9._dp-8._dp*Lu)*Ly+22._dp/
     & 9._dp*pisq-76._dp/9._dp+4._dp*Lu)*Lx
     & -56._dp/9._dp*Ly**3+(-80._dp/9._dp+8._dp*Lu
     & )*Ly**2+(152._dp/9._dp-8._dp*Lu+16._dp/9._dp*pisq
     & )*Ly+2._dp*pisq*(-5._dp+2._dp*Lu)
     & )

      AGTYE2u=(4._dp/3._dp*Li3x-4._dp/
     & 3._dp*Li2x*Lx-11._dp/36._dp*Lx**3+(4._dp/
     & 3._dp*Ly-13._dp/18._dp-1._dp/ 2._dp*Lu)*Lx**2
     & +(-7._dp/ 2._dp*Ly**2+(1._dp/3._dp+3._dp*Lu
     & )*Ly+1._dp/36._dp*pisq-31._dp/6._dp*Lu+40._dp/9._dp
     & )*Lx
     & +7._dp/3._dp*Ly**3+(-1._dp/3._dp-3._dp*Lu
     & )*Ly**2+(31._dp/3._dp*Lu-5._dp/9._dp*pisq-80._dp/
     & 9)*Ly
     & -25._dp/9._dp*zeta3-13._dp/9._dp*pisq*Lu-11._dp/
     & 3._dp*Lu**2+206._dp/ 27._dp*Lu+65._dp/81._dp+487._dp/108._dp*pisq
     & )*(t**2+s**2)/(s*t)
     & +(-11._dp/36._dp*Lx**3+(14._dp/9._dp-1._dp/
     & 2._dp*Lu+Ly)*Lx**2 +(-Ly**2+
     & (-28._dp/9._dp+LU)*Ly-11._dp/36._dp*pisq-1._dp/3._dp
     & )*Lx
     & +pisq*Ly-1._dp/18._dp*pisq*
     & (-28+9._dp*Lu))*(t**2-s**2)/(s*t)
     & -Lx**2._dp*(t**2+s**2)/u**2+2._dp*Lx*(t**2-s**2)/u**2
     & +(-11._dp/9._dp*Lx**3+(14._dp/
     & 9-2._dp*Lu+4._dp*Ly)*Lx**2
     & +(-14._dp/3._dp*Ly**2+(-40._dp/
     & 9+4._dp*Lu)*Ly+38._dp/9._dp-2._dp*Lu-11._dp/
     & 9._dp*pisq)*Lx+(-76._dp/9._dp+4._dp*Lu-8._dp/
     & 9._dp*pisq)*Ly
     & +28._dp/9._dp*Ly**3+(40._dp/9._dp-4._dp*Lu
     & )*Ly**2-pisq*(-5._dp+2._dp*Lu))


      AGTYF1u=(5._dp/36._dp*Lx**2+(1._dp/
     & 3._dp*Lu-10._dp/ 27-1._dp/3._dp*Ly)*Lx+1._dp/
     & 3._dp*Ly**2+(20._dp/ 27-2._dp/3._dp*Lu)*Ly

     & -20._dp/ 27._dp*Lu-13._dp/108._dp*pisq+1._dp/3._dp*Lu**2
     & )*(t**2+s**2)/(s*t)

      AGTYD2u=
     & (-6._dp*pisq-6._dp*Ly**2+12._dp*Lx*Ly-6._dp*Lx**2
     & )*t**2/s**2 -6._dp*Ly**2._dp*s**2/t**2
     & +6._dp*Lx**2._dp*(t**2+s**2)/u**2-12._dp*Lx*(t**2-s**2)/u**2
     & +(-4._dp*Li4y-4._dp*Li4z+
     & (-4._dp-4._dp*Lx+36._dp*Ly
     & )*Li3x+28._dp*Li3y*Lx
     & +(6._dp*Lx**2-12._dp*Lx*Ly+6._dp*pisq+4._dp*Lx
     & )*Li2x+20*Ly**3+(107._dp/
     & 3+3._dp*Lx**2-22._dp*Lu-30*Lx)*Ly**2
     & +(4._dp*Lx**3+24._dp*Lx**2+
     & (-4._dp*pisq+22._dp*Lu-107._dp/3._dp
     & )*Lx-57+36._dp*zeta3+22._dp*Lu-2._dp/3._dp*pisq
     & )*Ly
     & -Lx**4-19._dp/3._dp*Lx**3+(-11._dp*Lu+107._dp/6._dp+1._dp/
     & 3._dp*pisq)*Lx**2+
     & (-11._dp*Lu-32._dp*zeta3+57._dp/ 2-20._dp/3._dp*pisq
     & )*Lx
     & -43417._dp/324._dp-43._dp/
     & 6._dp*pisq*Lu-52._dp*zeta3*Lu+515._dp/9._dp*zeta3+251._dp/
     & 9._dp*pisq+1141._dp/ 27._dp*Lu+65._dp/36._dp*pisq**2
     & )*(t**2+s**2)/(s*t)
     & +
     & (-20*Li4y-48._dp*Li4x+20*Li4z+
     & (12._dp*Ly-16+12._dp*Lx)*Li3x
     & +
     & (4._dp*Ly**2+(-4._dp*Lx-32
     & )*Ly+2._dp*Lx**2+16._dp*Lx-10._dp/3._dp*pisq
     & )*Li2x
     & +(-12._dp*Lx+24._dp*Ly-32
     & )*Li3y+4._dp*Ly**3._dp*Lx+(-94._dp/
     & 3._dp*Lx+Lx**2+4._dp/3._dp*pisq)*Ly**2
     & +(46._dp/3._dp*Lx**2+(-20._dp/3._dp*pisq+22._dp/
     & 3._dp*Lu-251._dp/9._dp)*Lx-12._dp*zeta3+62._dp/
     & 3._dp*pisq)*Ly
     & -13._dp/9._dp*Lx**3+(251._dp/18._dp+3._dp*pisq-11._dp/
     & 3._dp*Lu)*Lx**2+(-16._dp/
     & 9._dp*pisq+24._dp*zeta3-57._dp/ 2+11._dp*Lu)*Lx
     & +16._dp*zeta3+251._dp/18._dp*pisq-11._dp/
     & 3._dp*pisq*Lu+94._dp/45._dp*pisq**2)*(t**2-s**2)/(s*t)
     & +(-16._dp*Li4y-16._dp*Li4z+
     & (-64+48._dp*Ly)*Li3x+48._dp*Li3y*Lx
     & +
     & (-16._dp*Lx*Ly+64._dp*Lx+8._dp*Lx**2+8._dp*pisq
     & )*Li2x+272._dp/9._dp*Ly**3
     & +(4._dp*Lx**2+344._dp/9._dp-136._dp/3._dp*Lx-88._dp/
     & 3._dp*Lu)*Ly**2+(8._dp*Lx**3+184._dp/
     & 3._dp*Lx**2+(88._dp/3._dp*Lu-344._dp/9._dp-16._dp/
     & 3._dp*pisq)*Lx
     & +88._dp/3._dp*Lu+32._dp/9._dp*pisq+48._dp*zeta3-1076._dp/
     & 9)*Ly-2._dp*Lx**4-112._dp/9._dp*Lx**3+
     & (-44._dp/3._dp*Lu+298._dp/9._dp)*Lx**2
     & +(-76._dp/9._dp*pisq-44._dp/3._dp*Lu+538._dp/
     & 9-48._dp*zeta3)*Lx+48._dp*pisq+48._dp*zeta3-44._dp/
     & 3._dp*pisq*Lu+92._dp/45._dp*pisq**2)

      AGTYE3u=(-8._dp/3._dp*Ly**3+(4._dp*Lu-26._dp/
     & 3+4._dp*Lx)*Ly**2+(-4._dp*Lx**2+
     & (26._dp/3._dp-4._dp*Lu)*Lx+4._dp/
     & 3._dp*pisq+10._dp-4._dp*Lu)*Ly
     & +4._dp/3._dp*Lx**3+(2._dp*Lu-13._dp/3._dp
     & )*Lx**2+(4._dp/3._dp*pisq+2._dp*Lu-5._dp
     & )*Lx
     & -2._dp/9._dp*zeta3+5._dp/3._dp*pisq*Lu+3401._dp/
     & 162-56._dp/9._dp*pisq-202._dp/ 27._dp*Lu)*(t**2+s**2)/(s*t)
     & +(4._dp/3._dp*Ly**2._dp*Lx+(-4._dp/
     & 3._dp*Lx**2+(-4._dp/3._dp*Lu+74._dp/9._dp
     & )*Lx-4._dp/3._dp*pisq)*Ly+4._dp/9._dp*Lx**3+
     & (-37._dp/9._dp+2._dp/3._dp*Lu)*Lx**2
     & +(-2._dp*Lu+5._dp+4._dp/9._dp*pisq)*Lx+1._dp/
     & 9._dp*pisq*(-37+6._dp*Lu))*(t**2-s**2)/(s*t)
     & +(-32._dp/9._dp*Ly**3+(16._dp/3._dp*Lu+16._dp/
     & 3._dp*Lx-104._dp/9._dp)*Ly**2+(-16._dp/
     & 3._dp*Lx**2+(-16._dp/3._dp*Lu+104._dp/9._dp)*Lx

     & +152._dp/9._dp-16._dp/3._dp*Lu+16._dp/9._dp*pisq
     & )*Ly+16._dp/9._dp*Lx**3+(8._dp/3._dp*Lu-52._dp/9._dp
     & )*Lx**2+(8._dp/3._dp*Lu-76._dp/9._dp+16._dp/
     & 9._dp*pisq)*Lx
     & +4._dp/3._dp*pisq*(2._dp*Lu-7._dp)
     & )

      AGTYE4u=(8._dp/3._dp*Lx**3+(8._dp*Lu-26._dp/
     & 3-8._dp*Ly)*Lx**2+(8._dp*Ly**2+(52._dp/
     & 3-16._dp*Lu)*Ly+8._dp*Lu+8._dp/3._dp*pisq-10._dp
     & )*Lx
     & -16._dp/3._dp*Ly**3+(-52._dp/3._dp+16._dp*Lu
     & )*Ly**2+(8._dp/3._dp*pisq-16._dp*Lu+20
     & )*Ly
     & -112._dp/9._dp*pisq-908._dp/ 27._dp*Lu-4._dp/
     & 9._dp*zeta3+22._dp/3._dp*pisq*Lu+3401._dp/81._dp
     & )*(t**2+s**2)/(s*t)
     & +(8._dp/9._dp*Lx**3+(8._dp/3._dp*Lu-74._dp/
     & 9-8._dp/3._dp*Ly)*Lx**2
     & +(8._dp/
     & 3._dp*Ly**2+(-16._dp/3._dp*Lu+148._dp/9._dp
     & )*Ly+8._dp/9._dp*pisq+10._dp-8._dp*Lu)*Lx
     & -8._dp/3._dp*pisq*Ly+2._dp/9._dp*pisq*
     & (12._dp*Lu-37))*(t**2-s**2)/(s*t)
     & +(32._dp/9._dp*Lx**3+(-104._dp/9._dp-32._dp/
     & 3._dp*Ly+32._dp/3._dp*Lu)*Lx**2
     & +(32._dp/3._dp*Ly**2+(-64._dp/3._dp*Lu+208._dp/
     & 9)*Ly-152._dp/9._dp+32._dp/3._dp*Lu+32._dp/
     & 9._dp*pisq)*Lx
     & -64._dp/9._dp*Ly**3+(-208._dp/9._dp+64._dp/3._dp*Lu
     & )*Ly**2+(-64._dp/3._dp*Lu+304._dp/9._dp+32._dp/
     & 9._dp*pisq)*Ly
     & +8._dp/3._dp*pisq*(-7._dp+4._dp*Lu))

      AGTYF2u=(4._dp/ 27._dp*pisq+32._dp/9._dp*Lu**2-160._dp/
     & 27._dp*Lu)*(t**2+s**2)/(s*t)



      AGTYG1s=
     & (14._dp*Lx**4+28._dp*Lx**3+8._dp*Lx**2._dp*Ly**2
     & +56._dp*Lx**2._dp*pisq-48._dp*Lx**2+12._dp*Lx**2._dp*Ly
     & +32._dp*Lx*Ly*pisq+80*pisq*Lx
     & +2._dp*Ly**4+12._dp*Ly**3-10._dp*Ly**2+8._dp*Ly**2._dp*pisq
     & +26._dp*pisq+24._dp*pisq*Ly-84._dp*Ly+102
     & )*t/u
     & +8._dp*Lx*
     & (Lx**3+Lx**2+4._dp*pisq*Lx+2._dp*pisq
     & )*t**2/u**2+2._dp*Lx**2._dp*(Lx**2+4._dp*pisq
     & )*t**3/u**3
     & +
     & (32._dp*Lx**3+8._dp*Lx**2._dp*Ly**2+80*pisq*Lx+32._dp*Lx*Ly*pisq
     & +8._dp*Ly**2._dp*Lx

     & +8._dp*Ly**4+32._dp*Ly**2._dp*pisq-32._dp*Ly**2-4._dp-56._dp*Ly
     & +24._dp*pisq
     & )
!      + t leftrightarrow u

      AGTYG2s=(-10._dp*Lx**4-106._dp/
     & 3._dp*Lx**3-8._dp*Lx**3._dp*Ly-2._dp*Lx**2._dp*Ly**2
     & -52._dp*Lx**2._dp*pisq-44._dp/3._dp*Lx**2._dp*Ls
     & +6._dp*Lx**2-40._dp/3._dp*Lx**2._dp*Ly
     & -4._dp*Ly**3._dp*Lx-56._dp/
     & 3._dp*Ly**2._dp*Lx-32._dp*Lx*Ly*pisq+8._dp*Lx*Ly-80*pisq*Lx
     & +140._dp/3._dp*Lx-20._dp/3._dp*Ly**3-22._dp/3._dp*Ly**2._dp*Ls
     & -20*Ly**2-6._dp*Ly**2._dp*pisq-18._dp*pisq*Ly-22._dp*Ly*Ls
     & +140._dp/3._dp*Ly-4._dp+154._dp/3._dp*Ls-40*pisq)*t/u
     & -8._dp*Lx*
     & (Lx**3+Lx**2+4._dp*pisq*Lx+2._dp*pisq
     & )*t**2/u**2-2._dp*Lx**2._dp*(Lx**2+4._dp*pisq
     & )*t**3/u**3
     & +
     & (-8._dp*Lx**3._dp*Ly-32._dp*Lx**3+4._dp*Lx**2._dp*pisq
     & -4._dp*Lx**2._dp*Ly**2-52._dp/
     & 3._dp*Lx**2._dp*Ly-8._dp*Ly**2._dp*Lx-40._dp/3._dp*Lx*Ly
     & -32._dp*Lx*Ly*pisq-76._dp*pisq*Lx-44._dp/
     & 3._dp*Lx*Ls-4._dp*Ly**4+8._dp/3._dp*Ly**3+8._dp/
     & 3._dp*Ly**2-44._dp/3._dp*Ly**2._dp*Ls-32._dp*Ly**2._dp*pisq
     & +28._dp*Ly+4._dp-24._dp*pisq)
!      + t leftrightarrow u

      AGTYG3s=(2._dp*Lx**4+2._dp*Lx**3._dp*Ly+22._dp/
     & 3._dp*Lx**3+11._dp/
     & 3._dp*Lx**2._dp*Ls+2._dp*Lx**2._dp*Ly**2+10._dp*Lx**2._dp*Ly
     & +68._dp/9._dp*Lx**2+13._dp*Lx**2._dp*pisq
     & +20._dp/3._dp*Ly**2._dp*Lx+6._dp*Lx*Ly*pisq+100._dp/
     & 9._dp*Lx*Ly+22._dp/3._dp*Lx*Ly*Ls+110._dp/
     & 9._dp*Lx*Ls+50._dp/3._dp*pisq*Lx+50._dp/
     & 9._dp*Ly**2
     & +2._dp*Ly**2._dp*pisq+8._dp/3._dp*pisq*Ly+110._dp/
     & 9._dp*Ly*Ls+2+121._dp/18._dp*Ls**2-11._dp/
     & 3._dp*pisq*Ls+1._dp/ 2._dp*pisq**2+13._dp/ 2._dp*pisq
     & )*t/u
     & +2._dp*Lx*
     & (Lx**3+Lx**2+4._dp*pisq*Lx+2._dp*pisq)*t**2/u**2+1._dp/
     & 2._dp*Lx**2._dp*(Lx**2+4._dp*pisq)*t**3/u**3
     & +(4._dp*Lx**3._dp*Ly+8._dp*Lx**3-2._dp*Lx**2._dp*pisq+26._dp/
     & 3._dp*Lx**2._dp*Ly+2._dp*Ly**2._dp*Lx+8._dp*Lx*Ly*pisq+20._dp/
     & 3._dp*Lx*Ly+22._dp/3._dp*Lx*Ls
     & +18._dp*pisq*Lx-4._dp/3._dp*Ly**3+20._dp/3._dp*Ly**2+22._dp/
     & 3._dp*Ly**2._dp*Ls+8._dp*Ly**2._dp*pisq+6._dp*pisq)
!      + t leftrightarrow u

      AGTYX1s=2._dp/3._dp*(3._dp*Ly+Ly**2-7._dp+2._dp*Lx**2
     & )*(2._dp*Ls+Ly+Lx)*t/u
     & +(4._dp/3._dp*Lx**2._dp*Ly+8._dp/3._dp*Lx*Ls+4._dp/
     & 3._dp*Lx*Ly+4._dp/3._dp*Ly**3+4._dp/3._dp*Ly**2+8._dp/
     & 3._dp*Ly**2._dp*Ls)
!      + t leftrightarrow u

      AGTYX2s=-1._dp/9._dp*(2._dp*Ls+Ly+Lx
     & )*
     & (11._dp*Ls+3._dp*Lx**2+6._dp*Lx*Ly+10._dp*Ly+10._dp*Lx-3._dp*pisq
     & )*t/u
     & +(-2._dp/3._dp*Lx*Ly-2._dp/3._dp*Ly**2-2._dp/
     & 3._dp*Lx**2._dp*Ly-4._dp/3._dp*Ly**2._dp*Ls-4._dp/
     & 3._dp*Lx*Ls-2._dp/3._dp*Ly**3)
!      + t leftrightarrow u

      AGTYX3s=1._dp/18._dp*(2._dp*Ls+Ly
     & +Lx)**2._dp*t/u
!      + t leftrightarrow u

      AGTYX4s= (-32._dp/3._dp*Ly*pisq+32._dp/
     & 3._dp*Ls*Lx**2+16._dp*Ls*Ly-64._dp/
     & 3._dp*Lx*pisq-16._dp*pisq+16._dp/3._dp*Ls*Ly**2-112._dp/
     & 3._dp*Ls)*t/u
     & +(32._dp/3._dp*Ls*Lx**2+32._dp/3._dp*Ls*Ly-64._dp/
     & 3._dp*Lx*pisq-32._dp/3._dp*pisq)
!      + t leftrightarrow u

      AGTYX5s= (32._dp/9._dp*Ls**2+32._dp/9._dp*pisq
     & )*t/u
!       + t leftrightarrow u





      AGTYG1u=(8._dp*Lx**4+(-32._dp*Ly+8._dp
     & )*Lx**3+(48._dp*Ly**2-24._dp*Ly+16._dp*pisq
     & )*Lx**2
     & +
     & (-32._dp*Ly**3+24._dp*Ly**2-32._dp*pisq*Ly+8._dp*pisq
     & )*Lx
     & +8._dp*Ly**4-8._dp*Ly**3+16._dp*Ly**2._dp*pisq-8._dp*pisq*Ly
     & +8._dp*pisq**2
     & )*t**2/s**2
     & +(8._dp**Ly**4-8._dp*Ly**3
     & +32._dp*pisq*Ly**2-16._dp*pisq*Ly)*s**2/t**2
     & +(2._dp*Lx**4-8._dp*Lx**3._dp*Ly+
     & (4._dp*pisq+12._dp*Ly**2)*Lx**2+
     & (-8._dp*Ly**3-8._dp*pisq*Ly
     & )*Lx
     & +2._dp*Ly**4+4._dp*Ly**2._dp*pisq+2._dp*pisq**2)*t**3/s**3
     & +(2._dp*Ly**4 + 8._dp*pisq*Ly**2) *s**3/t**3
     & +(8._dp*Lx**4+(-32._dp*Ly+20
     & )*Lx**3+(16._dp*pisq+56._dp*Ly**2-66._dp*Ly-29
     & )*Lx**2
     & +(-48._dp*Ly**3+78._dp*Ly**2+
     & (58-32._dp*pisq)*Ly-42+20*pisq)*Lx
     & +24._dp*Ly**4-52._dp*Ly**3+(56._dp*pisq-58
     & )*Ly**2+(-66._dp*pisq+84
     & )*Ly
     & +102+8._dp*pisq**2-29._dp*pisq)*(t**2+s**2)/(s*t)
     & +(6._dp*Lx**4+(-24._dp*Ly+8._dp
     & )*Lx**3+(36._dp*Ly**2-19-30*Ly+12._dp*pisq
     & )*Lx**2
     & +(-24._dp*Ly**3+30*Ly**2+
     & (38-24._dp*pisq)*Ly+8._dp*pisq+42
     & )*Lx
     & -12._dp*Ly**2._dp*pisq+2._dp*pisq*Ly+3._dp*pisq*
     & (2._dp*pisq-3._dp))*(t**2-s**2)/(s*t)
     & +(8._dp*Lx**4+(32-32._dp*Ly
     & )*Lx**3+(-104._dp*Ly+64._dp*Ly**2-32+16._dp*pisq
     & )*Lx**2
     & +(-64._dp*Ly**3+120*Ly**2+
     & (-32._dp*pisq+64)*Ly-56+32._dp*pisq)*Lx

     & +32._dp*Ly**4-80*Ly**3+64._dp*(pi-1)*
     & (pi+1)*Ly**2
     & +(-104._dp*pisq+112
     & )*Ly-8._dp+8._dp*pisq**2-32._dp*pisq)

      AGTYG2u=(-8._dp*Lx**4+(32._dp*Ly-8._dp
     & )*Lx**3+(-48._dp*Ly**2+24._dp*Ly-16._dp*pisq
     & )*Lx**2
     & +
     & (32._dp*Ly**3-24._dp*Ly**2+32._dp*pisq*Ly-8._dp*pisq
     & )*Lx
     & -8._dp*Ly**4+8._dp*Ly**3-16._dp*Ly**2._dp*pisq+8._dp*pisq*Ly
     & -8._dp*pisq**2)*t**2/s**2
     & -(8._dp**Ly**4-8._dp*Ly**3
     & +32._dp*pisq*Ly**2-16._dp*pisq*Ly)*s**2/t**2
     & +(-2._dp*Lx**4+8._dp*Lx**3._dp*Ly+
     & (-4._dp*pisq-12._dp*Ly**2)*Lx**2+
     & (8._dp*Ly**3+8._dp*pisq*Ly)*Lx
     & -2._dp*Ly**4-4._dp*Ly**2._dp*pisq-2._dp*pisq**2)*t**3/s**3
     & -(2._dp*Ly**4 + 8._dp*pisq*Ly**2) *s**3/t**3
     & +(-5._dp*Lx**4+(-21+26._dp*Ly
     & )*Lx**3+
     & (-11._dp*Lu-7._dp-50*Ly**2-13._dp*pisq+79._dp*Ly
     & )*Lx**2
     & +(48._dp*Ly**3-111._dp*Ly**2+
     & (6._dp+44._dp*pisq+22._dp*Lu)*Ly+140._dp/
     & 3-30*pisq-11._dp*Lu)*Lx
     & -24._dp*Ly**4+74._dp*Ly**3+
     & (-6._dp-56._dp*pisq-22._dp*Lu)*Ly**2+(-280._dp/
     & 3+85._dp*pisq+22._dp*Lu)*Ly
     & +154._dp/
     & 3._dp*Lu-4._dp+7._dp*pisq-11._dp*pisq*Lu-8._dp*pisq**2
     & )*(t**2+s**2)/(s*t)
     & +(-5._dp*Lx**4+(22._dp*Ly-43._dp/3._dp
     & )*Lx**3+(121._dp/
     & 3._dp*Ly-11._dp*pisq-36._dp*Ly**2+13-11._dp/3._dp*Lu
     & )*Lx**2
     & +(24._dp*Ly**3-121._dp/3._dp*Ly**2+(-26+22._dp/
     & 3._dp*Lu+20*pisq)*Ly-52._dp/
     & 3._dp*pisq+11._dp*Lu
     & )*Lx
     & +12._dp*Ly**2._dp*pisq-5._dp*pisq*Ly-1._dp/
     & 3._dp*pisq*(18._dp*pisq+11._dp*Lu-3._dp)
     & )*(t**2-s**2)/(s*t)
     & +(-4._dp*Lx**4+(-88._dp/3._dp+24._dp*Ly
     & )*Lx**3+(-44._dp/3._dp*Lu+8._dp/
     & 3-56._dp*Ly**2-12._dp*pisq+340._dp/3._dp*Ly)*Lx**2
     & +(64._dp*Ly**3-164._dp*Ly**2+(64._dp/
     & 3+48._dp*pisq+88._dp/3._dp*Lu)*Ly+28-124._dp/
     & 3._dp*pisq-44._dp/3._dp*Lu)*Lx
     & -32._dp*Ly**4+328._dp/3._dp*Ly**3+(-64._dp/
     & 3-64._dp*pisq-88._dp/3._dp*Lu)*Ly**2+
     & (-56+364._dp/3._dp*pisq+88._dp/3._dp*Lu)*Ly
     & +8._dp+8._dp/3._dp*pisq-44._dp/3._dp*pisq*Lu-8._dp*pisq**2
     & )

      AGTYG3u=(2._dp*Lx**4+(-8._dp*Ly+2
     & )*Lx**3+(12._dp*Ly**2-6._dp*Ly+4._dp*pisq
     & )*Lx**2
     & +
     & (-8._dp*Ly**3+6._dp*Ly**2-8._dp*pisq*Ly+2._dp*pisq
     & )*Lx
     & +2._dp*Ly**4-2._dp*Ly**3+4._dp*Ly**2._dp*pisq-2._dp*pisq*Ly
     & +2._dp*pisq**2
     & )*t**2/s**2
     & +(2._dp**Ly**4-2._dp*Ly**3
     & +8._dp*pisq*Ly**2-4._dp*pisq*Ly)*s**2/t**2
     & +(1._dp/ 2._dp*Lx**4-2._dp*Lx**3._dp*Ly+
     & (pisq+3._dp*Ly**2)*Lx**2+
     & (-2._dp*Ly**3-2._dp*pisq*Ly)*Lx
     & +1._dp/ 2._dp*Ly**4+Ly**2._dp*pisq+1._dp/ 2._dp*pisq**2
     & )*t**3/s**3 +(1._dp/ 2._dp*Ly**4 + 2._dp*pisq*Ly**2)
     & *s**3/t**3
     & +(Lx**4+(11._dp/3._dp-5._dp*Ly
     & )*Lx**3+(11._dp/6._dp*Lu+59._dp/9._dp+11._dp*Ly**2+9._dp/
     & 2._dp*pisq-58._dp/3._dp*Ly)*Lx**2
     & +(-12._dp*Ly**3+36._dp*Ly**2+
     & (-14._dp*pisq-218._dp/9._dp-11._dp*Lu)*Ly+110._dp/
     & 9._dp*Lu+41._dp/3._dp*pisq)*Lx
     & +6._dp*Ly**4-24._dp*Ly**3+(11._dp*Lu+218._dp/
     & 9+14._dp*pisq)*Ly**2+(-26._dp*pisq-220._dp/
     & 9._dp*Lu)*Ly
     & +121._dp/18._dp*Lu**2+2._dp*pisq**2+11._dp/ 2._dp*pisq*Lu+59._dp/
     & 9._dp*pisq+2)*(t**2+s**2)/(s*t)
     & +(Lx**4+(11._dp/3._dp-5._dp*Ly
     & )*Lx**3+(1-38._dp/3._dp*Ly+5._dp/
     & 2._dp*pisq+9._dp*Ly**2+11._dp/6._dp*Lu)*Lx**2
     & +(-6._dp*Ly**3+38._dp/3._dp*Ly**2+
     & (-4._dp*pisq-11._dp/3._dp*Lu-2)*Ly+11._dp/
     & 3._dp*pisq
     & )*Lx
     & -3._dp*Ly**2._dp*pisq+2._dp*pisq*Ly+1._dp/
     & 6._dp*pisq*(9._dp*pisq+11._dp*Lu-6._dp)
     & )*(t**2-s**2)/(s*t)
     & +((20._dp/3._dp-4._dp*Ly)*Lx**3+
     & (12._dp*Ly**2-92._dp/3._dp*Ly+20._dp/3._dp+22._dp/
     & 3._dp*Lu+2._dp*pisq)*Lx**2
     & +(-16._dp*Ly**3+52._dp*Ly**2+(-44._dp/
     & 3._dp*Lu-80._dp/3._dp-16._dp*pisq)*Ly+22._dp/
     & 3._dp*Lu+38._dp/3._dp*pisq)*Lx
     & +8._dp*Ly**4-104._dp/3._dp*Ly**3+(16._dp*pisq+80._dp/
     & 3+44._dp/3._dp*Lu)*Ly**2+(-104._dp/
     & 3._dp*pisq-44._dp/3._dp*Lu)*Ly
     & +2._dp/
     & 3._dp*pisq*(3._dp*pisq+10._dp+11._dp*Lu)
     & )

      AGTYX1u=(Lx**3+(2._dp*Lu-4._dp*Ly+1
     & )*Lx**2+(6._dp*Ly**2+(-4._dp-4._dp*Lu
     & )*Ly+pisq+2._dp*Lu-14._dp/3._dp)*Lx
     & -4._dp*Ly**3+(4._dp+4._dp*Lu)*Ly**2+
     & (-4._dp*Lu+28._dp/3._dp-4._dp*pisq
     & )*Ly+2._dp*pisq*Lu-28._dp/3._dp*Lu+pisq
     & )*(t**2+s**2)/(s*t)
     & +(1._dp/3._dp*Lx**3+(-4._dp/3._dp*Ly-1+2._dp/
     & 3._dp*Lu)*Lx**2+(4._dp/3._dp*Ly**2+
     & (-4._dp/3._dp*Lu+2)*Ly+1._dp/
     & 3._dp*pisq-2._dp*Lu)*Lx
     & +1._dp/3._dp*pisq*(2._dp*Lu+3._dp)
     & )*(t**2-s**2)/(s*t)
     & +(4._dp/3._dp*Lx**3+(-16._dp/3._dp*Ly+8._dp/
     & 3._dp*Lu+4._dp/3._dp)*Lx**2+(8._dp*Ly**2+
     & (-16._dp/3._dp*Lu-16._dp/3._dp)*Ly+4._dp/
     & 3._dp*pisq+8._dp/3._dp*Lu)*Lx
     & -16._dp/3._dp*Ly**3+(16._dp/3._dp+16._dp/3._dp*Lu
     & )*Ly**2+(-16._dp/3._dp*Lu-16._dp/3._dp*pisq
     & )*Ly+4._dp/3._dp*pisq*(2._dp*Lu+1)
     & )

      AGTYX2u=(-1._dp/6._dp*Lx**3+(-1._dp/
     & 3._dp*Lu+4._dp/3._dp*Ly-10._dp/9._dp)*Lx**2+
     & (-3._dp*Ly**2+(2._dp*Lu+40._dp/9._dp)*Ly-7._dp/
     & 6._dp*pisq-31._dp/9._dp*Lu)*Lx
     & +2._dp*Ly**3+(-2._dp*Lu-40._dp/9._dp
     & )*Ly**2+(2._dp*pisq+62._dp/9._dp*Lu
     & )*Ly-22._dp/9._dp*Lu**2-pisq*Lu-10._dp/9._dp*pisq
     & )*(t**2+s**2)/(s*t)
     & +(-1._dp/6._dp*Lx**3+(2._dp/3._dp*Ly-1._dp/
     & 3._dp*Lu)*Lx**2+(-2._dp/3._dp*Ly**2-1._dp/
     & 6._dp*pisq+2._dp/3._dp*Ly*Lu)*Lx-1._dp/
     & 3._dp*pisq*Lu)*(t**2-s**2)/(s*t)
     & +(-2._dp/3._dp*Lx**3+(8._dp/3._dp*Ly-4._dp/
     & 3._dp*Lu-2._dp/3._dp)*Lx**2+(-4._dp*Ly**2+
     & (8._dp/3._dp*Lu+8._dp/3._dp)*Ly-2._dp/3._dp*pisq-4._dp/
     & 3._dp*Lu)*Lx
     & +8._dp/3._dp*Ly**3+(-8._dp/3._dp-8._dp/3._dp*Lu
     & )*Ly**2+(8._dp/3._dp*Lu+8._dp/3._dp*pisq
     & )*Ly-2._dp/3._dp*pisq*(2._dp*Lu+1)
     & )

      AGTYX3u=(1._dp/18._dp*Lx**2+(-2._dp/9._dp*Ly+2._dp/9._dp*Lu)*Lx
     & -4._dp/9._dp*Ly*Lu+1._dp/18._dp*pisq+2._dp/9._dp*Lu**2
     & +2._dp/9._dp*Ly**2)*(t**2+s**2)/(s*t)


      AGTYX4u=32._dp/3._dp*Lu*
     & (Lx+Lx**2-2._dp*Ly-2._dp*Ly*Lx+2._dp*Ly**2+pisq
     & )
     & +8._dp/3._dp*Lu*
     & (3._dp*pisq-6._dp*Ly-6._dp*Ly*Lx-14+3._dp*Lx**2+6._dp*Ly**2
     & +3._dp*Lx)*(t**2+s**2)/(s*t)
     & -8._dp/3._dp*Lu*(-pisq+2._dp*Ly*Lx-Lx**2+3._dp*Lx)
     & *(t**2-s**2)/(s*t)

      AGTYX5u=32._dp/9._dp*Lu**2._dp*(t**2+s**2)/(s*t)

      stop
      end
