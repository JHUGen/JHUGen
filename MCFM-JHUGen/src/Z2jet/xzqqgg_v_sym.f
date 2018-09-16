      subroutine xzqqgg_v_sym(mqqb_ax)
      implicit none
      include 'types.f'

************************************************************************
*     Author J.M.Campbell, February 2000                               *
*                                                                      *
*     Supplemental to xzqqgg_v - just calculates the axial piece       *
*     with i1 and i4 swapped wrt that routine                          *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'lc.f'
      integer:: j,lh,h2,h3,hq,h(2:3)
      real(dp):: fac
      complex(dp):: m(2),mqqb_ax(2,2)
      complex(dp):: ml1_ax(2),ml2_ax(2)
      complex(dp):: a6treeg1,a64ax,a65ax
      integer,parameter::i1(2)=(/1,4/),i2(2)=(/2,3/),i3(2)=(/3,2/),
     &                   i4(2)=(/4,1/),i5(2)=(/6,5/),i6(2)=(/5,6/)
      character*9,parameter:: st1(2,2)=
     & reshape((/'q+g-g-qb-','q+g-g+qb-','q+g+g-qb-','q+g+g+qb-'/)
     & ,(/2,2/))
      character*9,parameter:: st3(2,2)=
     & reshape((/'q+qb-g-g-','q+qb-g-g+','q+qb-g+g-','q+qb-g+g+'/)
     & ,(/2,2/))

      fac=avegg*8._dp*gsq**2*esq**2*cf*xn**3*ason2pi
c--- no extra factor here since colour algebra is already done in (2.12)

      do hq=1,2
      do lh=1,2
      mqqb_ax(hq,lh)=0._dp

      if (colourchoice <= 1) then
      do h2=1,2
      do h3=1,2
        h(2)=h2
        h(3)=h3
        do j=1,2
        if (hq == 1) then
        m(j)=  a6treeg1(st1(3-h(i2(j)),3-h(i3(j))),
     &     i4(1),i2(j),i3(j),i1(1),i6(lh),i5(lh),zb,za)
c--- note: this symmetry relation (including minus sign) checked numerically
        ml1_ax(j)=-a64ax(st3(3-h(i2(j)),3-h(i3(j))),
     &     i4(1),i1(1),i2(j),i3(j),i6(lh),i5(lh),zb,za)
        ml2_ax(j)=-a65ax(st3(3-h(i2(j)),3-h(i3(j))),
     &     i4(1),i1(1),i2(j),i3(j),i6(lh),i5(lh),zb,za)
        else
        m(j)=  a6treeg1(st1(h(i2(j)),h(i3(j))),
     &     i4(1),i2(j),i3(j),i1(1),i5(lh),i6(lh),za,zb)
        ml1_ax(j)=a64ax(st3(h(i2(j)),h(i3(j))),
     &     i4(1),i1(1),i2(j),i3(j),i5(lh),i6(lh),za,zb)
        ml2_ax(j)=a65ax(st3(h(i2(j)),h(i3(j))),
     &     i4(1),i1(1),i2(j),i3(j),i5(lh),i6(lh),za,zb)
        endif
        enddo

      mqqb_ax(hq,lh)=mqqb_ax(hq,lh)+fac/xnsq*(
     &  conjg(m(1))*(
     &    (xn-2._dp/xn)*ml1_ax(1)-2._dp/xn*ml1_ax(2)+one/xn*ml2_ax(1))
     & +conjg(m(2))*(
     &    (xn-2._dp/xn)*ml1_ax(2)-2._dp/xn*ml1_ax(1)+one/xn*ml2_ax(2)))

      enddo
      enddo
      endif

      enddo
      enddo

      return
      end

