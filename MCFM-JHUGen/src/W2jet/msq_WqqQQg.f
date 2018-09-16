      subroutine msq_WqqQQg(i1,i2,i3,i4,i5,i6,i7,MN,ofac)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2001.                                                     *
*     Return matrix elements squared as a function of f1,f2,f3,f4      *
*     summed over helicity, using the formulae of                      *
*     Nagy and Trocsanyi, PRD59 014020 (1999)                          *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: Qh,hq,hg,f1,f2,f3,f4,i1,i2,i3,i4,i5,i6,i7
      real(dp):: A(5,5,5,5),B(5,5,5,5),C(5,5,5,5),D(5,5,5,5),
     & E(5,5,5,5),F(5,5,5,5),G(5,5,5,5)
      real(dp):: MN(5,5,5,5),
     & M0(5,5,5,5),Mx(5,5,5,5),My(5,5,5,5),
     & Mz(5,5,5,5),Mxx(5,5,5,5),Mxy(5,5,5,5),temp,ofac
      real(dp):: x,y,z
      parameter(x=xn/cf,y=half/cf,z=0.25_dp*(xn**2-two)/xn/cf**2)
      complex(dp)::
     &               mb1_1234(5,5,5,5,2,2,2),mb2_1234(5,5,5,5,2,2,2),
     &               mb1_3412(5,5,5,5,2,2,2),mb2_3412(5,5,5,5,2,2,2),
     &               mb1_3214(5,5,5,5,2,2,2),mb2_3214(5,5,5,5,2,2,2),
     &               mb1_1432(5,5,5,5,2,2,2),mb2_1432(5,5,5,5,2,2,2)
      integer:: d1234,d3412,d3214,d1432,hc,hd

C---set everything to zero
      do f1=1,5
      do f2=1,5
      do f3=1,5
      do f4=1,5
      MN(f1,f2,f3,f4)=zip
      A(f1,f2,f3,f4)=zip
      B(f1,f2,f3,f4)=zip
      C(f1,f2,f3,f4)=zip
      D(f1,f2,f3,f4)=zip
      E(f1,f2,f3,f4)=zip
      F(f1,f2,f3,f4)=zip
      G(f1,f2,f3,f4)=zip
      enddo
      enddo
      enddo
      enddo
      if (s(i6,i7) < 4*mbsq) return

C---exclude the photon pole, 4*mbsq choosen as a scale approx above upsilon

c---mb1_1234 etc, has 6 indices each with possible values 1 or 2
C---corresponding to f1,f3,hq,Qh,hg,lh
      call wmakemb(i1,i2,i3,i4,i5,i6,i7,mb1_1234,mb2_1234)
      call wmakemb(i3,i4,i1,i2,i5,i6,i7,mb1_3412,mb2_3412)
      call wmakemb(i3,i2,i1,i4,i5,i6,i7,mb1_3214,mb2_3214)
      call wmakemb(i1,i4,i3,i2,i5,i6,i7,mb1_1432,mb2_1432)

      do f1=1,5
      do f2=1,5
      do f3=1,5
      do f4=1,5

      temp=0._dp

      do hq=1,2
      do Qh=1,2
      do hg=1,2

      do hc=1,2
      do hd=1,2

      d1234=0
      d3412=0
      d3214=0
      d1432=0
      if ((hq == hc) .and. (Qh == hd)) then
        d1234=1
        d3412=1
      endif
      if ((hq == hd) .and. (Qh == hc)) then
        d3214=1
        d1432=1
      endif

      A(f1,f2,f3,f4)=A(f1,f2,f3,f4)
     & +abs(mb1_1234(f1,f2,f3,f4,hq,Qh,hg)*d1234)**2
     & +abs(mb2_1234(f1,f2,f3,f4,hq,Qh,hg)*d1234)**2
     & +abs(mb1_3412(f3,f4,f1,f2,Qh,hq,hg)*d3412)**2
     & +abs(mb2_3412(f3,f4,f1,f2,Qh,hq,hg)*d3412)**2

      D(f1,f2,f3,f4)=D(f1,f2,f3,f4)+two*real(
     .+mb1_1234(f1,f2,f3,f4,hq,Qh,hg)*d1234
     .*conjg(mb2_1234(f1,f2,f3,f4,hq,Qh,hg)*d1234)
     .+mb1_3412(f3,f4,f1,f2,Qh,hq,hg)*d3412
     .*conjg(mb2_3412(f3,f4,f1,f2,Qh,hq,hg)*d3412))

      A(f1,f2,f3,f4)=A(f1,f2,f3,f4)
     & +abs(mb1_3214(f3,f4,f1,f2,Qh,hq,hg)*d3214)**2
     & +abs(mb2_3214(f3,f4,f1,f2,Qh,hq,hg)*d3214)**2
     & +abs(mb1_1432(f1,f2,f3,f4,hq,Qh,hg)*d1432)**2
     & +abs(mb2_1432(f1,f2,f3,f4,hq,Qh,hg)*d1432)**2

      D(f1,f2,f3,f4)=D(f1,f2,f3,f4)+two*real(
     .+mb1_3214(f3,f4,f1,f2,Qh,hq,hg)*d3214
     .*conjg(mb2_3214(f3,f4,f1,f2,Qh,hq,hg)*d3214)
     .+mb1_1432(f1,f2,f3,f4,hq,Qh,hg)*d1432
     .*conjg(mb2_1432(f1,f2,f3,f4,hq,Qh,hg)*d1432))

      B(f1,f2,f3,f4)=B(f1,f2,f3,f4)-two*real(
     .+mb1_1234(f1,f2,f3,f4,hq,Qh,hg)*d1234
     .*conjg(mb1_1432(f1,f2,f3,f4,hq,Qh,hg)*d1432)
     .+mb2_1234(f1,f2,f3,f4,hq,Qh,hg)*d1234
     .*conjg(mb2_3214(f3,f4,f1,f2,Qh,hq,hg)*d3214)
     .+mb1_3412(f3,f4,f1,f2,Qh,hq,hg)*d3412
     .*conjg(mb1_3214(f3,f4,f1,f2,Qh,hq,hg)*d3214)
     .+mb2_3412(f3,f4,f1,f2,Qh,hq,hg)*d3412
     .*conjg(mb2_1432(f1,f2,f3,f4,hq,Qh,hg)*d1432))

      C(f1,f2,f3,f4)=C(f1,f2,f3,f4)-two*real(
     .+mb1_1234(f1,f2,f3,f4,hq,Qh,hg)*d1234
     .*conjg(mb1_3214(f3,f4,f1,f2,Qh,hq,hg)*d3214)
     .+mb2_1234(f1,f2,f3,f4,hq,Qh,hg)*d1234
     .*conjg(mb2_1432(f1,f2,f3,f4,hq,Qh,hg)*d1432)
     .+mb1_3412(f3,f4,f1,f2,Qh,hq,hg)*d3412
     .*conjg(mb1_1432(f1,f2,f3,f4,hq,Qh,hg)*d1432)
     .+mb2_3412(f3,f4,f1,f2,Qh,hq,hg)*d3412
     .*conjg(mb2_3214(f3,f4,f1,f2,Qh,hq,hg)*d3214))

      E(f1,f2,f3,f4)=E(f1,f2,f3,f4)-two*real(
     .+(mb1_1234(f1,f2,f3,f4,hq,Qh,hg)*d1234
     & +mb1_3412(f3,f4,f1,f2,Qh,hq,hg)*d3412)
     .*conjg(mb2_3214(f3,f4,f1,f2,Qh,hq,hg)*d3214
     .+mb2_1432(f1,f2,f3,f4,hq,Qh,hg)*d1432)
     .+(mb2_1234(f1,f2,f3,f4,hq,Qh,hg)*d1234
     & +mb2_3412(f3,f4,f1,f2,Qh,hq,hg)*d3412)
     .*conjg(mb1_3214(f3,f4,f1,f2,Qh,hq,hg)*d3214
     & +mb1_1432(f1,f2,f3,f4,hq,Qh,hg)*d1432))

       F(f1,f2,f3,f4)=F(f1,f2,f3,f4)+two*real(
     & +mb1_3214(f3,f4,f1,f2,Qh,hq,hg)*d3214
     & *conjg(mb1_1432(f1,f2,f3,f4,hq,Qh,hg)*d1432)
     & +mb2_3214(f3,f4,f1,f2,Qh,hq,hg)*d3214
     & *conjg(mb2_1432(f1,f2,f3,f4,hq,Qh,hg)*d1432))

       G(f1,f2,f3,f4)=G(f1,f2,f3,f4)+two*real(
     & +mb1_3214(f3,f4,f1,f2,Qh,hq,hg)*d3214
     & *conjg(mb2_1432(f1,f2,f3,f4,hq,Qh,hg)*d1432)
     & +mb2_3214(f3,f4,f1,f2,Qh,hq,hg)*d3214
     & *conjg(mb1_1432(f1,f2,f3,f4,hq,Qh,hg)*d1432))

      F(f1,f2,f3,f4)=F(f1,f2,f3,f4)+two*real(
     & +mb1_1234(f1,f2,f3,f4,hq,Qh,hg)*d1234
     & *conjg(mb1_3412(f3,f4,f1,f2,Qh,hq,hg)*d3412)
     & +mb2_1234(f1,f2,f3,f4,hq,Qh,hg)*d1234
     & *conjg(mb2_3412(f3,f4,f1,f2,Qh,hq,hg)*d3412))

      G(f1,f2,f3,f4)=G(f1,f2,f3,f4)+two*real(
     & +mb1_1234(f1,f2,f3,f4,hq,Qh,hg)*d1234
     & *conjg(mb2_3412(f3,f4,f1,f2,Qh,hq,hg)*d3412)
     & +mb2_1234(f1,f2,f3,f4,hq,Qh,hg)*d1234
     & *conjg(mb1_3412(f3,f4,f1,f2,Qh,hq,hg)*d3412))

      M0(f1,f2,f3,f4)=
     & B(f1,f2,f3,f4)+C(f1,f2,f3,f4)+E(f1,f2,f3,f4)
      Mx(f1,f2,f3,f4)=-0.5_dp
     & *(3._dp*C(f1,f2,f3,f4)+2._dp*E(f1,f2,f3,f4)+B(f1,f2,f3,f4))
      My(f1,f2,f3,f4)=A(f1,f2,f3,f4)+D(f1,f2,f3,f4)
      Mz(f1,f2,f3,f4)=F(f1,f2,f3,f4)+G(f1,f2,f3,f4)
      Mxx(f1,f2,f3,f4)=0.25_dp*(2._dp*C(f1,f2,f3,f4)+E(f1,f2,f3,f4))
      Mxy(f1,f2,f3,f4)=-0.5_dp*(F(f1,f2,f3,f4)+D(f1,f2,f3,f4))

      MN(f1,f2,f3,f4)=
     & ofac*CF**3*xn*(M0(f1,f2,f3,f4)+x*Mx(f1,f2,f3,f4)
     & +y*My(f1,f2,f3,f4)+z*Mz(f1,f2,f3,f4)
     & +x**2*Mxx(f1,f2,f3,f4)+x*y*Mxy(f1,f2,f3,f4))

      temp=MN(f1,f2,f3,f4)

      enddo
      enddo

      enddo
      enddo
      enddo

      M0(f1,f2,f3,f4)=
     & B(f1,f2,f3,f4)+C(f1,f2,f3,f4)+E(f1,f2,f3,f4)
      Mx(f1,f2,f3,f4)=-0.5_dp
     & *(3._dp*C(f1,f2,f3,f4)+2._dp*E(f1,f2,f3,f4)+B(f1,f2,f3,f4))
      My(f1,f2,f3,f4)=A(f1,f2,f3,f4)+D(f1,f2,f3,f4)
      Mz(f1,f2,f3,f4)=F(f1,f2,f3,f4)+G(f1,f2,f3,f4)
      Mxx(f1,f2,f3,f4)=0.25_dp*(2._dp*C(f1,f2,f3,f4)+E(f1,f2,f3,f4))
      Mxy(f1,f2,f3,f4)=-0.5_dp*(F(f1,f2,f3,f4)+D(f1,f2,f3,f4))

      MN(f1,f2,f3,f4)=
     & CF**3*xn*(M0(f1,f2,f3,f4)+x*Mx(f1,f2,f3,f4)
     & +y*My(f1,f2,f3,f4)+z*Mz(f1,f2,f3,f4)
     & +x**2*Mxx(f1,f2,f3,f4)+x*y*Mxy(f1,f2,f3,f4))

      enddo
      enddo
      enddo
      enddo

      return
      end

