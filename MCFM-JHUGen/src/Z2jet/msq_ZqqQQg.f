      subroutine msq_ZqqQQg(i1,i2,i3,i4,i5,i6,i7,MN,MI)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2001.                                                     *
*     Return matrix elements squared as a function of f1 and f2        *
*     summed over helicity, using the formulae of                      *
*     Nagy and Trocsnayi, PRD59 014020 (1999)                          *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      integer:: Qh,hq,hg,lh,f1,f3,i1,i2,i3,i4,i5,i6,i7,j
      real(dp):: A(2,2,2),B(2,2,2),C(2,2,2),D(2,2,2),E(2,2,2),
     & F(2,2,2),G(2,2,2)
      real(dp):: MI(2),MN(2,2),M0(2,2,2),Mx(2,2,2),My(2,2,2),
     & Mz(2,2,2),Mxx(2,2,2),Mxy(2,2,2)
      real(dp):: x,y,z
      parameter(x=xn/cf,y=half/cf,z=0.25_dp*(xn**2-two)/xn/cf**2)
      complex(dp)::
     &               mb1_1234(2,2,2,2,2,2),mb2_1234(2,2,2,2,2,2),
     &               mb1_3412(2,2,2,2,2,2),mb2_3412(2,2,2,2,2,2),
     &               mb1_3214(2,2,2,2,2,2),mb2_3214(2,2,2,2,2,2),
     &               mb1_1432(2,2,2,2,2,2),mb2_1432(2,2,2,2,2,2)


C---set everything to zero

      do f1=1,2
      MI(f1)=zip
      do f3=1,2
      MN(f1,f3)=zip

      do j=1,2
      A(j,f1,f3)=zip
      B(j,f1,f3)=zip
      C(j,f1,f3)=zip
      D(j,f1,f3)=zip
      E(j,f1,f3)=zip
      F(j,f1,f3)=zip
      G(j,f1,f3)=zip
      enddo
      enddo
      enddo

c---mb1_1234 etc, has 6 indices each with possible values 1 or 2
C---corresponding to f1,f3,hq,Qh,hg,lh
      call makemb(i1,i2,i3,i4,i5,i6,i7,mb1_1234,mb2_1234)
      call makemb(i3,i2,i1,i4,i5,i6,i7,mb1_3214,mb2_3214)
      call makemb(i3,i4,i1,i2,i5,i6,i7,mb1_3412,mb2_3412)
      call makemb(i1,i4,i3,i2,i5,i6,i7,mb1_1432,mb2_1432)

      do f1=1,2
      do f3=1,2

      do hq=1,2
      do Qh=1,2
      do hg=1,2
      do lh=1,2
      A(1,f1,f3)=A(1,f1,f3)
     & +abs(mb1_1234(f1,f3,hq,Qh,hg,lh))**2
     & +abs(mb2_1234(f1,f3,hq,Qh,hg,lh))**2
     & +abs(mb1_3412(f3,f1,Qh,hq,hg,lh))**2
     & +abs(mb2_3412(f3,f1,Qh,hq,hg,lh))**2
      D(1,f1,f3)=D(1,f1,f3)+two*real(
     .+mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_1234(f1,f3,hq,Qh,hg,lh))
     .+mb1_3412(f3,f1,Qh,hq,hg,lh)*conjg(mb2_3412(f3,f1,Qh,hq,hg,lh)))
      B(1,f1,f3)=zip
      C(1,f1,f3)=zip
      E(1,f1,f3)=zip
      F(1,f1,f3)=F(1,f1,f3)+two*real(
     & +mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb1_3412(f3,f1,Qh,hq,hg,lh))
     & +mb2_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_3412(f3,f1,Qh,hq,hg,lh)))
      G(1,f1,f3)=G(1,f1,f3)+two*real(
     & +mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_3412(f3,f1,Qh,hq,hg,lh))
     & +mb2_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb1_3412(f3,f1,Qh,hq,hg,lh)))


      A(2,f1,f3)=A(2,f1,f3)
     & +abs(mb1_1234(f1,f3,hq,Qh,hg,lh))**2
     & +abs(mb2_1234(f1,f3,hq,Qh,hg,lh))**2
     & +abs(mb1_3214(f3,f1,Qh,hq,hg,lh))**2
     & +abs(mb2_3214(f3,f1,Qh,hq,hg,lh))**2
     & +abs(mb1_1432(f1,f3,hq,Qh,hg,lh))**2
     & +abs(mb2_1432(f1,f3,hq,Qh,hg,lh))**2
     & +abs(mb1_3412(f3,f1,Qh,hq,hg,lh))**2
     & +abs(mb2_3412(f3,f1,Qh,hq,hg,lh))**2
      D(2,f1,f3)=D(2,f1,f3)+two*real(
     .+mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_1234(f1,f3,hq,Qh,hg,lh))
     .+mb1_3214(f3,f1,Qh,hq,hg,lh)*conjg(mb2_3214(f3,f1,Qh,hq,hg,lh))
     .+mb1_1432(f1,f3,hq,Qh,hg,lh)*conjg(mb2_1432(f1,f3,hq,Qh,hg,lh))
     .+mb1_3412(f3,f1,Qh,hq,hg,lh)*conjg(mb2_3412(f3,f1,Qh,hq,hg,lh)))
      if (hq == Qh) then
      B(2,f1,f3)=B(2,f1,f3)-two*real(
     .+mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb1_1432(f1,f3,hq,Qh,hg,lh))
     .+mb2_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_3214(f3,f1,Qh,hq,hg,lh))
     .+mb1_3412(f3,f1,Qh,hq,hg,lh)*conjg(mb1_3214(f3,f1,Qh,hq,hg,lh))
     .+mb2_3412(f3,f1,Qh,hq,hg,lh)*conjg(mb2_1432(f1,f3,hq,Qh,hg,lh)))

      C(2,f1,f3)=C(2,f1,f3)-two*real(
     .+mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb1_3214(f3,f1,Qh,hq,hg,lh))
     .+mb2_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_1432(f1,f3,hq,Qh,hg,lh))
     .+mb1_3412(f3,f1,Qh,hq,hg,lh)*conjg(mb1_1432(f1,f3,hq,Qh,hg,lh))
     .+mb2_3412(f3,f1,Qh,hq,hg,lh)*conjg(mb2_3214(f3,f1,Qh,hq,hg,lh)))

      E(2,f1,f3)=E(2,f1,f3)-two*real(
     .+(mb1_1234(f1,f3,hq,Qh,hg,lh)+mb1_3412(f3,f1,Qh,hq,hg,lh))
     .*conjg(mb2_3214(f3,f1,Qh,hq,hg,lh)+mb2_1432(f1,f3,hq,Qh,hg,lh))
     .+(mb2_1234(f1,f3,hq,Qh,hg,lh)+mb2_3412(f3,f1,Qh,hq,hg,lh))
     .*conjg(mb1_3214(f3,f1,Qh,hq,hg,lh)+mb1_1432(f1,f3,hq,Qh,hg,lh)))
      endif
      F(2,f1,f3)=F(2,f1,f3)+two*real(
     & +mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb1_3412(f3,f1,Qh,hq,hg,lh))
     & +mb1_3214(f3,f1,Qh,hq,hg,lh)*conjg(mb1_1432(f1,f3,hq,Qh,hg,lh))
     & +mb2_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_3412(f3,f1,Qh,hq,hg,lh))
     & +mb2_3214(f3,f1,Qh,hq,hg,lh)*conjg(mb2_1432(f1,f3,hq,Qh,hg,lh)))

      G(2,f1,f3)=G(2,f1,f3)+two*real(
     & +mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_3412(f3,f1,Qh,hq,hg,lh))
     & +mb1_3214(f3,f1,Qh,hq,hg,lh)*conjg(mb2_1432(f1,f3,hq,Qh,hg,lh))
     & +mb2_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb1_3412(f3,f1,Qh,hq,hg,lh))
     & +mb2_3214(f3,f1,Qh,hq,hg,lh)*conjg(mb1_1432(f1,f3,hq,Qh,hg,lh)))

      enddo
      enddo
      enddo
      enddo

      do j=1,2
      M0(j,f1,f3)=B(j,f1,f3)+C(j,f1,f3)+E(j,f1,f3)
      Mx(j,f1,f3)=-0.5_dp*(3._dp*C(j,f1,f3)+2._dp*E(j,f1,f3)+B(j,f1,f3))
      My(j,f1,f3)=A(j,f1,f3)+D(j,f1,f3)
      Mz(j,f1,f3)=F(j,f1,f3)+G(j,f1,f3)
      Mxx(j,f1,f3)=0.25_dp*(2._dp*C(j,f1,f3)+E(j,f1,f3))
      Mxy(j,f1,f3)=-0.5_dp*(F(j,f1,f3)+D(j,f1,f3))
      enddo
      MN(f1,f3)=CF**3*xn*(M0(1,f1,f3)+x*Mx(1,f1,f3)+y*My(1,f1,f3)
     & +z*Mz(1,f1,f3)+x**2*Mxx(1,f1,f3)+x*y*Mxy(1,f1,f3))
      if ((f1==f3)) MI(f1)=
     & CF**3*xn*(M0(2,f1,f3)+x*Mx(2,f1,f3)+y*My(2,f1,f3)
     & +z*Mz(2,f1,f3)+x**2*Mxx(2,f1,f3)+x*y*Mxy(2,f1,f3))
      enddo
      enddo

      return
      end

