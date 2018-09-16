      subroutine qqb_QQbdk_g(p,msq)
      implicit none
      include 'types.f'

C***********************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2002.                                                     *
*     calculate the element squared                                    *
*     for the process                                                  *
*----My notation                                                       *
*     q(-p1) +qbar(-p2)=t(nu(p3)+e^+(p4)+b(p5))                        *
*                      +t~(b~(p6)+e^-(p7)+nu~(p8))+g(p9)               *
*                                                                      *
*     Only five diagrams included, leading to 2 on-shell top quarks    *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'plabel.f'
      include 'zprods_com.f'
      integer:: j,k,h1,h2,h3,nu
      real(dp):: t(4),r(4),
     & msq(-nf:nf,-nf:nf),p(mxpart,4),ps(mxpart,4)
      real(dp):: ttbqqbg_sq,fac,
     & wtgg,wtqqb,wtqbq,wtqg,wtqbarg,wtgq,wtgqbar
      complex(dp)::
     & ttbgggppp,ttbgggmpp,ttbgggpmp,ttbgggppm,
     & ttbgggmmm,ttbgggpmm,ttbgggmpm,ttbgggmmp,
     & a123(2,2,2),a132(2,2,2),a213(2,2,2),a231(2,2,2),
     & a312(2,2,2),a321(2,2,2),
     & a6sum,a3sum1a,a3sum1b,a3sum2a,a3sum2b,a3sum3a,a3sum3b
      real(dp):: p3Dp5,p6Dp8,rDp7,tDp4,s34,s78
c--- these definitions are used for gauge check only
c      complex(dp):: a,
c     & ttbgggppp_full,ttbgggmpp_full,ttbgggpmp_full,ttbgggppm_full,
c     & ttbgggmmm_full,ttbgggpmm_full,ttbgggmpm_full,ttbgggmmp_full,
c     & ttbqqbsqpp_full,ttbqqbsqpm_full,ttbqqbsqmp_full,ttbqqbsqmm_full,
c     & ttbqqbtqpp_full,ttbqqbtqpm_full,ttbqqbtqmp_full,ttbqqbtqmm_full,
c     & ttbqqbqqpp_full,ttbqqbqqpm_full,ttbqqbqqmp_full,ttbqqbqqmm_full,
c     & ttbqqbrqpp_full,ttbqqbrqpm_full,ttbqqbrqmp_full,ttbqqbrqmm_full

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      p3Dp5=p(3,4)*p(5,4)-p(3,3)*p(5,3)-p(3,2)*p(5,2)-p(3,1)*p(5,1)
      p6Dp8=p(6,4)*p(8,4)-p(6,3)*p(8,3)-p(6,2)*p(8,2)-p(6,1)*p(8,1)

      s34=2._dp*(p(3,4)*p(4,4)-p(3,3)*p(4,3)-p(3,2)*p(4,2)-p(3,1)*p(4,1))
      s78=2._dp*(p(7,4)*p(8,4)-p(7,3)*p(8,3)-p(7,2)*p(8,2)-p(7,1)*p(8,1))

c      we will have no further need for p3 and p5
c      we will have no further need for p6 and p8

      do nu=1,4
      t(nu)=p(3,nu)+p(4,nu)+p(5,nu)
      r(nu)=p(6,nu)+p(7,nu)+p(8,nu)
      enddo
      tDp4=t(4)*p(4,4)-t(3)*p(4,3)-t(2)*p(4,2)-t(1)*p(4,1)
      rDp7=r(4)*p(7,4)-r(3)*p(7,3)-r(2)*p(7,2)-r(1)*p(7,1)
      do nu=1,4
      ps(1,nu)=p(1,nu)
      ps(2,nu)=p(2,nu)
      ps(3,nu)=p(9,nu)
c---rescaled positron
      ps(4,nu)=0.5_dp*mt**2/tDp4*p(4,nu)
c---demassifyied top
      ps(5,nu)=t(nu)-ps(4,nu)
c---rescaled electron
      ps(7,nu)=0.5_dp*mt**2/rDp7*p(7,nu)
c---demassifyied antitop
      ps(6,nu)=r(nu)-ps(7,nu)
      enddo

c      call writeout(ps)
c      pause

      call spinoru(7,ps,za,zb)

c--- gauge check of the ggg pieces (performed on 15/8/08)
c      j1=2
c      j2=1
c      write(6,*) 'ppp'
c      do j3=3,8
c      a=ttbgggppp_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,abs(a)
c      enddo
c      write(6,*) 'mmm'
c      do j3=3,8
c      a=ttbgggmmm_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,abs(a)
c      enddo
c      write(6,*) 'mpp'
c      do j3=3,8
c      a=ttbgggmpp_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,abs(a)
c      enddo
c      write(6,*) 'pmm'
c      do j3=3,8
c      a=ttbgggpmm_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,abs(a)
c      enddo
c      write(6,*) 'mpm'
c      do j3=3,8
c      a=ttbgggmpm_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,abs(a)
c      enddo
c      write(6,*) 'mmp'
c      do j3=3,8
c      a=ttbgggmmp_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,abs(a)
c      enddo
c      write(6,*) 'pmp'
c      do j3=3,8
c      a=ttbgggpmp_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,abs(a)
c      enddo
c      write(6,*) 'ppm'
c      do j3=3,8
c      a=ttbgggppm_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,abs(a)
c      enddo
c      pause

c--- gauge check of the qqbg pieces (performed on 15/8/08)
c      write(6,*) 'pp'
c      do j=3,8
c      write(6,*) j,abs(ttbqqbsqpp_full(1,2,9,5,3,6,8,j)),
c     &           abs(ttbqqbtqpp_full(1,2,9,5,3,6,8,j)),
c     &           abs(ttbqqbqqpp_full(1,2,9,5,3,6,8,j)),
c     &           abs(ttbqqbrqpp_full(1,2,9,5,3,6,8,j))
c      enddo
c      write(6,*) 'pm'
c      do j=3,8
c      write(6,*) j,abs(ttbqqbsqpm_full(1,2,9,5,3,6,8,j)),
c     &           abs(ttbqqbtqpm_full(1,2,9,5,3,6,8,j)),
c     &           abs(ttbqqbqqpm_full(1,2,9,5,3,6,8,j)),
c     &           abs(ttbqqbrqpm_full(1,2,9,5,3,6,8,j))
c      enddo
c      write(6,*) 'mp'
c      do j=3,8
c      write(6,*) j,abs(ttbqqbsqmp_full(1,2,9,5,3,6,8,j)),
c     &           abs(ttbqqbtqmp_full(1,2,9,5,3,6,8,j)),
c     &           abs(ttbqqbqqmp_full(1,2,9,5,3,6,8,j)),
c     &           abs(ttbqqbrqmp_full(1,2,9,5,3,6,8,j))
c      enddo
c      write(6,*) 'mm'
c      do j=3,8
c      write(6,*) j,abs(ttbqqbsqmm_full(1,2,9,5,3,6,8,j)),
c     &           abs(ttbqqbtqmm_full(1,2,9,5,3,6,8,j)),
c     &           abs(ttbqqbqqmm_full(1,2,9,5,3,6,8,j)),
c     &           abs(ttbqqbrqmm_full(1,2,9,5,3,6,8,j))
c      enddo
c      pause

c--- according to gfortify.frm and ttbggg.frm, the relationship between
c--- the indices of the amplitudes calls and the real momenta is:
c---  i1 -> p1                                         = ps(1)
c---  i2 -> p2                                         = ps(2)
c---  i3 -> p9                                         = ps(9)
c---  i4 -> pta                                        = ps(5)
c---  i5 -> ptb  with pta a rescaled version of p4     = ps(3)
c---  i6 -> pua                                        = ps(6)
c---  i7 -> pub  with pub a rescaled version of p7     = ps(8)

      a123(2,2,2)=ttbgggppp(1,2,3,4,5,6,7)
      a123(2,2,1)=ttbgggppm(1,2,3,4,5,6,7)
      a123(2,1,2)=ttbgggpmp(1,2,3,4,5,6,7)
      a123(2,1,1)=ttbgggpmm(1,2,3,4,5,6,7)
      a123(1,2,2)=ttbgggmpp(1,2,3,4,5,6,7)
      a123(1,2,1)=ttbgggmpm(1,2,3,4,5,6,7)
      a123(1,1,2)=ttbgggmmp(1,2,3,4,5,6,7)
      a123(1,1,1)=ttbgggmmm(1,2,3,4,5,6,7)

      a132(2,2,2)=ttbgggppp(1,3,2,4,5,6,7)
      a132(2,2,1)=ttbgggppm(1,3,2,4,5,6,7)
      a132(2,1,2)=ttbgggpmp(1,3,2,4,5,6,7)
      a132(2,1,1)=ttbgggpmm(1,3,2,4,5,6,7)
      a132(1,2,2)=ttbgggmpp(1,3,2,4,5,6,7)
      a132(1,2,1)=ttbgggmpm(1,3,2,4,5,6,7)
      a132(1,1,2)=ttbgggmmp(1,3,2,4,5,6,7)
      a132(1,1,1)=ttbgggmmm(1,3,2,4,5,6,7)

      a213(2,2,2)=ttbgggppp(2,1,3,4,5,6,7)
      a213(2,2,1)=ttbgggppm(2,1,3,4,5,6,7)
      a213(2,1,2)=ttbgggpmp(2,1,3,4,5,6,7)
      a213(2,1,1)=ttbgggpmm(2,1,3,4,5,6,7)
      a213(1,2,2)=ttbgggmpp(2,1,3,4,5,6,7)
      a213(1,2,1)=ttbgggmpm(2,1,3,4,5,6,7)
      a213(1,1,2)=ttbgggmmp(2,1,3,4,5,6,7)
      a213(1,1,1)=ttbgggmmm(2,1,3,4,5,6,7)

      a231(2,2,2)=ttbgggppp(2,3,1,4,5,6,7)
      a231(2,2,1)=ttbgggppm(2,3,1,4,5,6,7)
      a231(2,1,2)=ttbgggpmp(2,3,1,4,5,6,7)
      a231(2,1,1)=ttbgggpmm(2,3,1,4,5,6,7)
      a231(1,2,2)=ttbgggmpp(2,3,1,4,5,6,7)
      a231(1,2,1)=ttbgggmpm(2,3,1,4,5,6,7)
      a231(1,1,2)=ttbgggmmp(2,3,1,4,5,6,7)
      a231(1,1,1)=ttbgggmmm(2,3,1,4,5,6,7)

      a312(2,2,2)=ttbgggppp(3,1,2,4,5,6,7)
      a312(2,2,1)=ttbgggppm(3,1,2,4,5,6,7)
      a312(2,1,2)=ttbgggpmp(3,1,2,4,5,6,7)
      a312(2,1,1)=ttbgggpmm(3,1,2,4,5,6,7)
      a312(1,2,2)=ttbgggmpp(3,1,2,4,5,6,7)
      a312(1,2,1)=ttbgggmpm(3,1,2,4,5,6,7)
      a312(1,1,2)=ttbgggmmp(3,1,2,4,5,6,7)
      a312(1,1,1)=ttbgggmmm(3,1,2,4,5,6,7)

      a321(2,2,2)=ttbgggppp(3,2,1,4,5,6,7)
      a321(2,2,1)=ttbgggppm(3,2,1,4,5,6,7)
      a321(2,1,2)=ttbgggpmp(3,2,1,4,5,6,7)
      a321(2,1,1)=ttbgggpmm(3,2,1,4,5,6,7)
      a321(1,2,2)=ttbgggmpp(3,2,1,4,5,6,7)
      a321(1,2,1)=ttbgggmpm(3,2,1,4,5,6,7)
      a321(1,1,2)=ttbgggmmp(3,2,1,4,5,6,7)
      a321(1,1,1)=ttbgggmmm(3,2,1,4,5,6,7)

      wtgg=0._dp
      do h1=1,2
      do h2=1,2
      do h3=1,2
c--- NB: make sure to permute helicity labels appropriately too
        a3sum3a=a123(h1,h2,h3)+a132(h1,h3,h2)+a312(h3,h1,h2)
        a3sum1a=a231(h1,h2,h3)+a213(h1,h3,h2)+a123(h3,h1,h2)
        a3sum2a=a312(h1,h2,h3)+a321(h1,h3,h2)+a231(h3,h1,h2)
        a3sum3b=a213(h2,h1,h3)+a231(h2,h3,h1)+a321(h3,h2,h1)
        a3sum1b=a321(h2,h1,h3)+a312(h2,h3,h1)+a132(h3,h2,h1)
        a3sum2b=a132(h2,h1,h3)+a123(h2,h3,h1)+a213(h3,h2,h1)
      a6sum=a3sum1a+a3sum1b
      wtgg=wtgg+xn**3*cf*(
     &   (abs(a123(h1,h2,h3))**2+abs(a132(h1,h3,h2))**2
     &   +abs(a213(h2,h1,h3))**2+abs(a231(h2,h3,h1))**2
     &   +abs(a312(h3,h1,h2))**2+abs(a321(h3,h2,h1))**2)
     &  -(abs(a3sum1a)**2+abs(a3sum2a)**2+abs(a3sum3a)**2
     &   +abs(a3sum1b)**2+abs(a3sum2b)**2+abs(a3sum3b)**2)/xn**2
     &  +(abs(a6sum)**2)*(xn**2+1._dp)/xn**4
     &                     )
      enddo
      enddo
      enddo

      wtqqb=ttbqqbg_sq(1,2,3,4,5,6,7)
      wtqbq=ttbqqbg_sq(2,1,3,4,5,6,7)
      wtqg=ttbqqbg_sq(1,3,2,4,5,6,7)
      wtgq=ttbqqbg_sq(2,3,1,4,5,6,7)
      wtqbarg=ttbqqbg_sq(3,1,2,4,5,6,7)
      wtgqbar=ttbqqbg_sq(3,2,1,4,5,6,7)

c--- overall factor, starting with couplings
      fac=2._dp*(gwsq/2._dp)**4*gsq**3
c--- include top and anti-top propagators
      fac=fac/(mt*twidth)**4
c--- include W decays
      fac=fac
     & *8._dp*p3Dp5/((s34-wmass**2)**2+(wmass*wwidth)**2)
     & *8._dp*p6Dp8/((s78-wmass**2)**2+(wmass*wwidth)**2)
c--- correct normalization for p4 and p7
      fac=fac
     & /(0.5_dp*mt**2/tDp4)
     & /(0.5_dp*mt**2/rDp7)
C--include factor for hadronic decays
      if (plabel(3)== 'pp') fac=2._dp*xn*fac
      if (plabel(7)== 'pp') fac=2._dp*xn*fac

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if     (j > 0) then
          msq(j,-j)=aveqq*fac*wtqqb
          msq(j, 0)=aveqg*fac*wtqg
          msq(0, j)=aveqg*fac*wtgq
      elseif (j == 0) then
          msq(j,j)=avegg*fac*wtgg
      elseif (j < 0) then
          msq(j,-j)=aveqq*fac*wtqbq
          msq(j, 0)=aveqg*fac*wtqbarg
          msq(0, j)=aveqg*fac*wtgqbar
      endif
      enddo
      return
      end



