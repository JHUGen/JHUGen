      subroutine gen9dk_rap(r,p,wt9,*)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'decay1q2a.f'
      integer nu
      double precision r(mxdim)
      double precision wt9,p(mxpart,4),tp(4),tm(4),bp(4),bm(4),n(4),e(4)
      double precision bmg(4),g(4),ep(4),em(4),nn(4),nb(4),wp(4),wm(4)
      double precision wtttw,wttp,wttm,wtbmg,wtwp,wtwm,s3min,wt0
      parameter(wt0=1d0/twopi**5)
      
c--- alternate radiation between decay of top (=1) and antitop (=2) quarks;
c--- use additional (always uniform) variable to determine choice
      if (r(24) .lt. 0.5d0) then
        decay1q2a=1
      else
        decay1q2a=2
      endif

*     q(-p1) +qbar(-p2)=t(nu(p3)+e^+(p4)+b(p5))                        *
*                       +t~(b~(p6)+e^-(p7)+nu(p8)+g(p11))              *
*                       +W^+(nu(p9),mu^+(p10))                         *

      wt9=0d0
C---call gen4 that uses r(1)....r(10)
      call gen4(r,p,wtttw,*999)
      wtttw=(mt*twidth*pi)**2*wtttw
      do nu=1,4
      tp(nu)=p(5,nu)
      tm(nu)=p(6,nu)
      n(nu)=p(3,nu)
      e(nu)=p(4,nu)      
      enddo
      
      s3min=0d0
c      n2=0
c      n3=1
c      mass3=wmass      
c      width3=wwidth      

c--- radiation from top (not W)
      call phi1_2(r(11),r(12),r(13),r(14),tp,bmg,wp,wttp,*999)
      call phi3m(r(15),r(16),bmg,bm,g,mb,zip,wtbmg,*999)
      call phi3m0(r(17),r(18),wp,nn,ep,wtwp,*999)
c--- anti-top decay
      call phi1_2m(mb,r(19),r(20),r(21),s3min,tm,bp,wm,wttm,*999)
      call phi3m0(r(22),r(23),wm,em,nb,wtwm,*999)
      wt9=wt0*wtttw*wttp*wtbmg*wtwp*wttm*wtwm 
c--- multiply by a factor of two since we are including radiation
c--- from top and anti-top quarks at the same time
      wt9=wt9*2d0

      do nu=1,4
        if (decay1q2a .eq. 1) then
          p(3,nu)=nn(nu)
          p(4,nu)=ep(nu)
          p(5,nu)=bm(nu)
          p(6,nu)=bp(nu)
          p(7,nu)=em(nu)
          p(8,nu)=nb(nu)
        else
          p(3,nu)=nb(nu)
          p(4,nu)=em(nu)
          p(5,nu)=bp(nu)
          p(6,nu)=bm(nu)
          p(7,nu)=ep(nu)
          p(8,nu)=nn(nu)
        endif
        p(9,nu)=n(nu)
        p(10,nu)=e(nu)
        p(11,nu)=g(nu)
      enddo

      return

 999  wt9=0d0
      p(:,:)=0d0
      return 1
      
      end

