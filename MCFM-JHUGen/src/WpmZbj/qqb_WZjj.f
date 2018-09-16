      subroutine qqb_WZjj(p,msq)
      implicit none
      include 'types.f'
C-----Matrix squared for
C-----p(p1)+p(p2)-->W(l(p3)+a(p4))+Z(l(p5)+a(p6)+p(p7)+p(p8)
C-----Author: R.K.Ellis February 2013

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'nwz.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
     & ,WZbbmsq,WZccmsq,WZggmsq,WZddidmsq,WZuuidmsq,
     & ms(-nf:nf,-nf:nf,-nf:nf,-nf:nf)

      integer:: i,j,k,l,ba,ac,sa,ua,da,gg,dq,uq,sq,cq,bq
      parameter(ba=-5,ac=-4,sa=-3,ua=-2,da=-1,
     &     gg=0,bq=+5,cq=+4,sq=+3,uq=+2,dq=+1)

      msq(:,:)=0._dp
      ms(:,:,:,:)=0._dp
c      open(unit=66,file='madpoints.f',status='unknown')
c      call writeformadgraph(p)
c      write (*,*)
c      write (*,*) " Phase space point:"
c      write (*,*)
c      write (*,*) "-----------------------------------------------------------------------------"
c      write (*,*)  "n        E             px             py              pz               m "
c      do i=1,8
c      if (i < 3)
c     &    write (*,'(i2,1x,5e15.7)') i, -P(i,4),-P(i,1),-P(i,2),-P(i,3),
c     & p(i,4)**2-p(i,1)**2-p(i,2)**2-p(i,3)**2
c      if (i >= 3)
c     &    write (*,'(i2,1x,5e15.7)') i, P(i,4),P(i,1),P(i,2),P(i,3),
cc     & sqrt(abs(DOT(p(0,i),p(0,i))))
c     & p(i,4)**2-p(i,1)**2-p(i,2)**2-p(i,3)**2
c      enddo
c      write (*,*) "-----------------------------------------------------------------------------"
c      write (*,*) "-----------------------------------------------------------------------------"

      call spinoru(8,p,za,zb)

      if (nwz == 1) then

C***********************************************************
C-----Identity tables antiquark-quark
c--- d~  d > d  u~ and s  s~ > s  c~
      ms(da,dq,dq,ua)=WZddidmsq(8,7,3,4,5,6,1,2)
      ms(sa,sq,sq,ac)=ms(da,dq,dq,ua)
c--- u~  u > d  u~ and c~  c > s  c~
      ms(ua,uq,dq,ua)=WZuuidmsq(2,7,3,4,5,6,1,8)
      ms(ac,cq,sq,ac)=ms(ua,uq,dq,ua)
c--- d~  u > d  d~ and s~  c > s  s~
      ms(da,uq,dq,da)=WZddidmsq(2,1,3,4,5,6,7,8)
      ms(sa,cq,sq,sa)=ms(da,uq,dq,da)
c--- s~  s > d  u~ and d~  d > s  c~ and b~  b > d  u~ and b~  b > s  c~
      ms(da,dq,sq,ac)=WZbbmsq(8,7,3,4,5,6,1,2)
      ms(sa,sq,dq,ua)=ms(da,dq,sq,ac)
      ms(ba,bq,dq,ua)=ms(da,dq,sq,ac)
      ms(ba,bq,sq,ac)=ms(da,dq,sq,ac)
c--- d~  u > s  s~ and s~  c > d  d~ and s~  c > b  b~ and d~  u > b  b~
      ms(da,uq,sq,sa)=WZbbmsq(2,1,3,4,5,6,7,8)
      ms(sa,cq,dq,da)=ms(da,uq,sq,sa)
      ms(sa,cq,bq,ba)=ms(da,uq,sq,sa)
      ms(da,uq,bq,ba)=ms(da,uq,sq,sa)
c--- d~  u > u  u~ and s~  c > c  c~
      ms(da,uq,uq,ua)=WZuuidmsq(2,1,3,4,5,6,7,8)
      ms(sa,cq,cq,ac)=ms(da,uq,uq,ua)
c--- s~  u  > d  s~
      ms(sa,uq,dq,sa)=WZbbmsq(2,7,3,4,5,6,1,8)
c--- b~  u   > d  b~
      ms(ba,uq,dq,ba)=ms(sa,uq,dq,sa)
c--- d~  c  > s  d~
      ms(da,cq,sq,da)=ms(sa,uq,dq,sa)
c--- s~  d > d  c~
      ms(sa,dq,dq,ac)=WZbbmsq(8,1,3,4,5,6,7,2)
c--- d~  s  > s  u~
      ms(da,sq,sq,ua)=ms(sa,dq,dq,ac)
c--- d~  b  > b  u~
      ms(da,bq,bq,ua)=ms(sa,dq,dq,ac)
c--- s~  b  > b  c~
      ms(sa,bq,bq,ac)=ms(sa,dq,dq,ac)
c--- u~  u > s  c~ and c~  c > d  u~
      ms(ua,uq,sq,ac)=WZccmsq(8,7,3,4,5,6,1,2)
      ms(ac,cq,dq,ua)=ms(ua,uq,sq,ac)
c--- d~  u  > c  c~
      ms(da,uq,cq,ac)=WZccmsq(2,1,3,4,5,6,7,8)
c--- s~  c  > u  u~
      ms(sa,cq,uq,ua)=ms(da,uq,cq,ac)
c--- c~  u  > d  c~
      ms(ac,uq,dq,ac)=WZccmsq(2,7,3,4,5,6,1,8)
c--- u~  c  > s  u~
      ms(ua,cq,sq,ua)=ms(ac,uq,dq,ac)
c--- s~ u  > u  c~
      ms(sa,uq,uq,ac)=WZccmsq(8,1,3,4,5,6,7,2)
c--- d~  c  > c  u~
      ms(da,cq,cq,ua)=ms(sa,uq,uq,ac)
c--- b~  c > s  b~
      ms(ba,cq,sq,ba)=ms(sa,uq,dq,sa)

C***********************************************************
C-----Identity tables quark-antiquark
c--- d  d~ > d  u~ and s  s~ > s  c~
      ms(dq,da,dq,ua)=WZddidmsq(8,7,3,4,5,6,2,1)
      ms(sq,sa,sq,ac)=ms(dq,da,dq,ua)
c--- u  d~ > d  d~ and c  s~ > s  s~
      ms(uq,da,dq,da)=WZddidmsq(1,2,3,4,5,6,7,8)
      ms(cq,sa,sq,sa)=ms(uq,da,dq,da)
c--- u  d~ > u  u~ and c  s~ > c  c~
      ms(uq,da,uq,ua)=WZuuidmsq(1,2,3,4,5,6,7,8)
      ms(cq,sa,cq,ac)=ms(uq,da,uq,ua)
c--- u  u~ > d  u~ and c  c~ > s  c~
      ms(uq,ua,dq,ua)=WZuuidmsq(1,7,3,4,5,6,2,8)
      ms(cq,ac,sq,ac)=ms(uq,ua,dq,ua)
c--- s  s~ > d  u~ and d  d~ > s  c~ and b  b~ > d  u~ and b  b~ > s  c~
      ms(dq,da,sq,ac)=WZbbmsq(8,7,3,4,5,6,2,1)
      ms(sq,sa,dq,ua)=ms(dq,da,sq,ac)
      ms(bq,ba,dq,ua)=ms(dq,da,sq,ac)
      ms(bq,ba,sq,ac)=ms(dq,da,sq,ac)
c--- u  d~ > s  s~ and c  s~ > d  d~ and c  s~ > b  b~ and u  d~ > b  b~
      ms(uq,da,sq,sa)=WZbbmsq(1,2,3,4,5,6,7,8)
      ms(cq,sa,dq,da)=ms(uq,da,sq,sa)
      ms(cq,sa,bq,ba)=ms(uq,da,sq,sa)
      ms(uq,da,bq,ba)=ms(uq,da,sq,sa)
c--- u  s~ > d  s~
      ms(uq,sa,dq,sa)=WZbbmsq(1,7,3,4,5,6,2,8)
c--- u  b~ > d  b~
      ms(uq,ba,dq,ba)=ms(uq,sa,dq,sa)
c--- c  d~ > s  d~
      ms(cq,da,sq,da)=ms(uq,sa,dq,sa)
c--- d  s~ > d  c~
      ms(dq,sa,dq,ac)=WZbbmsq(8,2,3,4,5,6,7,1)
c--- s  d~ > s  u~
      ms(sq,da,sq,ua)=ms(dq,sa,dq,ac)
c--- b  d~ > b  u~
      ms(bq,da,bq,ua)=ms(dq,sa,dq,ac)
c--- b  s~ > b  c~
      ms(bq,sa,bq,ac)=ms(dq,sa,dq,ac)
c--- u  u~ > s  c~ and c  c~ > d  u~
      ms(uq,ua,sq,ac)=WZccmsq(8,7,3,4,5,6,2,1)
      ms(cq,ac,dq,ua)=ms(uq,ua,sq,ac)
c--- u  d~ > c  c~
      ms(uq,da,cq,ac)=WZccmsq(1,2,3,4,5,6,7,8)
c--- c  s~ > u  u~
      ms(cq,sa,uq,ua)=ms(uq,da,cq,ac)
c--- u  c~ > d  c~
      ms(uq,ac,dq,ac)=WZccmsq(1,7,3,4,5,6,2,8)
c--- c  u~ > s  u~
      ms(cq,ua,sq,ua)=ms(uq,ac,dq,ac)
c--- u  s~ > u  c~
      ms(uq,sa,uq,ac)=WZccmsq(8,2,3,4,5,6,7,1)
c--- c  d~ > c  u~
      ms(cq,da,cq,ua)=ms(uq,sa,uq,ac)
c--- c  b~ > s  b~
      ms(cq,ba,sq,ba)=ms(uq,sa,dq,sa)
C***********************************************************
C-----Identity tables antiquark-antiquark
c--- u~ d~ > u~ u~
      ms(ua,da,ua,ua)=half*WZuuidmsq(7,2,3,4,5,6,1,8)
c--- c~ s~ > c~ c~
      ms(ac,sa,ac,ac)=ms(ua,da,ua,ua)
c--- d~ u~ > u~ u~
      ms(da,ua,ua,ua)=half*WZuuidmsq(7,1,3,4,5,6,2,8)
c--- s~ c~ > c~ c~
      ms(sa,ac,ac,ac)=ms(da,ua,ua,ua)
c--- d~ d~ > d~ u~
      ms(da,da,da,ua)=WZddidmsq(8,1,3,4,5,6,2,7)
c--- s~ s~ > s~ c~
      ms(sa,sa,sa,ac)=ms(da,da,da,ua)
c--- d~ c~ > u~ c~
      ms(da,ac,ua,ac)=WZccmsq(7,1,3,4,5,6,2,8)
c--- s~ u~ > c~ u~
      ms(sa,ua,ac,ua)=ms(da,ac,ua,ac)
c--- u~ s~ > c~  u~
      ms(ua,sa,ac,ua)=WZccmsq(7,2,3,4,5,6,1,8)
c--- c~ d~ > u~ c~
      ms(ac,da,ua,ac)=ms(ua,sa,ac,ua)
c--- s~ d~ u~ s~ and s~ d~ > c~ d~
c--- d~ s~ > u~ s~ and d~ s~ > c~ d~
      ms(sa,da,ua,sa)=WZbbmsq(7,2,3,4,5,6,1,8)
      ms(sa,da,ac,da)=WZbbmsq(7,1,3,4,5,6,2,8)
      ms(da,sa,ua,sa)=ms(sa,da,ac,da)
      ms(da,sa,ac,da)=ms(sa,da,ua,sa)
c--- d~ b~ > u~ b~
      ms(da,ba,ua,ba)=ms(da,sa,ua,sa)
c--- s~ b~ > c~ b~
      ms(sa,ba,ac,ba)=ms(da,sa,ua,sa)
c--- b~ d~ > u~ b~
      ms(ba,da,ua,ba)=ms(sa,da,ua,sa)
c--- b~ s~ > c~ b~
      ms(ba,sa,ac,ba)=ms(sa,da,ua,sa)
C***********************************************************
C-----Identity tables quark-quark
c--- u  d  > d  d
      ms(uq,dq,dq,dq)=half*WZddidmsq(1,7,3,4,5,6,8,2)
c--- c  s  > s  s
      ms(cq,sq,sq,sq)=ms(uq,dq,dq,dq)
c--- d  u  > d  d
      ms(dq,uq,dq,dq)=half*WZddidmsq(2,7,3,4,5,6,8,1)
c--- s  c  > s  s
      ms(sq,cq,sq,sq)=ms(dq,uq,dq,dq)
c--- u  u  > d  u
      ms(uq,uq,dq,uq)=WZuuidmsq(1,7,3,4,5,6,8,2)
c--- c  c  > s  c
      ms(cq,cq,sq,cq)=ms(uq,uq,dq,uq)
c--- u  b  > d  b
      ms(uq,bq,dq,bq)=WZbbmsq(1,7,3,4,5,6,8,2)
c--- u  s  > d  s
      ms(uq,sq,dq,sq)=ms(uq,bq,dq,bq)
c--- c  b  > s  b
      ms(cq,bq,sq,bq)=ms(uq,bq,dq,bq)
c--- c  d  > s  d
      ms(cq,dq,sq,dq)=ms(uq,bq,dq,bq)
c--- b  u  > d  b
      ms(bq,uq,dq,bq)=WZbbmsq(2,7,3,4,5,6,8,1)
c--- s  u  > d  s
      ms(sq,uq,dq,sq)=ms(bq,uq,dq,bq)
c--- b  c  > s  b
      ms(bq,cq,sq,bq)=ms(bq,uq,dq,bq)
c--- d  c  > s  d
      ms(dq,cq,sq,dq)=ms(bq,uq,dq,bq)
c--- u  c  > d  c
      ms(uq,cq,dq,cq)=WZccmsq(1,7,3,4,5,6,8,2)
c--- u  c  > u  s
      ms(uq,cq,uq,sq)=WZccmsq(2,8,3,4,5,6,7,1)
c--- c  u  > c  d
      ms(cq,uq,cq,dq)=ms(uq,cq,uq,sq)
c--- c  u  > s  u
      ms(cq,uq,sq,uq)=ms(uq,cq,dq,cq)

C***********************2q-2g-processes*********************************
c--- u   d~ > g  g
      ms(uq,da,gg,gg)=aveqq*half*WZggmsq(1,2,3,4,5,6,7,8)
c--- c   s~ > g  g
      ms(cq,sa,gg,gg)=ms(uq,da,gg,gg)
c--- d~  u  > g  g
      ms(da,uq,gg,gg)=aveqq*half*WZggmsq(2,1,3,4,5,6,7,8)
c--- s~  c  > g  g
      ms(sa,cq,gg,gg)=ms(da,uq,gg,gg)

c--- g  g  > u~  d
      ms(gg,gg,ua,dq)=avegg*WZggmsq(7,8,3,4,5,6,2,1)
c--- g  g  > c~  s
      ms(gg,gg,ac,sq)=ms(gg,gg,ua,dq)

c--- u  g  > d g
      ms(uq,gg,dq,gg)=aveqg*WZggmsq(1,7,3,4,5,6,2,8)
c--- c  g  > s g
      ms(cq,gg,sq,gg)=ms(uq,gg,dq,gg)
c--- g  u > d g
      ms(gg,uq,dq,gg)=aveqg*WZggmsq(2,7,3,4,5,6,1,8)
c--- g  c > s g
      ms(gg,cq,sq,gg)=ms(gg,uq,dq,gg)

c--- g  d~  > u~  g
      ms(gg,da,ua,gg)=aveqg*WZggmsq(7,2,3,4,5,6,1,8)
c--- g  s~  > c~  g
      ms(gg,sa,ac,gg)=ms(gg,da,ua,gg)
c--- d~ g  > u~  g
      ms(da,gg,ua,gg)=aveqg*WZggmsq(7,1,3,4,5,6,2,8)
c--- s~ g > c~  g
      ms(sa,gg,ac,gg)=ms(da,gg,ua,gg)

c      call madcheckp(ms)
c      write(6,*) 'qq'
c      include 'qq+'
c      write(6,*) 'aa'
c      include 'aa+'
c      write(6,*) 'qa'
c      include 'qa+'
c      write(6,*) 'aq'
c      include 'aq+'

C***********************************************************************

      elseif (nwz == -1) then

C***********************************************************
C-----Identity table antiquark-quark
c--- u~  u > u  d~ and c  c~ > c  s~
      ms(ua,uq,uq,da)=WZuuidmsq(8,7,3,4,5,6,1,2)
      ms(ac,cq,cq,sa)=ms(ua,uq,uq,da)
c--- u~  d > u  u~ and c~  s > c  c~
      ms(ua,dq,uq,ua)=WZuuidmsq(2,1,3,4,5,6,7,8)
      ms(ac,sq,cq,ac)=ms(ua,dq,uq,ua)
c--- d~  d > u  d~ and s~  s > c  s~
      ms(da,dq,uq,da)=WZddidmsq(2,7,3,4,5,6,1,8)
      ms(sa,sq,cq,sa)=ms(da,dq,uq,da)
c--- u~  d > d  d~ and c~  s > s  s~
      ms(ua,dq,dq,da)=WZddidmsq(2,1,3,4,5,6,7,8)
      ms(ac,sq,sq,sa)=ms(ua,dq,dq,da)

c--- d~  d > c  s~ and s~  s > u  d~
      ms(da,dq,cq,sa)=WZbbmsq(8,7,3,4,5,6,1,2)
      ms(sa,sq,uq,da)=ms(da,dq,cq,sa)
      ms(ba,bq,uq,da)=ms(da,dq,cq,sa)
      ms(ba,bq,cq,sa)=ms(da,dq,cq,sa)
c--- u~  d  > s  s~ and c~  s  > d  d~ and c~  s > b  b~ and u~  d > b  b~
      ms(ua,dq,sq,sa)=WZbbmsq(2,1,3,4,5,6,7,8)
      ms(ac,sq,dq,da)=ms(ua,dq,sq,sa)
      ms(ua,dq,bq,ba)=ms(ua,dq,sq,sa)
      ms(ac,sq,bq,ba)=ms(ua,dq,sq,sa)
c--- s~  d  > u  s~
      ms(sa,dq,uq,sa)=WZbbmsq(2,7,3,4,5,6,1,8)
c--- d~  s  > c  d~
      ms(da,sq,cq,da)=ms(sa,dq,uq,sa)
c--- b~  s > c  b~
      ms(ba,sq,cq,ba)=ms(sa,dq,uq,sa)
c--- b~  d   > u  b~
      ms(ba,dq,uq,ba)=ms(sa,dq,uq,sa)
c--- c~ d  > d  s~
      ms(ac,dq,dq,sa)=WZbbmsq(8,1,3,4,5,6,7,2)
c--- u~  s  > s  d~
      ms(ua,sq,sq,da)=ms(ac,dq,dq,sa)
c--- u~  b  > b  d~
      ms(ua,bq,bq,da)=ms(ua,sq,sq,da)
c--- c~  b  > b  s~
      ms(ac,bq,bq,sa)=ms(ua,sq,sq,da)
c--- c~  c > u  d~ and u~  u > c  s~ and b~  b > u  d~ and b~  b > c  s~
      ms(ua,uq,cq,sa)=WZccmsq(8,7,3,4,5,6,1,2)
      ms(ac,cq,uq,da)=ms(ua,uq,cq,sa)
c--- u~  d > c  c~ and c~  s > u  u~ and
      ms(ua,dq,cq,ac)=WZccmsq(2,1,3,4,5,6,7,8)
      ms(ac,sq,uq,ua)=ms(ua,dq,cq,ac)
      ms(ac,sq,bq,ba)=ms(ua,dq,cq,ac)
c--- c~  u > u  s~
      ms(ac,uq,uq,sa)=WZccmsq(8,1,3,4,5,6,7,2)
c--- u~  c  > c  d~
      ms(ua,cq,cq,da)=ms(ac,uq,uq,sa)
c--- c~  d  > u  c~
      ms(ac,dq,uq,ac)=WZccmsq(2,7,3,4,5,6,1,8)
c--- u~  s  > c  u~
      ms(ua,sq,cq,ua)=ms(ac,dq,uq,ac)

C***********************************************************
C-----Identity table antiquark-antiquark
c--- d~ u~ > d~ d~
      ms(da,ua,da,da)=half*WZddidmsq(7,2,3,4,5,6,1,8)
c--- s~ c~ > s~ s~
      ms(sa,ac,sa,sa)=ms(da,ua,da,da)
c--- u~ d~ > d~ d~
      ms(ua,da,da,da)=half*WZddidmsq(7,1,3,4,5,6,2,8)
c--- c~ s~ > s~ s~
      ms(ac,sa,sa,sa)=ms(ua,da,da,da)
c--- u~ u~ > u~ d~
      ms(ua,ua,ua,da)=WZuuidmsq(8,1,3,4,5,6,2,7)
c--- c~ c~ > c~ s~
      ms(ac,ac,ac,sa)=ms(ua,ua,ua,da)
c--- u~ s~ > d~ s~
      ms(ua,sa,da,sa)=WZbbmsq(7,1,3,4,5,6,2,8)
c--- c~ b~ > s~ b~
      ms(ac,ba,sa,ba)=ms(ua,sa,da,sa)
c--- u~ b~ > d~ b~
      ms(ua,ba,da,ba)=ms(ua,sa,da,sa)
c--- c~ d~ > s~ u~
      ms(ac,da,sa,da)=ms(ua,sa,da,sa)
c--- d~ c~ > s~  d~
      ms(da,ac,sa,da)=WZbbmsq(7,2,3,4,5,6,1,8)
c--- b~ u~ > d~ b~
      ms(ba,ua,da,ba)=ms(da,ac,sa,da)
c--- b~ c~ > s~ b~
      ms(ba,ac,sa,ba)=ms(da,ac,sa,da)
c--- s~ u~ > d~ s~
      ms(sa,ua,da,sa)=ms(da,ac,sa,da)
c--- c~ u~ > d~ c~ and c~ u~ > s~ u~
c--- u~ c~ > d~ c~ and u~ c~ > s~ u~
      ms(ac,ua,da,ac)=WZccmsq(7,2,3,4,5,6,1,8)
      ms(ac,ua,sa,ua)=WZccmsq(7,1,3,4,5,6,2,8)
      ms(ua,ac,da,ac)=ms(ac,ua,sa,ua)
      ms(ua,ac,sa,ua)=ms(ac,ua,da,ac)

C***********************************************************
C-----Identity tablec quark-quark
c--- d  u  > u  u
      ms(dq,uq,uq,uq)=half*WZuuidmsq(1,7,3,4,5,6,8,2)
c--- s  c  > c  c
      ms(sq,cq,cq,cq)=ms(dq,uq,uq,uq)
c--- u  d  > u  u
      ms(uq,dq,uq,uq)=half*WZuuidmsq(2,7,3,4,5,6,8,1)
c--- c  s  > c  c
      ms(cq,sq,cq,cq)=ms(uq,dq,uq,uq)
c--- d  d  > u  d
      ms(dq,dq,uq,dq)=WZddidmsq(1,7,3,4,5,6,8,2)
c--- s  s  > c  s
      ms(sq,sq,cq,sq)=ms(dq,dq,uq,dq)
c--- d  b  > u  b
      ms(dq,bq,uq,bq)=WZbbmsq(1,7,3,4,5,6,8,2)
c--- d  c  > u  c
      ms(dq,cq,uq,cq)=WZccmsq(1,7,3,4,5,6,8,2)
c--- s  b  > c  b
      ms(sq,bq,cq,bq)=ms(dq,bq,uq,bq)
c--- s  u  > c  u
      ms(sq,uq,cq,uq)=ms(dq,cq,uq,cq)
c--- b  d  > u  b
      ms(bq,dq,uq,bq)=WZbbmsq(2,7,3,4,5,6,8,1)
c--- b  s  > c  b
      ms(bq,sq,cq,bq)=ms(bq,dq,uq,bq)
c--- c  d  > u  c
      ms(cq,dq,uq,cq)=WZccmsq(2,7,3,4,5,6,8,1)
c--- u  s  > c  u
      ms(uq,sq,cq,uq)=ms(cq,dq,uq,cq)
c--- d  s  > u  s
      ms(dq,sq,uq,sq)=WZbbmsq(1,7,3,4,5,6,8,2)
c--- d  s  > d  s
      ms(dq,sq,dq,cq)=WZbbmsq(2,8,3,4,5,6,7,1)
c--- s  d  > s  u
      ms(sq,dq,sq,uq)=ms(dq,sq,dq,cq)
c--- s  d  > c  u
      ms(sq,dq,cq,dq)=ms(dq,sq,uq,sq)

C***********************************************************
C-----Identity tablec quark-antiquark
c--- u  u~ > u  d~ and c  c~ > c  s~
      ms(uq,ua,uq,da)=WZuuidmsq(8,7,3,4,5,6,2,1)
      ms(cq,ac,cq,sa)=ms(uq,ua,uq,da)
c--- d  d~ > u  d~ and s  s~ > c  s~
      ms(dq,da,uq,da)=WZddidmsq(1,7,3,4,5,6,2,8)
      ms(sq,sa,cq,sa)=ms(dq,da,uq,da)
c--- d  d~ > c  s~ and s  s~ > u  d~
      ms(dq,da,cq,sa)=WZbbmsq(8,7,3,4,5,6,2,1)
      ms(sq,sa,uq,da)=ms(dq,da,cq,sa)
      ms(bq,ba,uq,da)=ms(dq,da,cq,sa)
      ms(bq,ba,cq,sa)=ms(dq,da,cq,sa)
c--- d  u~ > u  u~ and s  c~ > c  c~
      ms(dq,ua,uq,ua)=WZuuidmsq(1,2,3,4,5,6,7,8)
      ms(sq,ac,cq,ac)=ms(dq,ua,uq,ua)
c--- d  u~ > c  c~ and s  c~ > u  u~ and s  c~ > b  b~ and d  u~ > b  b~
      ms(dq,ua,cq,ac)=WZccmsq(1,2,3,4,5,6,7,8)
      ms(sq,ac,uq,ua)=ms(dq,ua,cq,ac)
c--- d  u~ > d  d~ and s  c~ > s  s~
      ms(dq,ua,dq,da)=WZddidmsq(1,2,3,4,5,6,7,8)
      ms(sq,ac,sq,sa)=ms(dq,ua,dq,da)
c--- d  u~ > s  s~
      ms(dq,ua,sq,sa)=WZbbmsq(1,2,3,4,5,6,7,8)
c--- s  c~ > d  d~
      ms(sq,ac,dq,da)=ms(dq,ua,sq,sa)
      ms(sq,ac,bq,ba)=ms(dq,ua,sq,sa)
      ms(dq,ua,bq,ba)=ms(dq,ua,sq,sa)
c--- d  s~ > u  s~
      ms(dq,sa,uq,sa)=WZbbmsq(1,7,3,4,5,6,2,8)
c--- d  b~ > u  b~
      ms(dq,ba,uq,ba)=ms(dq,sa,uq,sa)
c--- s  d~ > c  d~
      ms(sq,da,cq,da)=ms(dq,sa,uq,sa)
c--- s  b~ > c  b~
      ms(sq,ba,cq,ba)=ms(dq,sa,uq,sa)
c--- d  c~ > d  s~
      ms(dq,ac,dq,sa)=WZbbmsq(8,2,3,4,5,6,7,1)
c--- s  u~ > s  d~
      ms(sq,ua,sq,da)=ms(dq,ac,dq,sa)
c--- b  u~ > b  d~
      ms(bq,ua,bq,da)=ms(dq,ac,dq,sa)
c--- b  c~ > b  s~
      ms(bq,ac,bq,sa)=ms(dq,ac,dq,sa)
c--- d  c~ > u  c~
      ms(dq,ac,uq,ac)=WZccmsq(1,7,3,4,5,6,2,8)
c--- s  u~ > c  u~
      ms(sq,ua,cq,ua)=ms(dq,ac,uq,ac)
c--- c  c~ > u  d~ and u  u~ > c  s~ and b  b~ > u  d~ and b  b~ > c  s~
      ms(uq,ua,cq,sa)=WZccmsq(8,7,3,4,5,6,2,1)
      ms(cq,ac,uq,da)=ms(uq,ua,cq,sa)
c--- u  c~ > u  s~
      ms(uq,ac,uq,sa)=WZccmsq(8,2,3,4,5,6,7,1)
c--- c  u~ > c  d~
      ms(cq,ua,cq,da)=ms(uq,ac,uq,sa)


C***********************2q-2g-processes*********************************
c--- d   u~ > g  g
      ms(dq,ua,gg,gg)=aveqq*half*WZggmsq(1,2,3,4,5,6,7,8)
c--- s   c~ > g  g
      ms(sq,ac,gg,gg)=ms(dq,ua,gg,gg)
c--- u~  d  > g  g
      ms(ua,dq,gg,gg)=aveqq*half*WZggmsq(2,1,3,4,5,6,7,8)
c--- c~  s  > g  g
      ms(ac,sq,gg,gg)=ms(ua,dq,gg,gg)

c--- g  g  > d~  u
      ms(gg,gg,da,uq)=avegg*WZggmsq(7,8,3,4,5,6,2,1)
c--- g  g  > s~  c
      ms(gg,gg,sa,cq)=ms(gg,gg,da,uq)

c--- d  g  > u g
      ms(dq,gg,uq,gg)=aveqg*WZggmsq(1,7,3,4,5,6,2,8)
c--- s  g  > c g
      ms(sq,gg,cq,gg)=ms(dq,gg,uq,gg)
c--- g  d > u g
      ms(gg,dq,uq,gg)=aveqg*WZggmsq(2,7,3,4,5,6,1,8)
c--- g  s > c g
      ms(gg,sq,cq,gg)=ms(gg,dq,uq,gg)

c--- g  u~  > d~  g
      ms(gg,ua,da,gg)=aveqg*WZggmsq(7,2,3,4,5,6,1,8)
c--- g  c~  > s~  g
      ms(gg,ac,sa,gg)=ms(gg,ua,da,gg)
c--- u~ g  > d~  g
      ms(ua,gg,da,gg)=aveqg*WZggmsq(7,1,3,4,5,6,2,8)
c--- c~ g > s~  g
      ms(ac,gg,sa,gg)=ms(ua,gg,da,gg)

C***********************************************************************

c      write(6,*) 'qq'
c      include 'qq-'
c      write(6,*) 'aa'
c      include 'aa-'
c      write(6,*) 'qa'
c      include 'qa-'
c      write(6,*) 'aq'
c      include 'aq-'
c      call madcheckm(ms)

      endif

      msq(:,:)=0._dp
      do i=-nf,nf
      do j=-nf,nf
      do k=-nf,nf
      do l=-nf,nf
      msq(i,j)=msq(i,j)+ms(i,j,k,l)
      enddo
      enddo
      enddo
      enddo

      return
      end






