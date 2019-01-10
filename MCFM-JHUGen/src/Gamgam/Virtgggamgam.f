      double precision function virtgamgam(s,t,u)
      implicit none
      include 'epinv.f'
      include 'constants.f'
      include 'scale.f'
      include 'qcdcouple.f'
      include 'nflav.f'
      integer h1,h2,h3,h4
      double precision temp,s,t,u
      double complex m1(2,2,2,2),m2fin(2,2,2,2),lnrat,
     & I1ggtogamgam,xlog,FL(2,2,2,2),FSL(2,2,2,2)

      xlog=lnrat(musq,-s)
c--- This is taken from hep-ph/0109078 Eq.(2.11); note however that the log
c--- proportional to the beta-function coefficient is added in Eq. (2.11)
c--- and subtracted again in Eq. (4.5), therefore we omit it here.
      I1ggtogamgam=
     & -xn*((epinv**2+epinv*xlog+0.5d0*xlog**2)
     & +(11d0/6d0-dfloat(nflav)/(3d0*xn))*epinv)

      temp=0d0
      call M1fill(s,t,u,m1)

c--- testing the 1-loop matrix element contribution up to Order(ep^2)      
c      call oneloopep(s,t,u,IxM1ep)

c--- pure 2-loop contribution
      call twoloop(s,t,u,FL,FSL)

c--- Compare m1 with m1ep
c      do h1=1,2
c      do h2=1,2
c      do h3=1,2
c      do h4=1,2
c      write(6,*) h1,h2,h3,h4,m1(h1,h2,h3,h4)*I1ggtogamgam,
c     &           IxM1ep(h1,h2,h3,h4)
c      enddo
c      enddo
c      enddo
c      enddo
c      pause
 
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      m2fin(h1,h2,h3,h4)=xn*FL(h1,h2,h3,h4)-FSL(h1,h2,h3,h4)/xn
      temp=temp+2d0*ason2pi*dble(dconjg(m1(h1,h2,h3,h4))
     & *(m1(h1,h2,h3,h4)*I1ggtogamgam+m2fin(h1,h2,h3,h4)))
c      temp=temp+2d0*ason2pi*dble(dconjg(m1(h1,h2,h3,h4))
c     & *(IxM1ep(h1,h2,h3,h4)+m2fin(h1,h2,h3,h4)))
      enddo
      enddo
      enddo
      enddo
      
      virtgamgam=temp
      
      return
      end
