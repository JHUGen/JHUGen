      double complex function ehsva4(s,t,u)
C     ehsv:EqnA.8
      implicit none
      double precision s,t,u
      double complex ehsvb4
      ehsva4=ehsvb4(s,t,u)+ehsvb4(u,s,t)+ehsvb4(t,u,s)
      return 
      end

      double complex function ehsva2(s,t,u)
C     ehsv:EqnA.9
      implicit none
      double precision s,t,u
      double complex ehsvb2
      ehsva2=ehsvb2(s,t,u)+ehsvb2(s,u,t)
      return 
      end

      double complex function ehsvb4(s,t,u)
      implicit none
C     ehsv:EqnA.10
      include 'masses.f'
      double precision hmass2,s,t,u
      double complex w2,w3
      hmass2=s+t+u
c--- The Fermilab preprint has w2(s), but it makes no difference due
c--- to symmetrization in ehsva4 above      
      ehsvb4=mbsq/hmass2*(-2d0/3d0
     . +(mbsq/hmass2-0.25d0)*(w2(t)-w2(hmass2)+w3(s,t,u,hmass2)))
      return 
      end

      double complex function ehsvb2(s,t,u)
C     ehsv:EqnA.11
      implicit none
      include 'masses.f'
      double precision hmass2,s,t,u
      double complex w1,w2,w3
      hmass2=s+t+u
      ehsvb2=mbsq/hmass2**2*(s*(u-s)/(s+u)
     . +2d0*u*t*(u+2d0*s)/(s+u)**2*(w1(t)-w1(hmass2))
     . +(mbsq-0.25d0*s)
     . *(0.5d0*w2(s)+0.5d0*w2(hmass2)-w2(t)+w3(s,t,u,hmass2))
     . +s**2*(2d0*mbsq/(s+u)**2-0.5d0/(s+u))*(w2(t)-w2(hmass2))
     . +0.5d0*u*t/s*(w2(hmass2)-2d0*w2(t))
     . +0.125d0*(s-12d0*mbsq-4d0*u*t/s)*w3(t,s,u,hmass2))
      return 
      end

      double complex function ehsva5(s,t,u)
C     ehsv:EqnA.14
      implicit none
      include 'masses.f'
      double precision hmass2,s,t,u
      double complex w1,w2
      hmass2=s+t+u
      ehsva5=mbsq/hmass2*(4d0+4d0*s/(u+t)*(w1(s)-w1(hmass2))
     . +(1d0-4d0*mbsq/(u+t))*(w2(s)-w2(hmass2)))
      return 
      end


      double complex function w1(s)
C     ehsv:EqnA.19
      implicit none
      include 'constants.f'
      include 'masses.f'
      double precision rat,s,temp,acosh,asinh
      rat=4d0*mbsq/s
      temp=dsqrt(dabs(1d0/rat))
      if (rat .lt. 0d0) then
          w1=2d0*dsqrt(1d0-rat)*asinh(temp)
      elseif (rat .gt. 1d0) then
          w1=2d0*dsqrt(rat-1d0)*asin(temp)
      else 
          temp=2d0*acosh(temp)
          w1=dsqrt(1d0-rat)*dcmplx(temp,-pi)
       endif
      return
      end

      double complex function w2(s)
C     ehsv:EqnA.20
      implicit none
      include 'constants.f'
      include 'masses.f'
      double precision rat,s,tempr,tempi,acosh,asinh
      rat=s/(4d0*mbsq)
      tempr=dsqrt(dabs(rat))
      if (rat .lt. 0d0) then
          tempr=asinh(tempr)
          w2=4d0*tempr**2
      elseif (rat .gt. 1d0) then
          tempr=acosh(tempr)
          tempi=-4d0*tempr*pi
          tempr=+4d0*tempr**2-pi**2
          w2=dcmplx(tempr,tempi)
      else 
          tempr=asin(tempr)
          w2=-4d0*tempr**2
      endif
      return
      end

      double complex function w3(s,t,u,varg)
C     ehsv:EqnA.17
      implicit none
      double complex i3
      double precision s,t,u,varg
      w3=i3(s,t,u,varg)-i3(s,t,u,s)-i3(s,t,u,u)
      return
      end


      double complex function i3(s,t,u,varg)
C     ehsv:EqnA.21
      implicit none
      include 'constants.f'
      include 'masses.f'
      double precision s,t,u,varg,rat,al,be,ga,r,theta,phi
      double precision arg1,arg2,arg3,arg4,ddilog,arg
      double complex cli2,zth,zph
      rat=4d0*mbsq/varg
      if (rat .lt. 0d0) then
           be=0.5d0*(1d0+dsqrt(1d0+4d0*t*mbsq/(u*s)))
           ga=0.5d0*(1d0+dsqrt(1d0-rat))
           arg1=ga/(ga+be-1d0)
           arg2=(ga-1d0)/(ga+be-1d0)
           arg3=(be-ga)/be
           arg4=(be-ga)/(be-1d0)
           i3=2d0/(2d0*be-1d0)
     .     *(-ddilog(arg1)+ddilog(arg2)+ddilog(arg3)-ddilog(arg4)
     .     +0.5d0*(dlog(be)**2-dlog(be-1d0)**2)
     .     +dlog(ga)*dlog((ga+be-1d0)/be)
     .     +dlog(ga-1d0)*dlog((be-1d0)/(ga+be-1d0)))
      elseif (rat .gt. 1d0) then
           be=0.5d0*(1d0+dsqrt(1d0+4d0*t*mbsq/(u*s)))
           al=dsqrt(rat-1d0)           
           r=dsqrt((al**2+1d0)/(al**2+(2d0*be-1d0)**2))
           arg=r*(al**2+2d0*be-1d0)/(1d0+al**2)
           if (arg .ge. 1d0) then
             phi=0d0
           else
             phi=dacos(arg)
           endif
           arg=r*(al**2-2d0*be+1d0)/(1d0+al**2)
           if (arg .ge. 1d0) then
             theta=0d0
           else
             theta=dacos(arg)
           endif
           zth=r*dcmplx(cos(theta),sin(theta))
           zph=r*dcmplx(cos(phi),sin(phi))
           i3=2d0/(2d0*be-1d0)
     .     *(2d0*dble(cli2(zth))-2d0*dble(cli2(zph))
     .     +(phi-theta)*(phi+theta-pi))
      else
           be=0.5d0*(1d0+dsqrt(1d0+4d0*t*mbsq/(u*s)))
           ga=0.5d0*(1d0+dsqrt(1d0-rat))
           arg1=ga/(ga+be-1d0)
           arg2=(ga-1d0)/(ga+be-1d0)
           arg3=ga/(ga-be)
           arg4=(ga-1d0)/(ga-be)
      
           i3=2d0/(2d0*be-1d0)
     .     *(-ddilog(arg1)+ddilog(arg2)+ddilog(arg3)-ddilog(arg4)
     .     +dlog(ga/(1d0-ga))*dlog((ga+be-1d0)/(be-ga))
     .     -im*pi*dlog((ga+be-1d0)/(be-ga)))
      endif

      return
      end

      double precision function acosh(y)
      implicit none
      double precision y      
      acosh=dlog(y+dsqrt(y**2-1d0))
      return
      end

      double precision function asinh(y)
      implicit none
      double precision y      
      asinh=dlog(y+dsqrt(y**2+1d0))
      return
      end

