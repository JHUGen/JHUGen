      function ehsva4(s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: ehsva4
C     ehsv:EqnA.8

      real(dp):: s,t,u
      complex(dp):: ehsvb4
      ehsva4=ehsvb4(s,t,u)+ehsvb4(u,s,t)+ehsvb4(t,u,s)
      return
      end

      function ehsva2(s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: ehsva2
C     ehsv:EqnA.9

      real(dp):: s,t,u
      complex(dp):: ehsvb2
      ehsva2=ehsvb2(s,t,u)+ehsvb2(s,u,t)
      return
      end

      function ehsvb4(s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: ehsvb4

C     ehsv:EqnA.10
      include 'masses.f'
      real(dp):: hmass2,s,t,u
      complex(dp):: w2,w3
      hmass2=s+t+u
c--- The Fermilab preprint has w2(s), but it makes no difference due
c--- to symmetrization in ehsva4 above
      ehsvb4=mbsq/hmass2*(-2._dp/3._dp
     & +(mbsq/hmass2-0.25_dp)*(w2(t)-w2(hmass2)+w3(s,t,u,hmass2)))
      return
      end

      function ehsvb2(s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: ehsvb2
C     ehsv:EqnA.11

      include 'masses.f'
      real(dp):: hmass2,s,t,u
      complex(dp):: w1,w2,w3
      hmass2=s+t+u
      ehsvb2=mbsq/hmass2**2*(s*(u-s)/(s+u)
     & +2._dp*u*t*(u+2._dp*s)/(s+u)**2*(w1(t)-w1(hmass2))
     & +(mbsq-0.25_dp*s)
     & *(0.5_dp*w2(s)+0.5_dp*w2(hmass2)-w2(t)+w3(s,t,u,hmass2))
     & +s**2*(2._dp*mbsq/(s+u)**2-0.5_dp/(s+u))*(w2(t)-w2(hmass2))
     & +0.5_dp*u*t/s*(w2(hmass2)-2._dp*w2(t))
     & +0.125_dp*(s-12._dp*mbsq-4._dp*u*t/s)*w3(t,s,u,hmass2))
      return
      end

      function ehsva5(s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: ehsva5
C     ehsv:EqnA.14

      include 'masses.f'
      real(dp):: hmass2,s,t,u
      complex(dp):: w1,w2
      hmass2=s+t+u
      ehsva5=mbsq/hmass2*(4._dp+4._dp*s/(u+t)*(w1(s)-w1(hmass2))
     & +(1._dp-4._dp*mbsq/(u+t))*(w2(s)-w2(hmass2)))
      return
      end


      function w1(s)
      implicit none
      include 'types.f'
      complex(dp):: w1
C     ehsv:EqnA.19

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      real(dp):: rat,s,temp,acosh,asinh
      rat=4._dp*mbsq/s
      temp=sqrt(abs(1._dp/rat))
      if (rat < 0._dp) then
          w1=2._dp*sqrt(1._dp-rat)*asinh(temp)
      elseif (rat > 1._dp) then
          w1=2._dp*sqrt(rat-1._dp)*asin(temp)
      else
          temp=2._dp*acosh(temp)
          w1=sqrt(1._dp-rat)*cplx2(temp,-pi)
       endif
      return
      end

      function w2(s)
      implicit none
      include 'types.f'
      complex(dp):: w2
C     ehsv:EqnA.20

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      real(dp):: rat,s,tempr,tempi,acosh,asinh
      rat=s/(4._dp*mbsq)
      tempr=sqrt(abs(rat))
      if (rat < 0._dp) then
          tempr=asinh(tempr)
          w2=4._dp*tempr**2
      elseif (rat > 1._dp) then
          tempr=acosh(tempr)
          tempi=-4._dp*tempr*pi
          tempr=+4._dp*tempr**2-pi**2
          w2=cplx2(tempr,tempi)
      else
          tempr=asin(tempr)
          w2=-4._dp*tempr**2
      endif
      return
      end

      function w3(s,t,u,varg)
      implicit none
      include 'types.f'
      complex(dp):: w3
C     ehsv:EqnA.17

      complex(dp):: i3
      real(dp):: s,t,u,varg
      w3=i3(s,t,u,varg)-i3(s,t,u,s)-i3(s,t,u,u)
      return
      end


      function i3(s,t,u,varg)
      implicit none
      include 'types.f'
      complex(dp):: i3
C     ehsv:EqnA.21

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      real(dp):: s,t,u,varg,rat,al,be,ga,r,theta,phi
      real(dp):: arg1,arg2,arg3,arg4,ddilog,arg
      complex(dp):: cli2,zth,zph
      rat=4._dp*mbsq/varg
      if (rat < 0._dp) then
           be=0.5_dp*(1._dp+sqrt(1._dp+4._dp*t*mbsq/(u*s)))
           ga=0.5_dp*(1._dp+sqrt(1._dp-rat))
           arg1=ga/(ga+be-1._dp)
           arg2=(ga-1._dp)/(ga+be-1._dp)
           arg3=(be-ga)/be
           arg4=(be-ga)/(be-1._dp)
           i3=2._dp/(2._dp*be-1._dp)
     &     *(-ddilog(arg1)+ddilog(arg2)+ddilog(arg3)-ddilog(arg4)
     &     +0.5_dp*(log(be)**2-log(be-1._dp)**2)
     &     +log(ga)*log((ga+be-1._dp)/be)
     &     +log(ga-1._dp)*log((be-1._dp)/(ga+be-1._dp)))
      elseif (rat > 1._dp) then
           be=0.5_dp*(1._dp+sqrt(1._dp+4._dp*t*mbsq/(u*s)))
           al=sqrt(rat-1._dp)
           r=sqrt((al**2+1._dp)/(al**2+(2._dp*be-1._dp)**2))
           arg=r*(al**2+2._dp*be-1._dp)/(1._dp+al**2)
           if (arg >= 1._dp) then
             phi=0._dp
           else
             phi=acos(arg)
           endif
           arg=r*(al**2-2._dp*be+1._dp)/(1._dp+al**2)
           if (arg >= 1._dp) then
             theta=0._dp
           else
             theta=acos(arg)
           endif
           zth=r*cplx2(cos(theta),sin(theta))
           zph=r*cplx2(cos(phi),sin(phi))
           i3=2._dp/(2._dp*be-1._dp)
     &     *(2._dp*real(cli2(zth))-2._dp*real(cli2(zph))
     &     +(phi-theta)*(phi+theta-pi))
      else
           be=0.5_dp*(1._dp+sqrt(1._dp+4._dp*t*mbsq/(u*s)))
           ga=0.5_dp*(1._dp+sqrt(1._dp-rat))
           arg1=ga/(ga+be-1._dp)
           arg2=(ga-1._dp)/(ga+be-1._dp)
           arg3=ga/(ga-be)
           arg4=(ga-1._dp)/(ga-be)

           i3=2._dp/(2._dp*be-1._dp)
     &     *(-ddilog(arg1)+ddilog(arg2)+ddilog(arg3)-ddilog(arg4)
     &     +log(ga/(1._dp-ga))*log((ga+be-1._dp)/(be-ga))
     &     -im*pi*log((ga+be-1._dp)/(be-ga)))
      endif

      return
      end

      function acosh(y)
      implicit none
      include 'types.f'
      real(dp):: acosh

      real(dp):: y
      acosh=log(y+sqrt(y**2-1._dp))
      return
      end

      function asinh(y)
      implicit none
      include 'types.f'
      real(dp):: asinh

      real(dp):: y
      asinh=log(y+sqrt(y**2+1._dp))
      return
      end

