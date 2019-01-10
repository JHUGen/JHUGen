      double complex function LLs1(xsq,s,t,qsq,ms2,mc2)
c===============================c
c refers to bx1:                c
c   ms2                         c
c      z     s     /            c
c       z         /             c
c        z-------/              c
c        z       |              c
c   ms2  z       |     t        c
c        z       |              c
c        y=======\\             c
c       y         \\            c
c      y           \\           c
c  qsq                          c
c                               c
c                               c
c===============================c
      implicit none
      double precision xsq,s,t,qsq,ms2,mc2,cf0,cf1,cf2,cf3,cf4
      double complex bx1,tr1,tr3,box1,tri1,tri2,tri3,tri4

      tri1=tr1(xsq,t,mc2)
      tri2=tr3(s,qsq,mc2,ms2)
      tri3=tr3(t,qsq,ms2,mc2)
      tri4=tr1(xsq,s,ms2)
      box1=bx1(xsq,s,t,qsq,ms2,mc2)

      cf0=(-4*(mc2**2*ms2 + t*(ms2*(qsq - s) + s*(-qsq + s + t)) + 
     &      mc2*(ms2**2 + s*(qsq - t) - ms2*(qsq + s + t))))/
     &  ((ms2 - s)**2*(mc2 - t)**2)

      cf1=2/(ms2 - s)
      cf2=(2*(-mc2**2 + (-qsq + s)*t + mc2*(-2*ms2 + qsq + s + t)))/
     &  ((ms2 - s)*(mc2 - t)**2)
      cf3=(2*(-2*mc2*ms2 - ms2**2 + s*(-qsq + t) + ms2*(qsq + s + t)))/
     &  ((ms2 - s)**2*(mc2 - t))
      cf4=2/(mc2 - t)

      LLs1=(+2d0*box1
     &        +cf1*tri1
     &        +cf2*tri2
     &        +cf3*tri3
     &        +cf4*tri4)/2d0/cf0

      return
      end


      double complex function LLs2(xsq,s,t,qsq,ms2,mc2)
c===============================c
c   qsq=-s-t-msq                c
c                      ms2      c
c      $            #           c
c       $          #            c
c        $#########             c
c        ||       |             c
c        ||       |     t       c
c        ||       |             c
c        /========\\            c
c       /          \\           c
c      /            \\          c
c                      mc2      c
c            s                  c
c                               c
c===============================c
      implicit none
      double precision xsq,s,t,qsq,ms2,mc2,cf0,cf1,cf2,cf3,cf4
      double complex bx2,box1,tr2,tr3,tr4,tr5,tri1,tri2,tri3,tri4

      tri1=tr2(xsq,t,ms2,mc2)
      tri2=tr3(s,qsq,ms2,mc2)
      tri3=tr5(t,qsq,mc2,ms2)
      tri4=tr4(s,mc2)
      box1=bx2(xsq,s,t,qsq,ms2,mc2)

c      write(6,*) 'tri4=',tri4

      cf0=(-4*mc2*qsq*(mc2 - ms2 + qsq) + 
     &    4*(mc2 - ms2)*qsq*s + 
     &    4*((ms2 + qsq - s)*s + mc2*(-ms2 + qsq + s))*t - 
     &    4*s*t**2)/
     &  ((mc2 - s)**2*(mc2**2 + (ms2 - t)**2 - 
     &      2*mc2*(ms2 + t)))

      cf1=2/(mc2 - s)

      cf2=(2*(-ms2**2 + mc2*(ms2 + qsq - s) + (-qsq + s)*t + 
     &      ms2*(qsq + s + t)))/
     &  ((mc2 - s)*(mc2**2 + (ms2 - t)**2 - 
     &      2*mc2*(ms2 + t)))

      cf3=(-2*(qsq - t)*(mc2**2 + s*(ms2 - t) - 
     &      mc2*(ms2 - 2*qsq + s + t)))/
     &  ((mc2 - s)**2*(mc2**2 + (ms2 - t)**2 - 
     &      2*mc2*(ms2 + t)))

      cf4=(2*(-mc2**2 + s*(-ms2 + t) + 
     &      mc2*(ms2 - 2*qsq + s + t)))/
     &  ((mc2 - s)*(mc2**2 + (ms2 - t)**2 - 
     &      2*mc2*(ms2 + t)))

      LLs2=(+2d0*box1
     &        +cf1*tri1
     &        +cf2*tri2
     &        +cf3*tri3
     &        +cf4*tri4
     &        )/2d0/cf0

c      write(6,*) 'tri2=',tri2
c      write(6,*) 'tri3=',tri3
c      write(6,*) 'tri4=',tri4

      return
      end


      double complex function bx1(xsq,s,t,qsq,ms2,mc2)
c===============================c
c                               c
c   ms2                         c
c      z     s     /            c
c       z         /             c
c        z-------/              c
c        z       |              c
c   ms2  z       |     t        c
c        z       |              c
c        y=======\\             c
c       y         \\            c
c      y           \\           c
c  qsq                          c
c                               c
c                               c
c===============================c
      implicit none
      include 'constants.f'
      include 'epinv.f'
      double precision xsq,mu,ms,mc,s,t,qsq,ms2,mc2,A,B,lam
      double complex lnrat,logq

      mu=dsqrt(xsq)
      ms=dsqrt(ms2)
      mc=dsqrt(mc2)

      lam=(mc2+ms2-qsq)/(dsqrt(mc2*ms2))
      A=+1d0/(lam**2-4d0)
      B=+2d0/(lam*(dsqrt(lam**2-4d0)-lam)+4d0)

      if (qsq.gt.((ms+mc)**2)) then
      logq=
     &   -lnrat(-B,1d0)**2
     &   -lnrat(-A,-B)*lnrat(A+B,B)
      else
      logq=
     &   -lnrat(B,1d0)**2
     &   -lnrat(A,B)*lnrat(A+B,B)
      endif

      bx1= - (
     &         -dcmplx(epinv**2)
     &         +epinv*lnrat(mc2-t,mc2)
     &         +epinv*lnrat(ms2-s,ms2)
     &         -epinv*lnrat(xsq,mc*ms)
     &         +lnrat(xsq,mc2)*lnrat(mc2-t,mc2)
     &         +lnrat(xsq,mc2)*lnrat(ms2-s,ms2)
     &         -2d0*lnrat(mu,mc)**2
     &         -2d0*lnrat(mu,mc)*lnrat(mc,ms)
     &   +logq
     &   +0.25d0*lnrat(A,1d0)**2
     &   +0.25d0*lnrat(A+B,1d0)**2
     &   -2d0*lnrat(ms2-s,dsqrt(mc2*ms2))*lnrat(mc2-t,mc2)
     &   +dcmplx(pisq/2d0)
     &   )/(ms2-s)/(mc2-t)

      return
      end


      double complex function bx2(xsq,s,t,qsq,ms2,mc2)
c===============================c
c   qsq=-s-t-msq                c
c                      ms2      c
c      $            #           c
c       $          #            c
c        $#########             c
c        ||       |             c
c        ||       |     t       c
c        ||       |             c
c        /========\\            c
c       /          \\           c
c      /            \\          c
c                      mc2      c
c            s                  c
c                               c
c===============================c
      implicit none
      include 'constants.f'
      include 'epinv.f'

      double precision xsq,s,t,qsq,ms2,mc2,ddilog
      double precision ss,tt,ms,mc,mu,xs,x3
      double precision xsn,xsd,x3n,x3d,xx3
      double complex lnrat
      
      ss=t
      tt=s
      ms=dsqrt(ms2)
      mc=dsqrt(mc2)
      mu=dsqrt(xsq)

      xsn=(1d0-dsqrt(1d0-4d0*ms*mc/(ss-(ms-mc)**2)))
      xsd=(1d0+dsqrt(1d0-4d0*ms*mc/(ss-(ms-mc)**2)))

      xs=-xsn/xsd

      x3n=(1d0-dsqrt(1d0-4d0*ms*mc/(qsq-(ms-mc)**2)))
      x3d=(1d0+dsqrt(1d0-4d0*ms*mc/(qsq-(ms-mc)**2)))

      x3=-x3n/x3d

c      write(6,*) 'xs=',xs
c      write(6,*) 'x3=',x3
c      write(6,*) 'xs*x3=',xs*x3
c      write(6,*) 'xs/x3=',xs/x3

      if (x3.lt.0d0) then
      xx3=-x3
      else
      xx3=x3
      endif

      bx2= + (
     &   -epinv*lnrat(xsn,-xsd)
     &   -2d0*lnrat(xsn,-xsd)*lnrat(mc*mu,mc2-tt)
     &   +2d0*lnrat(xsn,-xsd)*log(1d0-xs**2)
     &   +dcmplx(pisq/2d0)
     &   +dcmplx(ddilog(xs**2))
     &   +lnrat(x3n,-x3d)**2
     &   - 2d0*(+dcmplx(ddilog(xs*x3))
     &          +(lnrat(xsn,-xsd)+lnrat(x3n,-x3d))
     &                          *(log(1d0-xs*x3)))
     &   - 2d0*(+dcmplx(ddilog(xs/x3))
     &          +(lnrat(xsn,-xsd)+lnrat(-x3d,x3n))
     &                          *(lnrat(x3-xs,xx3)))
     &   )*xs/(ms*mc)/(tt-mc2)/(1d0-xs**2)

c      write(6,*) 'x3=',x3
c      write(6,*) 'xs=',xs,x3-xs

      return
      end


      double complex function tr1(xsq,t,mc2)
c===============================c
c                               c
c                               c
c            |                  c
c            |                  c
c            |                  c
c           / \                 c
c          /   \                c
c         /     \               c
c        x=======\\             c
c       x         \\            c
c      x           \\           c
c                               c
c    t                          c
c                               c
c===============================c
      implicit none
      include 'constants.f'
      include 'epinv.f'
      double precision xsq,t,mc2,ddilog,mc,mu
      double complex lnrat,dd1

      mc=dsqrt(mc2)
      mu=dsqrt(xsq)

      if (t.gt.mc2) then
      dd1=-impi*dlog(1d0+mc2/(t-mc2))
      else
      dd1=czip
      endif
      dd1=dd1+dcmplx(ddilog(1d0-mc2/(mc2-t)))

      tr1= (
     &       -dcmplx(epinv**2)
     &       -2d0*epinv*lnrat(mc*mu,mc2-t)
     &       +lnrat(mc2,mc2-t)**2
     &       -2d0*lnrat(mu*mc,mc2-t)**2
     &       +2d0*dd1
     &       -dcmplx(pisqo6)
     &     )/2d0/(mc2-t)

      return
      end


      double complex function tr2(xsq,s,ms2,mc2)
c===============================c
c                               c
c         s                     c
c            $                  c
c            $                  c
c            $                  c
c           x \\                c
c          x   \\               c
c         x     \\              c
c        x-------\\             c
c       x         \\            c
c      x           \\           c
c                               c
c    ms2           mc2          c
c                               c
c===============================c
      implicit none
      include 'constants.f'
      include 'epinv.f'

      double precision xsq,s,ms2,mc2,ddilog
      double precision ms,mc,xs
      double precision xsn,xsd
      double complex lnrat
      
      ms=dsqrt(ms2)
      mc=dsqrt(mc2)

      xsn=(1d0-dsqrt(1d0-4d0*ms*mc/(s-(ms-mc)**2)))
      xsd=(1d0+dsqrt(1d0-4d0*ms*mc/(s-(ms-mc)**2)))

      xs=-xsn/xsd

      tr2= + (
     &   -epinv*lnrat(xsn,-xsd)
     &   -lnrat(xsn,-xsd)*lnrat(xsq,mc*ms)
     &   -0.5d0*lnrat(xsn,-xsd)**2
     &   +2d0*lnrat(xsn,-xsd)*log(1d0-xs**2)
     &   -dcmplx(pisqo6)
     &   +dcmplx(ddilog(xs**2))
     &   +dcmplx(0.5d0*log(mc/ms)**2)
     &   +dcmplx(ddilog(1d0-xs*mc/ms))
     &   -impi*lnrat(ms-xs*mc,ms)  ! only if xs < 0
     &   +dcmplx(ddilog(1d0-xs*ms/mc))
     &   -impi*lnrat(mc-xs*ms,mc)  ! only if xs < 0
     &   )*xs/(ms*mc)/(1d0-xs**2)

      return
      end

      double complex function tr1f(xsq,t,mc2)
c===============================c
c                               c
c                               c
c            |                  c
c            |                  c
c            |                  c
c           / \                 c
c          /   \                c
c         /     \               c
c        x=======\\             c
c       x         \\            c
c      x           \\           c
c                               c
c    t                          c
c                               c
c===============================c
      implicit none
      include 'constants.f'
      double precision xsq,t,mc2,ddilog
      double complex lnrat,dd1

c      mc=dsqrt(mc2)
c      mu=dsqrt(xsq)

      if (t.gt.mc2) then
      dd1=-impi*dlog(1d0+mc2/(t-mc2))
      else
      dd1=czip
      endif
      dd1=dd1+dcmplx(ddilog(1d0-mc2/(mc2-t)))

      tr1f= (
c     &       -epinv**2
c     &       -2d0*epinv*lnrat(mc*mu,mc2-t)
c     &       -2d0*lnrat(mu*mc,mc2-t)**2
     &       +lnrat(mc2,mc2-t)**2
     &       +2d0*dd1
     &       -dcmplx(pisqo6)
     &     )/2d0/(mc2-t)

      return
      end


      double complex function tr2f(xsq,s,ms2,mc2)
c===============================c
c                               c
c         s                     c
c            $                  c
c            $                  c
c            $                  c
c           x \\                c
c          x   \\               c
c         x     \\              c
c        x-------\\             c
c       x         \\            c
c      x           \\           c
c                               c
c    ms2           mc2          c
c                               c
c===============================c
      implicit none
      include 'constants.f'

      double precision xsq,s,ms2,mc2,ddilog
      double precision ms,mc,xs
      double precision xsn,xsd
      double complex lnrat
      
      ms=dsqrt(ms2)
      mc=dsqrt(mc2)

      xsn=(1d0-dsqrt(1d0-4d0*ms*mc/(s-(ms-mc)**2)))
      xsd=(1d0+dsqrt(1d0-4d0*ms*mc/(s-(ms-mc)**2)))

      xs=-xsn/xsd

      tr2f= + (
c     &   -epinv*lnrat(xsn,-xsd)
c     &   -lnrat(xsn,-xsd)*lnrat(xsq,mc*ms)
     &   -0.5d0*lnrat(xsn,-xsd)**2
     &   +2d0*lnrat(xsn,-xsd)*log(1d0-xs**2)
     &   -dcmplx(pisqo6)
     &   +dcmplx(ddilog(xs**2))
     &   +dcmplx(0.5d0*log(mc/ms)**2)
     &   +dcmplx(ddilog(1d0-xs*mc/ms))
     &   -impi*lnrat(ms-xs*mc,ms)  ! only if xs < 0
     &   +dcmplx(ddilog(1d0-xs*ms/mc))
     &   -impi*lnrat(mc-xs*ms,mc)  ! only if xs < 0
     &   )*xs/(ms*mc)/(1d0-xs**2)

      return
      end

      double complex function tr4(s,mc2)
      implicit none
c 
c                          $  s
c                          $
c                          $
c                          y\
c                    mc2  y  \ zip
c                        y    \
c             zip  -----yyyyyyyyyyy  mc2  
c
c
c
      double precision s,ss,mc2,ddilog
      double complex lnrat

      if (s.lt.0d0) then
      ss = -s
      else
      ss = s
      endif

      tr4= -(lnrat(ss,mc2)*lnrat(mc2-s,mc2)+dcmplx(ddilog(1d0-s/mc2)))
     &     /(mc2-s)

      return
      end


      double complex function tr5(msc,qsq,mc2,ms2)
      implicit none
      include 'constants.f'

c 
c                          | 
c                          |
c                          |
c                          yy
c                    mc2  y  y
c                        y    y
c             msc   $$$$$xxxxxx*****  qsq
c
c                          ms2
c
      double precision msc,qsq,mc2,ms2,ddilog,arg1,arg2,arg3,arg4
      double complex logs

      arg1=(2*msc)/(mc2-ms2+msc-Sqrt(-4*mc2*msc+(mc2-ms2+msc)**2))
      arg2=(2*msc)/(mc2-ms2+msc+Sqrt(-4*mc2*msc+(mc2-ms2+msc)**2))
      arg3=(2*qsq)/(mc2-ms2+qsq-Sqrt(-4*mc2*qsq+(mc2-ms2+qsq)**2))
      arg4=(2*qsq)/(mc2-ms2+qsq+Sqrt(-4*mc2*qsq+(mc2-ms2+qsq)**2))

      if (qsq.gt.((sqrt(ms2)+sqrt(mc2))**2)) then
      logs=
     & -impi*log(arg3)
     & +impi*log(arg4)
      else
      logs=czip
      endif

      tr5= - (
     & +dcmplx(ddilog(arg1))
     & +impi*log(arg1)
     & +dcmplx(ddilog(arg2))
     & -impi*log(arg2)
     & -dcmplx(ddilog(arg3))
     & -dcmplx(ddilog(arg4))
     & +logs
     &    )/(msc-qsq)

c      write(6,*) 'arg1=',arg1
c      write(6,*) 'arg2=',arg2
c      write(6,*) 'arg3=',arg3
c      write(6,*) 'arg4=',arg4

      return
      end


      double complex function B0x(xsq,s,mc2,ms2)
c 
c                           
c                         mc2   
c                          
c                         yyyy
c                        y    y
c                   *****x    x*****  s
c                         xxxx
c                          
c                          
c                          ms2
c
      implicit none
      include 'constants.f'
      include 'epinv.f'
      double precision xsq,mc2,ms2,s,mc,ms
      double complex lam2,lam,lll,prod

      mc=sqrt(mc2)
      ms=sqrt(ms2)
      lam2=dcmplx(s**2+ms2**2+mc2**2-2*s*ms2-2*s*mc2-2*ms2*mc2)

      if   ( s .lt. 0d0 ) then
      lam=sqrt(lam2)
      lll=log((dcmplx(s-ms2-mc2)-lam)/(dcmplx(s-ms2-mc2)+lam))
      prod=lam*lll
      elseif ( ( s .gt. 0d0 ) .and.
     &         ( s .lt. (ms+mc)**2 ) ) then
      lam=sqrt(lam2)
      lll=log((dcmplx(s-ms2-mc2)-lam)/(dcmplx(s-ms2-mc2)+lam))+2*impi
      prod=lam*lll
      else
      lam=sqrt(lam2)
      lll=log((dcmplx(s-ms2-mc2)-lam)/(dcmplx(s-ms2-mc2)+lam))+2*impi
      prod=lam*lll
      endif

      B0x=
     &     + dcmplx(epinv
     &     + dlog(xsq/mc/ms)
     &     + 2d0)
     &     +(+dcmplx((ms2-mc2)*dlog(mc2/ms2))
     &       +prod
     &      )/2/s

      return
      end



      double complex function B0xf(xsq,s,mc2,ms2)
c 
c                           
c                         mc2   
c                          
c                         yyyy
c                        y    y
c                   *****x    x*****  s
c                         xxxx
c                          
c                          
c                          ms2
c
      implicit none
      include 'constants.f'
c      include 'epinv.f'
      double precision xsq,mc2,ms2,s,mc,ms
      double complex lam2,lam,lll,prod

      mc=sqrt(mc2)
      ms=sqrt(ms2)
      lam2=dcmplx(s**2+ms2**2+mc2**2-2*s*ms2-2*s*mc2-2*ms2*mc2)

      if   ( s .lt. 0d0 ) then
      lam=sqrt(lam2)
      lll=log((dcmplx(s-ms2-mc2)-lam)/(dcmplx(s-ms2-mc2)+lam))
      prod=lam*lll
      elseif ( ( s .gt. 0d0 ) .and.
     &         ( s .lt. (ms+mc)**2 ) ) then
      lam=sqrt(lam2)
      lll=log((dcmplx(s-ms2-mc2)-lam)/(dcmplx(s-ms2-mc2)+lam))+2*impi
      prod=lam*lll
      else
      lam=sqrt(lam2)
      lll=log((dcmplx(s-ms2-mc2)-lam)/(dcmplx(s-ms2-mc2)+lam))+2*impi
      prod=lam*lll
      endif

      B0xf=
c     &     + epinv
c     &     + 2d0
     &     + dcmplx(dlog(xsq/mc/ms))
     &     +(+dcmplx((ms2-mc2)*dlog(mc2/ms2))
     &       +prod
     &      )/2/s

      return
      end


      double complex function tr3(s,qsq,mc2,ms2)
      implicit none
      include 'constants.f'
c     Triangle function with one internal mass and two spacelike
c     offshell lines and one timelike of mass k1sq
c 
c                          $  s
c                          $
c                          $
c                          %\
c                    ms2  %  \
c                        %    \
c             qsq   *****==========  mc2  
c
c
c
c   Adapted from
c   %\cite{Denner:kt}
c   \bibitem{Denner:kt}
c   A.~Denner,
c   %``Techniques For Calculation Of Electroweak Radiative 
c   Corrections At The One Loop Level And Results For W Physics At Lep-200,''
c   Fortsch.\ Phys.\  {\bf 41}, 307 (1993).
c   %%CITATION = FPYKA,41,307;%%

      integer i,j,k
      double precision lambda,alpha,y0(0:2),xp(0:2),xm(0:2),
     & yip(0:2),yim(0:2),alphai(0:2),x,y,z,s,qsq,mc2,ms2
      double precision psq(0:2,0:2),msq(0:2),I3mer,I3mei

      double precision ddilog,polylog
      double precision aa,bb,cc,dd

      integer,parameter:: jj(0:2)=(/1,2,0/)
      integer,parameter:: kk(0:2)=(/2,0,1/)

      polylog(y,z)=ddilog((y-1d0)/z)-ddilog(y/z)
      lambda(x,y,z)=sqrt(x**2+y**2+z**2-2d0*(x*y+y*z+z*x))

      psq(0,1)=qsq
      psq(1,2)=mc2
      psq(2,0)=s
      msq(0)=ms2
      msq(1)=mc2
      msq(2)=0d0
      alpha=psq(0,1)**2+psq(1,2)**2+psq(2,0)**2
     & -2d0*(psq(0,1)*psq(1,2)+psq(1,2)*psq(2,0)+psq(2,0)*psq(0,1))
      if (alpha .lt. 0d0) then
      write(6,*) 'I3me:Formula not implemented for alpha imaginary'
      stop
      endif
      alpha=sqrt(alpha)
      I3mer=0d0
      I3mei=0d0

      do i=0,2
c Do some checks for numerical stability.
         if (i.eq.2 .and. -qsq/max(mc2,ms2).lt.10d-4) then
            if (mc2.gt.ms2*10) then
               yip(i) = (mc2*(ms2 - s)**2)/((mc2 - ms2)*(mc2 - s)**2) + 
     &              (mc2*qsq*(ms2 - s)**2*(3*mc2**3 - 3*mc2**2*ms2 
     &              + mc2*ms2*(ms2 - 4*s) + 
     &              ms2*s*(2*ms2 + s)))/((mc2 - ms2)**3*(mc2 - s)**4)
               yim(i) =  (mc2 - ms2)/qsq - (ms2/(mc2 - ms2)) +
     &              (mc2*(mc2 - ms2))/(mc2 - s)**2 + 
     &              (2*mc2)/(-mc2 + s)
            elseif (ms2.gt.mc2*10) then
               yip(i) =  (mc2 - ms2)/qsq - (ms2/(mc2 - ms2)) +
     &              (mc2*(mc2 - ms2))/(mc2 - s)**2 + 
     &              (2*mc2)/(-mc2 + s)
               yim(i)= (mc2*(ms2 - s)**2)/((mc2 - ms2)*(mc2 - s)**2) + 
     &              (mc2*qsq*(ms2 - s)**2*(3*mc2**3 - 3*mc2**2*ms2 +
     &              mc2*ms2*(ms2 - 4*s) + 
     &              ms2*s*(2*ms2 + s)))/((mc2 - ms2)**3*(mc2 - s)**4)
            else
c Also this approximation does not work fine.
c Let's hope that the original numbers are okay               
            endif
         else
            j=jj(i)
            k=kk(i)
            alphai(i)=lambda(psq(j,k),msq(j),msq(k))
            y0(i)=0.5d0/alpha/psq(j,k)
     &           *(psq(j,k)*(psq(j,k)-psq(k,i)-psq(i,j)
     &           +2d0*msq(i)-msq(j)-msq(k))
     &           -(psq(k,i)-psq(i,j))*(msq(j)-msq(k)))
     &           +0.5d0*(psq(j,k)-msq(j)+msq(k))/psq(j,k)
            xp(i)=0.5d0/psq(j,k)*(psq(j,k)-msq(j)+msq(k)+alphai(i))
            xm(i)=0.5d0/psq(j,k)*(psq(j,k)-msq(j)+msq(k)-alphai(i))
            yip(i)=y0(i)-xp(i)
            yim(i)=y0(i)-xm(i)
            yip(i)=y0(i)-xp(i)
            yim(i)=y0(i)-xm(i)
         endif
         I3mer=I3mer+polylog(y0(i),yip(i))+polylog(y0(i),yim(i))
      enddo 

      if (qsq.gt.((sqrt(ms2)+sqrt(mc2))**2)) then
      aa=(y0(2)-1d0)/yim(2)
      bb=y0(1)/yip(1)
      cc=(y0(2)-1d0)/yip(2)
      dd=y0(1)/yim(1)
      I3mei=pi*dlog(aa*dd/bb/cc)
      else
      I3mei=0d0
      endif

C---- Minus sign for (-1)^n
      tr3=Dcmplx(I3mer,I3mei)/alpha

      return
      end

