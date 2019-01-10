      double precision function vol_wt(W)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'process.f'
      logical first
      double precision W,eval_int,volsave
      data first/.true./
      save first,volsave

      if (case .eq. 'vlchwh') then
c--- volume for production of Wtg
      if (first) then
        volsave=eval_int(W,mt,wmass)/twopi**5
        first=.false.
      endif
      vol_wt=volsave
      vol_wt=vol_wt*(wmass*wwidth/2d0)
      vol_wt=vol_wt/8d0/pi           
      vol_wt=vol_wt*(mt*twidth/2d0)
      vol_wt=vol_wt*pi*(mt**2-wmass**2)/2d0/mt**2
     .    /twopi**2
      vol_wt=vol_wt*(wmass*wwidth/2d0)
      vol_wt=vol_wt/8d0/pi
      return
      endif
      
c--- volume for W+t production      
      vol_wt=pi*dsqrt((W-mt**2-wmass**2)**2-4d0*mt**2*wmass**2)/2d0/W
     .      /twopi**2
      vol_wt=vol_wt*(wmass*wwidth/2d0)
      vol_wt=vol_wt/8d0/pi
      
      if (case .eq. 'vlchwn') return

      if (case .eq. 'vlchwt') then
c--- extra volume for decay t->Wb
      vol_wt=vol_wt*(mt*twidth/2d0)
      vol_wt=vol_wt*pi*(mt**2-wmass**2)/2d0/mt**2
     .    /twopi**2
      vol_wt=vol_wt*(wmass*wwidth/2d0)
      vol_wt=vol_wt/8d0/pi
      endif
      
      if (case .eq. 'vlchwg') then
c--- extra volume for decay t->Wbg
      vol_wt=vol_wt*(mt*twidth/2d0)
      vol_wt=vol_wt*pisq/8d0/mt**2/twopi**5
     .   *(mt**4-wmass**4-2d0*mt**2*wmass**2*dlog(mt**2/wmass**2))
      vol_wt=vol_wt*(wmass*wwidth/2d0)
      vol_wt=vol_wt/8d0/pi
      endif
      
      return
      end
      
      double precision function eval_int(s,mt,mw)
      implicit none
      integer j,n
      double precision s,mt,mw,dx,fun_int,part,min

      n=1000000
      dx=(s-(mt+mw)**2)/dfloat(n)
      min=(mt+mw)**2

      eval_int=0d0
      do j=1,n
      part=dx*(fun_int(s,mt,mw,min+dfloat(j)*dx)
     .        +fun_int(s,mt,mw,min+dfloat(j-1)*dx))/2d0
      eval_int=eval_int+part
      enddo

      return
      end

      double precision function fun_int(s,mt,mw,x)
      implicit none
      include 'constants.f'
      double precision s,mt,mw,x,root

      root=(x-mt**2-mw**2)**2-4d0*mt**2*mw**2
      if (root .lt. 0d0) root=0d0
      
      fun_int=pi*pion4/s*(s-x)*dsqrt(root)/x

      return
      end
