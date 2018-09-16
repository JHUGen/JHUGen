      function vol_wt(W)
      implicit none
      include 'types.f'
      real(dp):: vol_wt
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'kprocess.f'
      logical:: first
      real(dp):: W,eval_int,volsave
      data first/.true./
      save first,volsave

      if (kcase==kvlchwh) then
c--- volume for production of Wtg
      if (first) then
        volsave=eval_int(W,mt,wmass)/twopi**5
        first=.false.
      endif
      vol_wt=volsave
      vol_wt=vol_wt*(wmass*wwidth/2._dp)
      vol_wt=vol_wt/8._dp/pi           
      vol_wt=vol_wt*(mt*twidth/2._dp)
      vol_wt=vol_wt*pi*(mt**2-wmass**2)/2._dp/mt**2
     &    /twopi**2
      vol_wt=vol_wt*(wmass*wwidth/2._dp)
      vol_wt=vol_wt/8._dp/pi
      return
      endif
      
c--- volume for W+t production      
      vol_wt=pi*sqrt((W-mt**2-wmass**2)**2-4._dp*mt**2*wmass**2)/2._dp/W
     &      /twopi**2
      vol_wt=vol_wt*(wmass*wwidth/2._dp)
      vol_wt=vol_wt/8._dp/pi
      
      if (kcase==kvlchwn) return

      if (kcase==kvlchwt) then
c--- extra volume for decay t->Wb
      vol_wt=vol_wt*(mt*twidth/2._dp)
      vol_wt=vol_wt*pi*(mt**2-wmass**2)/2._dp/mt**2
     &    /twopi**2
      vol_wt=vol_wt*(wmass*wwidth/2._dp)
      vol_wt=vol_wt/8._dp/pi
      endif
      
      if (kcase==kvlchwg) then
c--- extra volume for decay t->Wbg
      vol_wt=vol_wt*(mt*twidth/2._dp)
      vol_wt=vol_wt*pisq/8._dp/mt**2/twopi**5
     &   *(mt**4-wmass**4-2._dp*mt**2*wmass**2*log(mt**2/wmass**2))
      vol_wt=vol_wt*(wmass*wwidth/2._dp)
      vol_wt=vol_wt/8._dp/pi
      endif
      
      return
      end
      
      function eval_int(s,mt,mw)
      implicit none
      include 'types.f'
      real(dp):: eval_int
      
      integer:: j,n
      real(dp):: s,mt,mw,dx,fun_int,part,min

      n=1000000
      dx=(s-(mt+mw)**2)/real(n,dp)
      min=(mt+mw)**2

      eval_int=0._dp
      do j=1,n
      part=dx*(fun_int(s,mt,mw,min+real(j,dp)*dx)
     &        +fun_int(s,mt,mw,min+real(j-1,dp)*dx))/2._dp
      eval_int=eval_int+part
      enddo

      return
      end

      function fun_int(s,mt,mw,x)
      implicit none
      include 'types.f'
      real(dp):: fun_int
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: s,mt,mw,x,root

      root=(x-mt**2-mw**2)**2-4._dp*mt**2*mw**2
      if (root < 0._dp) root=0._dp
      
      fun_int=pi*pion4/s*(s-x)*sqrt(root)/x

      return
      end
