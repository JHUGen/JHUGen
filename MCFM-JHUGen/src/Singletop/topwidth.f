c      program main
c      implicit none
c      double precision mt,mw,tw,topwidth
c      include 'constants.f'
c      include 'qcdcouple.f'
c      
c      ason2pi=0.119d0/2d0/pi
c      
c      mw=82d0
c      mt=180d0 
c      tw=topwidth(mt,mw)
c      
c      write(6,*) 'twidth/born ',tw
c      
c      end


      double precision function topwidth(mt,mw)
      implicit none
C     alpha_s correction to top width
C     Formula taken from
C %\cite{Kuhn:1996ug}
c \bibitem{Kuhn:1996ug}
c J.~H.~Kuhn,
c %``Theory of top quark production and decay,''
c arXiv:hep-ph/9707321.
c %%CITATION = HEP-PH 9707321;%%
      double precision y,f0,f1,mt,mw,omy,ddilog
      include 'constants.f'
      include 'qcdcouple.f'
      
      y=(mw/mt)**2
      omy=1d0-y
      F0=2d0*omy**2*(1d0+2d0*y)
      F1=F0*(pisq+2d0*(ddilog(y)-ddilog(omy)))
     . +4d0*y*(1d0-y*(1d0+2d0*y))*log(y)
     . +2d0*omy**2*(5d0+4d0*y)*log(omy)
     . -omy*(5d0+y*(9d0-6d0*y))
     
c      F1=F0*(2d0/3d0*pisq-2.5d0-3d0*y+4.5d0*y**2
c     . -3d0*y**2*dlog(y))
c      write(6,*) F1/F0
     
      topwidth=1d0-cf*ason2pi*F1/F0
      
c      write(6,*) 'dQCD ',-cf/2d0*F1/F0
      return       
      end

c      double precision function virtbit(mt,mw)
c      implicit none
C     alpha_s correction to top width, virtual+dipole integration only
c      double precision y,f0,d1,mt,mw,omy,ddilog
c      include 'constants.f'
c      include 'qcdcouple.f'

c      y=(mw/mt)**2
c      omy=1d0-y
c      F0=2d0*omy**2*(1d0+2d0*y)


c      d1=-F0*(pisq+2d0*(ddilog(y)-ddilog(omy)))
c     .  +log(omy)*(-10d0+12d0*y+6d0*y**2-8d0*y**3)
c     .  +y*log(y)*(-6d0-5d0*y+14d0*y**2)
c     .  +5d0/2d0+y-13d0/2d0*y**2+3d0*y**3

c      virtbit=1d0+cf*ason2pi*d1/F0
c      return
c      end

c      double precision function realbit(mt,mw)
c      implicit none
C     alpha_s correction to top width
C     Formula taken from
c      double precision mt,mw,rsq,omrsq,wlog,rlog,ddilog,ct,Kfun
c      include 'constants.f'
c      include 'qcdcouple.f'
c      include 'epinv.f'
c      include 'epinv2.f'

c      rsq=(mw/mt)**2
c      omrsq=1d0-rsq
c      wlog=log(omrsq)
c      rlog=log(rsq)

c      ct=epinv2*epinv
c     . +epinv*(2.5d0-2d0*wlog)
c     . +25d0/4d0+0.5d0*(1d0/omrsq**2-8d0/omrsq+7d0)*rlog
c     . +0.5d0/omrsq+2d0*ddilog(omrsq)-5d0*pisqo6
c     . -5d0*wlog+2d0*wlog**2

c      Kfun=1d0/rsq*wlog
       
c      ct=ct-epinv*epinv2
c     . -epinv*(2.5d0-2d0*wlog)
c     . -0.5d0*(11d0)-pisqo6-2d0*ddilog(rsq)
c     .  +3d0*wlog-2d0*wlog**2-Kfun
c
c
c      realbit=1d0+cf*ason2pi*ct
c      
c      return
c      end
