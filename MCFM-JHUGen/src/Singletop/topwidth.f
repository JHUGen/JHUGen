c      program main
c      
c      real(dp):: mt,mw,tw,topwidth
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'qcdcouple.f'
c      
c      ason2pi=0.119_dp/2._dp/pi
c      
c      mw=82._dp
c      mt=180._dp 
c      tw=topwidth(mt,mw)
c      
c      write(6,*) 'twidth/born ',tw
c      
c      end


      function topwidth(mt,mw)
      implicit none
      include 'types.f'
      real(dp):: topwidth
      
C     alpha_s correction to top width
C     Formula taken from
C %\cite{Kuhn:1996ug}
c \bibitem{Kuhn:1996ug}
c J.~H.~Kuhn,
c %``Theory of top quark production and decay,''
c arXiv:hep-ph/9707321.
c %%CITATION = HEP-PH 9707321;%%
      real(dp):: y,f0,f1,mt,mw,omy,ddilog
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      
      y=(mw/mt)**2
      omy=1._dp-y
      F0=2._dp*omy**2*(1._dp+2._dp*y)
      F1=F0*(pisq+2._dp*(ddilog(y)-ddilog(omy)))
     & +4._dp*y*(1._dp-y*(1._dp+2._dp*y))*log(y)
     & +2._dp*omy**2*(5._dp+4._dp*y)*log(omy)
     & -omy*(5._dp+y*(9._dp-6._dp*y))
     
c      F1=F0*(2._dp/3._dp*pisq-2.5_dp-3._dp*y+4.5_dp*y**2
c     & -3._dp*y**2*log(y))
c      write(6,*) F1/F0
     
      topwidth=1._dp-cf*ason2pi*F1/F0
      
c      write(6,*) 'dQCD ',-cf/2._dp*F1/F0
      return       
      end

c      function virtbit(mt,mw)
c      implicit none
c      include 'types.f'
c      real(dp):: virtbit
c      
C     alpha_s correction to top width, virtual+dipole integration only
c      real(dp):: y,f0,d1,mt,mw,omy,ddilog
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'qcdcouple.f'

c      y=(mw/mt)**2
c      omy=1._dp-y
c      F0=2._dp*omy**2*(1._dp+2._dp*y)


c      d1=-F0*(pisq+2._dp*(ddilog(y)-ddilog(omy)))
c     &  +log(omy)*(-10._dp+12._dp*y+6._dp*y**2-8._dp*y**3)
c     &  +y*log(y)*(-6._dp-5._dp*y+14._dp*y**2)
c     &  +5._dp/2._dp+y-13._dp/2._dp*y**2+3._dp*y**3

c      virtbit=1._dp+cf*ason2pi*d1/F0
c      return
c      end

c      function realbit(mt,mw)
c      implicit none
c      include 'types.f'
c      real(dp):: realbit
c      
C     alpha_s correction to top width
C     Formula taken from
c      real(dp):: mt,mw,rsq,omrsq,wlog,rlog,ddilog,ct,Kfun
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'qcdcouple.f'
c      include 'epinv.f'
c      include 'epinv2.f'

c      rsq=(mw/mt)**2
c      omrsq=1._dp-rsq
c      wlog=log(omrsq)
c      rlog=log(rsq)

c      ct=epinv2*epinv
c     & +epinv*(2.5_dp-2._dp*wlog)
c     & +25._dp/4._dp+0.5_dp*(1._dp/omrsq**2-8._dp/omrsq+7._dp)*rlog
c     & +0.5_dp/omrsq+2._dp*ddilog(omrsq)-5._dp*pisqo6
c     & -5._dp*wlog+2._dp*wlog**2

c      Kfun=1._dp/rsq*wlog
       
c      ct=ct-epinv*epinv2
c     & -epinv*(2.5_dp-2._dp*wlog)
c     & -0.5_dp*(11._dp)-pisqo6-2._dp*ddilog(rsq)
c     &  +3._dp*wlog-2._dp*wlog**2-Kfun
c
c
c      realbit=1._dp+cf*ason2pi*ct
c      
c      return
c      end
