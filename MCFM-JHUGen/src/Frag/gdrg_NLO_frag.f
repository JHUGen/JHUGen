!===== functions for NLO frag set of hep-ph/980631 
!===== C. Williams Feb 14 

      double precision function GdRG_NLO_frag(x,i) 
      implicit none
      include 'scale.f' 
      include 'frag.f' 
      include 'constants.f' 
      include 'ewcouple.f'
      include 'ewcharge.f' 
      include 'qcdcouple.f' 
      double precision P1qgam,conv_P0qqDqgam,conv_P0qqP0qgam
      external P1qgam,conv_P0qqDqgam,conv_P0qqP0qgam
      double precision x
      integer i 
      double precision mu_zero 
      parameter(mu_zero=0.64d0)
      double precision lmufmu0 
      double precision DNPqgam,P0qgam
      external DNPqgam,P0qgam
      double precision aewon2pi 

      aewon2pi=esq/(fourpi*twopi) 
      lmufmu0=dlog(frag_scale**2/mu_zero**2) 
      GdRG_NLO_frag=DNPqgam(x)+lmufmu0*P0qgam(x)
     & +ason2pi*lmufmu0*P1qgam(x) 
     & +one/2d0*ason2pi*lmufmu0**2*conv_P0qqP0qgam(x) 
     & +ason2pi*lmufmu0*conv_P0qqDqgam(x) 
     
      GdRG_NLO_frag=GdRG_NLO_frag*aewon2pi*Q(i)**2

      return 
      end


      
      double precision function P1qgam(x) 
!+===== one loop splitting for q=>gam 
      implicit none 
      include 'constants.f' 
      double precision x,omx,lx,lomx,dilmx 
      double precision ddilog 
      double precision P0qgam
      external P0qgam 

      omx=one-x
      lx=dlog(x) 
      lomx=dlog(omx) 
      dilmx=ddilog(omx) 

!======piece not prop to LO split 
      
      P1qgam=-one/2d0+9d0/2d0*x+(-8d0+x/2d0)*lx+2d0*x*lomx
     &     +(one-x/2d0)*lx**2

!====== piece propto LO split 
      P1qgam=P1qgam+P0qgam(x)*(lomx**2+4d0*lx*lomx+8d0*dilmx
     &     -4d0/3d0*pi**2)

!===== Casimir 
      P1qgam=V/(2d0*xn)*P1qgam 
      return 
      end


      double precision function conv_P0qqDqgam(x) 
      implicit none 
      include 'constants.f' 
      double precision x,xi
      double precision b_cdPD,i_cdPD,int_cPD
      double precision DNPqgam
      double precision dgauss,test
      common/i_cdPD_x/xi
      external P0qgam,dgauss
      external i_cdPD
!==== b_cdPP is boundary term arising from delta(1-t) piece in P_0(q=>q) 
      xi=x
      b_cdPD=3d0/2d0*DNPqgam(x)
!====== also boundary term from g(1)*log(1-x) in plus distribution 
      b_cdPD=b_cdPD+2d0*DNPqgam(x)*dlog(one-x) 
!======
      int_cPD=dgauss(i_cdPD,x,1d0,0.001d0) 
      
      conv_P0qqDqgam=b_cdPD+int_cPD
!===== overall casimir 
      conv_P0qqDqgam=V/(2d0*xn)*conv_P0qqDqgam 
      return 
      end 
      


      double precision function conv_P0qqP0qgam(x) 
!==== function for the convoultion of P_0(q=>q) with P_0(q=>gam) i.e.
!==== int(x,1) dt/t P_0(q=>q)(t) P_0(q=>gam)(x/t)  

      implicit none 
      include 'constants.f' 
      double precision x,xi
      double precision b_cdPP,i_cdPP,int_cPP
      double precision P0qgam 
      double precision dgauss,test
      common/i_cdPP_x/xi
      external P0qgam,dgauss
      external i_cdPP
!==== b_cdPP is boundary term arising from delta(1-t) piece in P_0(q=>q) 
      xi=x
      b_cdPP=3d0/2d0*P0qgam(x)
!====== also boundary term from g(1)*log(1-x) in plus distribution 
      b_cdPP=b_cdPP+2d0*P0qgam(x)*dlog(one-x) 
!======
      int_cPP=dgauss(i_cdPP,x,1d0,0.001d0) 
      
      conv_P0qqP0qgam=b_cdPP+int_cPP
!===== overall casimir 
      conv_P0qqP0qgam=V/(2d0*xn)*conv_P0qqP0qgam 
      return 
      end 
      
      double precision function i_cdPP(t) 
      implicit none 
      include 'constants.f' 
      double precision x,t
      double precision P0qqreg,P0qgam 
      external P0qqreg,P0qgam 
      common/i_cdPP_x/x 

!===== integrand for PP convo 

      i_cdPP=((one+t**2)*P0qgam(x/t)/t-2d0*P0qgam(x))/(one-t)
     
      return 
      end
      
      double precision function P0qgam(x) 
!====== eps=>0 part of splitting function 
      implicit none
      include 'constants.f' 
      double precision x 
      
      P0qgam=(one+(one-x)**2)/x

      return 
      end

      double precision function DNPqgam(x) 
!====== fragmentation function NP  Note that factor
!===== of (alpha_ew*eq^2/twopi) extracted from 2.14
      implicit none
      include 'constants.f' 
      double precision x,omx,omxsq
      double precision P0qgam 
      external P0qgam
      omx=one-x
      omxsq=omx**2
      DNPqgam=-P0qgam(x)*dlog(omxsq)+20.8d0*omx-11.07d0

      return 
      end

      double precision function i_cdPD(t) 
      implicit none 
      include 'constants.f' 
      double precision x,t
      double precision P0qqreg,DNPqgam 
      external P0qqreg,DNPqgam 
      common/i_cdPD_x/x         !

!===== integrand for PD convo 

       i_cdPD=((one+t**2)*DNPqgam(x/t)/t-2d0*DNPqgam(x))/(one-t)
     
       return 
       end
      
