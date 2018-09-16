!===== functions for NLO frag set of hep-ph/980631 
!===== C. Williams Feb 14 

      function GdRG_NLO_frag(x,i) 
      implicit none
      include 'types.f'
      real(dp):: GdRG_NLO_frag
      
      include 'scale.f' 
      include 'frag.f' 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'ewcouple.f'
      include 'ewcharge.f' 
      include 'qcdcouple.f' 
      real(dp):: P1qgam,conv_P0qqDqgam,conv_P0qqP0qgam
      external P1qgam,conv_P0qqDqgam,conv_P0qqP0qgam
      real(dp):: x
      integer:: i 
      real(dp):: mu_zero 
      parameter(mu_zero=0.64_dp)
      real(dp):: lmufmu0 
      real(dp):: DNPqgam,P0qgam
      external DNPqgam,P0qgam
      real(dp):: aewon2pi 

      aewon2pi=esq/(fourpi*twopi) 
      lmufmu0=log(frag_scale**2/mu_zero**2) 
      GdRG_NLO_frag=DNPqgam(x)+lmufmu0*P0qgam(x)
     & +ason2pi*lmufmu0*P1qgam(x) 
     & +one/2._dp*ason2pi*lmufmu0**2*conv_P0qqP0qgam(x) 
     & +ason2pi*lmufmu0*conv_P0qqDqgam(x) 
     
      GdRG_NLO_frag=GdRG_NLO_frag*aewon2pi*Q(i)**2

      return 
      end


      
      function P1qgam(x) 
      implicit none
      include 'types.f'
      real(dp):: P1qgam
!+===== one loop splitting for q=>gam 
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      real(dp):: x,omx,lx,lomx,dilmx 
      real(dp):: ddilog 
      real(dp):: P0qgam
      external P0qgam 

      omx=one-x
      lx=log(x) 
      lomx=log(omx) 
      dilmx=ddilog(omx) 

!======piece not prop to LO split 
      
      P1qgam=-one/2._dp+9._dp/2._dp*x+(-8._dp+x/2._dp)*lx+2._dp*x*lomx
     &     +(one-x/2._dp)*lx**2

!====== piece propto LO split 
      P1qgam=P1qgam+P0qgam(x)*(lomx**2+4._dp*lx*lomx+8._dp*dilmx
     &     -4._dp/3._dp*pi**2)

!===== Casimir 
      P1qgam=V/(2._dp*xn)*P1qgam 
      return 
      end


      function conv_P0qqDqgam(x) 
      implicit none
      include 'types.f'
      real(dp):: conv_P0qqDqgam
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      real(dp):: x,xi
      real(dp):: b_cdPD,i_cdPD,int_cPD
      real(dp):: DNPqgam
      real(dp):: dgauss,test
      common/i_cdPD_x/xi
      external P0qgam,dgauss
      external i_cdPD
!==== b_cdPP is boundary term arising from delta(1-t) piece in P_0(q=>q) 
      xi=x
      b_cdPD=3._dp/2._dp*DNPqgam(x)
!====== also boundary term from g(1)*log(1-x) in plus distribution 
      b_cdPD=b_cdPD+2._dp*DNPqgam(x)*log(one-x) 
!======
      int_cPD=dgauss(i_cdPD,x,1._dp,0.001_dp) 
      
      conv_P0qqDqgam=b_cdPD+int_cPD
!===== overall casimir 
      conv_P0qqDqgam=V/(2._dp*xn)*conv_P0qqDqgam 
      return 
      end 
      


      function conv_P0qqP0qgam(x) 
      implicit none
      include 'types.f'
      real(dp):: conv_P0qqP0qgam
!==== function for the convoultion of P_0(q=>q) with P_0(q=>gam) i.e.
!==== int(x,1) dt/t P_0(q=>q)(t) P_0(q=>gam)(x/t)  

       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      real(dp):: x,xi
      real(dp):: b_cdPP,i_cdPP,int_cPP
      real(dp):: P0qgam 
      real(dp):: dgauss,test
      common/i_cdPP_x/xi
      external P0qgam,dgauss
      external i_cdPP
!==== b_cdPP is boundary term arising from delta(1-t) piece in P_0(q=>q) 
      xi=x
      b_cdPP=3._dp/2._dp*P0qgam(x)
!====== also boundary term from g(1)*log(1-x) in plus distribution 
      b_cdPP=b_cdPP+2._dp*P0qgam(x)*log(one-x) 
!======
      int_cPP=dgauss(i_cdPP,x,1._dp,0.001_dp) 
      
      conv_P0qqP0qgam=b_cdPP+int_cPP
!===== overall casimir 
      conv_P0qqP0qgam=V/(2._dp*xn)*conv_P0qqP0qgam 
      return 
      end 
      
      function i_cdPP(t) 
      implicit none
      include 'types.f'
      real(dp):: i_cdPP
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      real(dp):: x,t
      real(dp):: P0qqreg,P0qgam 
      external P0qqreg,P0qgam 
      common/i_cdPP_x/x 

!===== integrand for PP convo 

      i_cdPP=((one+t**2)*P0qgam(x/t)/t-2._dp*P0qgam(x))/(one-t)
     
      return 
      end
      
      function P0qgam(x) 
      implicit none
      include 'types.f'
      real(dp):: P0qgam
!====== eps=>0 part of splitting function 
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      real(dp):: x 
      
      P0qgam=(one+(one-x)**2)/x

      return 
      end

      function DNPqgam(x) 
      implicit none
      include 'types.f'
      real(dp):: DNPqgam
!====== fragmentation function NP  Note that factor
!===== of (alpha_ew*eq^2/twopi) extracted from 2.14
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      real(dp):: x,omx,omxsq
      real(dp):: P0qgam 
      external P0qgam
      omx=one-x
      omxsq=omx**2
      DNPqgam=-P0qgam(x)*log(omxsq)+20.8_dp*omx-11.07_dp

      return 
      end

      function i_cdPD(t) 
      implicit none
      include 'types.f'
      real(dp):: i_cdPD
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      real(dp):: x,t
      real(dp):: P0qqreg,DNPqgam 
      external P0qqreg,DNPqgam 
      common/i_cdPD_x/x         !

!===== integrand for PD convo 

       i_cdPD=((one+t**2)*DNPqgam(x/t)/t-2._dp*DNPqgam(x))/(one-t)
     
       return 
       end
      
