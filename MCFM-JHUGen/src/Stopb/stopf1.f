      function tr3c00f(xsqV,s,qsq,mc2,ms2)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: tr3c00f
      
      real(dp):: xsqV,s,qsq,ms2,mc2,Denc
      include 'stopf1inc.f'

      Denc = -(mc2**2+s**2+qsq**2-2*(mc2*s+mc2*qsq+s*qsq))/4._dp
      Denc = 1/Denc
      tr3c00f =
     &        (4*(cplx1(3._dp)+ B0cgsf)*s + 
     &    Denc*(ms2 - s)*(ls*(ms2 - s)*(mc2 - qsq + s) + 
     &       s*(-2*lVc*mc2 + B0cgsf*(mc2 + qsq - s) + 
     &          lVs*(mc2 - qsq + s) + 2*mc2*(-ms2 + s)*tr3Xc)))/
     &  (16._dp*s)

      return
      end



      function tr3c001f(xsqV,s,qsq,mc2,ms2)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: tr3c001f
      
      real(dp):: xsqV,s,qsq,ms2,mc2,Denc
      include 'stopf1inc.f'

      Denc = -(mc2**2+s**2+qsq**2-2*(mc2*s+mc2*qsq+s*qsq))/4._dp
      Denc = 1/Denc
      tr3c001f =
     &        -(8*(cplx1(19._dp)+ 6*B0cgsf + 3*BfunX)*s**2 + 
     &     9*Denc**2*mc2*(ms2 - s)**2*s*
     &      (ls*(ms2 - s)*(mc2 - qsq + s) + 
     &        s*(-2*lVc*mc2 + B0cgsf*(mc2 + qsq - s) + 
     &           lVs*(mc2 - qsq + s) + 2*mc2*(-ms2 + s)*tr3Xc)) + 
     &     6*Denc*(ms2 - s)*
     &      (ls*(ms2 - s)*(mc2 - qsq + s)*(ms2 + 2*s) + 
     &        s*(-cplx1(ms2*qsq)+(cplx1(3._dp)+2*B0cgsf-2*lVs)*ms2*s - 
     &           s*(-2*B0cgsf*qsq + BfunX*qsq + 2*lVs*qsq + 
     &              2*B0cgsf*s + BfunX*s - 2*lVs*s) + 
     &           mc2*(cplx1(ms2) - 6*ms2*s*tr3Xc + 
     &              s*(2*B0cgsf + BfunX - 4*lVc + 2*lVs + 6*s*tr3Xc)
     &              ))))/(288._dp*s**2)

      return
      end



      function tr3c002f(xsqV,s,qsq,mc2,ms2)
      implicit none
      include 'types.f'
      complex(dp):: tr3c002f
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: xsqV,s,qsq,ms2,mc2,Denc
      include 'stopf1inc.f'

      Denc = -(mc2**2+s**2+qsq**2-2*(mc2*s+mc2*qsq+s*qsq))/4._dp
      Denc = 1/Denc
      tr3c002f =
     &        (-8*(cplx1(8._dp)+ 3*B0cgsf + 3*BfunX)*s + 
     &    6*Denc*(ms2 - s)*(ls*(ms2 - s)**2 + 
     &       s*(-((cone + BfunX)*mc2) - cplx1(2*ms2) - 3*B0cgsf*ms2 + 
     &          3*lVs*ms2 + cplx1(qsq) + BfunX*qsq + 
     &          (-cone + B0cgsf + BfunX - lVs)*s)) + 
     &    9*Denc**2*mc2*(ms2 - s)**2*s*
     &     (-2*lVs*s + 2*ls*(-ms2 + s) + lVc*(mc2 - qsq + s) + 
     &       B0cgsf*(-mc2 + qsq + s) + 
     &       (ms2 - s)*(mc2 - qsq + s)*tr3Xc))/(288._dp*s)

      return
      end



      function tr3s00f(xsqV,t,qsq,ms2,mc2)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: tr3s00f
      
      real(dp):: xsqV,t,qsq,ms2,mc2,Dens
      include 'stopf1inc.f'

      Dens = -(ms2**2+t**2+qsq**2-2*(ms2*t+ms2*qsq+t*qsq))/4._dp
      Dens = 1/Dens
      tr3s00f =
     &        (4*(cplx1(3._dp)+ B0cgsf)*t + 
     &    Dens*(mc2 - t)*(lc*(mc2 - t)*(ms2 - qsq + t) + 
     &       t*(B0cgsf*(ms2 + qsq - t) + lVc*(ms2 - qsq + t) - 
     &          2*ms2*(lVs + mc2*tr3Xs - t*tr3Xs))))/(16._dp*t)

      return
      end



      function tr3s001f(xsqV,t,qsq,ms2,mc2)
      implicit none
      include 'types.f'
      complex(dp):: tr3s001f
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: xsqV,t,qsq,ms2,mc2,Dens
      include 'stopf1inc.f'

      Dens = -(ms2**2+t**2+qsq**2-2*(ms2*t+ms2*qsq+t*qsq))/4._dp
      Dens = 1/Dens
      tr3s001f =
     &        -(8*(cplx1(19._dp)+ 6*B0cgsf + 3*BfunX)*t + 
     &     9*Dens**2*ms2*(mc2 - t)**2*t*
     &      (-2*lc*(mc2 - t) - 2*lVc*t + B0cgsf*(-ms2 + qsq + t) + 
     &        (ms2 - qsq + t)*(lVs + mc2*tr3Xs - t*tr3Xs)) + 
     &     6*Dens*(mc2 - t)*
     &      (lc*(mc2 - t)*(mc2 + 3*ms2 - 3*qsq + 2*t) + 
     &        t*((-cplx1(3._dp)+ BfunX)*ms2 + 
     &           B0cgsf*(-mc2 + ms2 + 3*qsq - 2*t) + 
     &           (cone + BfunX)*(qsq - t) + 
     &           lVc*(mc2 + 3*ms2 - 3*qsq + 2*t) - 
     &           2*ms2*(2*lVs + 3*(mc2 - t)*tr3Xs))))/(288._dp*t)

      return
      end



      function tr3s002f(xsqV,t,qsq,ms2,mc2)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: tr3s002f
      
      real(dp):: xsqV,t,qsq,ms2,mc2,Dens
      include 'stopf1inc.f'

      Dens = -(ms2**2+t**2+qsq**2-2*(ms2*t+ms2*qsq+t*qsq))/4._dp
      Dens = 1/Dens
      tr3s002f =
     &        (-8*(cplx1(8._dp)+ 3*B0cgsf + 3*BfunX)*t**2 + 
     &    6*Dens*(mc2 - t)*(lc*(mc2 - t)**2*(ms2 - qsq + t) + 
     &       t*(((cplx1(2._dp) + B0cgsf - BfunX - lVc)*ms2 - 
     &             (B0cgsf + BfunX - lVc)*(qsq - t))*t + 
     &          cplx1(mc2*(ms2 - qsq + t)))) + 
     &    9*Dens**2*ms2*(mc2 - t)**2*t*
     &     (lc*(mc2 - t)*(ms2 - qsq + t) + 
     &       t*(B0cgsf*(ms2 + qsq - t) + lVc*(ms2 - qsq + t) - 
     &          2*ms2*(lVs + mc2*tr3Xs - t*tr3Xs))))/(288._dp*t**2)

      return
      end



      function Bfun(xsqV,qsq,ms2,mc2)
      implicit none
      include 'types.f'
      complex(dp):: Bfun
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: xsqV,qsq,ms2,mc2
      include 'stopf1inc.f'

      Bfun =
     &((-cone - B0cgsf + lVc)*mc2 + (cone + B0cgsf - lVs)*ms2)/qsq

      return
      end



