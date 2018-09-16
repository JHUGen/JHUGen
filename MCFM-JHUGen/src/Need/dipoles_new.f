************************************************************************
*     Author: J. M. Campbell                                           *
*     August, 1999  (updated April, 2001)                              *
*                                                                      *
*     Comments added October 15th 2001.                                *
*                                                                      *
*     Revised by R.K. Ellis, November 9th, 16th 2001.                  *
*                                                                      *
*     Routines which return various pieces of the integrated           *
*     subtraction terms, used in both _v and _z routines               *
************************************************************************

************************************************************************
*                                                                      *
*     The labelling of the routines is as follows:                     *
*     The collinear pair is assumed to be incoming,                    *
*     so a reversal has to be made for the final state cases           *
*                                                                      *
*              -------->------------>--------                          *
*                  j        /         i                                *
*                          /                                           *
*                         /                                            *
*                                                                      *
*                represented by {ii/if}_ij                             *
*                                                                      *
************************************************************************

***********************************************************************
*************************** INITIAL-INITIAL ***************************
***********************************************************************

***************************** Quark-Quark *****************************
      function iin_qq(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: iin_qq
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx,Pqqreg,alfax
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--TH-V
c-- Id,aqq=
c--  [delta(1-x)]*(epinv*(epinv-L)+1/2*L^2+3/2*epinv-[pi]^2/6)
c--  +(1-x)-(1+x)*(L+2*[ln(1-x)])-(1+x^2)*[ln(x)]/[1-x]
c--  +4*[ln(1-x)/(1-xp)]+2*L/[1-xp]

      
      if (vorz == 1) then
        iin_qq=epinv*(epinv2-L)+0.5_dp*L**2+1.5_dp*epinv-pisqo6
     &   -epinv*1.5_dp
        if (scheme == 'tH-V') then
           return
        elseif (scheme == 'dred') then
           iin_qq=iin_qq-half
           return
        endif
      endif
      
      omx=one-x
      lomx=log(omx)
      lx=log(x)
      
      if (vorz == 2) then
        Pqqreg=-one-x
        iin_qq=omx+Pqqreg*(two*lomx+L-epinv)-(one+x**2)/omx*lx
        alfax=alfa/omx
        if (alfax < 1._dp) iin_qq=iin_qq+(two/omx+Pqqreg)*log(alfax)
        return
      endif
      
      iin_qq=two/omx*(two*lomx+L-epinv)
      
      return
      end

***************************** Quark-Gluon *****************************
      function iin_qg(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: iin_qg
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx,Pqgreg,alfax
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--TH-V
c-- Id,aqg=(1-2*x*(1-x))*(-[ln(x)]+L+2*[ln(1-x)])+2*x*(1-x)
      iin_qg=0._dp
      if ((vorz == 1) .or. (vorz == 3)) return
      
      omx=one-x
      lomx=log(omx)
      lx=log(x)
      
      if (vorz == 2) then
        Pqgreg=one-two*x*omx
        iin_qg=Pqgreg*(two*lomx-lx+L-epinv)+two*x*omx
        alfax=alfa/omx
        if (alfax < 1._dp) iin_qg=iin_qg+Pqgreg*log(alfax)
      endif
      return
      end
      
***************************** Gluon-Quark *****************************
      function iin_gq(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: iin_gq
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx,Pgqreg,alfax
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial quark-quark (--> gluon) antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--TH-V
c-- Id,agq=(1+(1-x)^2)/x*(-[ln(x)]+L+2*[ln(1-x)])+x

      
      iin_gq=0._dp
      if ((vorz == 1) .or. (vorz == 3)) return
      
      omx=one-x
      lomx=log(omx)
      lx=log(x)
      
      if (vorz == 2) then
        Pgqreg=(one+omx**2)/x
        iin_gq=Pgqreg*(two*lomx-lx+L-epinv)+x
        alfax=alfa/omx
        if (alfax < 1._dp) iin_gq=iin_gq+Pgqreg*log(alfax)
        return
      endif

      return
      end

***************************** Gluon-Gluon *****************************
      function iin_gg(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: iin_gg
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx,Pggreg,alfax
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--TH-V
c-- Id,agg=(epinv*(epinv-L)+1/2*L^2+epinv*11/6-[pi]^2/6
c--  -nf/3/xn*epinv)*[delta(1-x)]
c--  -2*[ln(x)]/[1-x]
c--  +2*(-1+x*(1-x)+(1-x)/x)*(-[ln(x)]+L+2*[ln(1-x)])
c--  +(4*[ln(1-x)/(1-xp)]+2*L/[1-xp])
      
      if (vorz == 1) then
        iin_gg=epinv*(epinv2-L)+half*L**2-pisqo6
     &       +(11._dp-two*nf/xn)/6._dp*epinv
     &       -(11._dp-two*nf/xn)/6._dp*epinv
        if (scheme == 'tH-V') then
        return
        elseif (scheme == 'dred') then
        iin_gg=iin_gg-1._dp/6._dp
        return
        endif
      endif
      
      omx=one-x
      lomx=log(omx)
      
      if (vorz == 2) then
        Pggreg=omx/x+x*omx-one
        lx=log(x)
        iin_gg=two*Pggreg*(two*lomx-lx+L-epinv)-two*lx/omx
        alfax=alfa/omx
        if (alfax < 1._dp) iin_gg=iin_gg+Pggreg*log(alfax)
        return
      endif
      
      iin_gg=two*(two*lomx+L-epinv)/omx
      
      return
      end

***********************************************************************
**************************** INITIAL-FINAL ****************************
***********************************************************************

***************************** Quark-Quark *****************************
      function ifn_qq(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: ifn_qq
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx,ltmx,Pqqpr,Pqqreg
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-final quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c-- TH-V
c-- Id,aqq=(epinv*(epinv-L)+1/2*L^2+3/2*epinv+[pi]^2/6)*[delta(1-x)]
c--  +(1-x-2/[1-x]*[ln(2-x)]
c--  -(1+x)*(L+[ln(1-x)])-(1+x^2)*[ln(x)]/[1-x]
c--  +4*[ln(1-x)/(1-xp)]+2*L/[1-xp]
      
      if (vorz == 1) then
        ifn_qq=epinv*(epinv2-L)+half*L**2+pisqo6
        if (scheme == 'tH-V') then
        return
        elseif (scheme == 'dred') then
        ifn_qq=ifn_qq-half
        return
        endif
      endif
      
      omx=one-x
      lomx=log(omx)
      
      if (vorz == 2) then
        Pqqreg=-one-x
        Pqqpr=omx
        ltmx=log(two-x)
        lx=log(x)
        ifn_qq=Pqqpr+Pqqreg*(lomx+log(alfa)+L-epinv-lx)
     &    -two/omx*(lx+log((omx+alfa)/alfa))
        return
      endif
      
      ifn_qq=two/omx*(two*lomx+L-epinv)
      
      return
      end

***************************** Gluon-Gluon *****************************
      function ifn_gg(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: ifn_gg
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx,ltmx,Pggreg
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-final gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c-- TH-V
c-- Id,agg=[delta(1-x)]*(
c--  epinv*(epinv-L)+1/2*L^2+11/6*epinv+[pi]^2/6-1/3*epinv*nf/xn)
c--  +2*(-1+(1-x)/x+x*(1-x))*(L-[ln(x)]+[ln(1-x)])
c--  -2*[ln(2-x)]/[1-x]-2*[ln(x)]/[1-x]
c--  +4*[ln(1-x)/(1-xp)]+2*L/[1-xp]
      
      if (vorz == 1) then
        ifn_gg=epinv*(epinv2-L)+half*L**2+pisqo6
c     &       +(11._dp-two*nf/xn)*epinv/6._dp
        if (scheme == 'tH-V') then
        return
        elseif (scheme == 'dred') then 
        ifn_gg=ifn_gg-1._dp/6._dp
        return
        endif
      endif
      
      omx=one-x
      lomx=log(omx)
      
      if (vorz == 2) then
        Pggreg=two*(omx/x+x*omx-one)
        ltmx=log(two-x)
        lx=log(x)
        ifn_gg=Pggreg*(lomx+log(alfa)+L-epinv-lx)
     &    -two/omx*(lx+log((omx+alfa)/alfa))
        return
      endif
      
       ifn_gg=two/omx*(two*lomx+L-epinv)
      
      return
      end

C----Not necessary because non-soft singular initial states
C----are included as initial-initial dipoles
***************************** Quark-Gluon *****************************
c      function ifn_qg(x,L,vorz)
c      implicit none
      include 'types.f'
      real(dp):: ifn_qg
c      
c      integer:: vorz
c      real(dp):: x,L,omx,lx,lomx
c      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
c      include 'epinv.f'
c--- returns the integral of the subtraction term for an
c--- initial-final gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c-- TH-V
c-- Id,aqg=(1-2*x*(1-x))*(L-[ln(x)]+[ln(1-x)])+2*x*(1-x)
c      
c      ifn_qg=0._dp
c      if ((vorz == 1).or.(vorz == 3)) return
c      
c      omx=one-x
c      lomx=log(omx)
c      lx=log(x)
c      
c      if (vorz == 2) then
c        ifn_qg=(one-two*x*omx)*(lomx-lx+L-epinv)+two*x*omx
c      endif
      
c      return
c      end
      
***************************** Gluon-Quark *****************************
c      function ifn_gq(x,L,vorz)
c      implicit none
      include 'types.f'
      real(dp):: ifn_gq
c      
c      integer:: vorz
c      real(dp):: x,L,omx,lx,lomx
c      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
c      include 'epinv.f'
c--- returns the integral of the subtraction term for an
c--- initial-final gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c-- TH-V
c-- Id,agq=(1+(1-x)^2)/x*(L-[ln(x)]+[ln(1-x)])+x
c      
c      ifn_gq=0._dp
c      if ((vorz == 1).or.(vorz == 3)) return
      
c      omx=one-x
c      lomx=log(omx)
c      lx=log(x)
      
c      if (vorz == 2) then
c        ifn_gq=(one+omx**2)/x*(lomx-lx+L-epinv)+x
c      endif
      
c      return
c      end

***********************************************************************
**************************** FINAL-INITIAL ****************************
***********************************************************************

***************************** Quark-Quark *****************************
      function fin_qq(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: fin_qq
      
      integer:: vorz
      real(dp):: x,L,omx,gamq
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- final-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--MSbar
c-- Id,aqq=(epinv*(epinv-L)+1/2*L^2+3/2*(epinv-L)+7/2-[pi]^2/2)*[delta(1-x)]
c--  +2/[1-x]*[ln(2-x)]
c--  +(-2*[ln(1-x)/(1-xp)]-3/2/[1-xp])
      
      gamq=1.5_dp
      if (vorz == 1) then
         fin_qq=epinv*(epinv2-L)+half*L**2+gamq*(epinv-L)
     &   +3.5_dp-half*pisq

         if (scheme == 'tH-V') then
           return
         elseif (scheme == 'dred') then
           fin_qq=fin_qq-half
         return
         endif
      endif
      
      omx=one-x
      
      if (vorz == 2) then
        fin_qq=two*log(two-x)/omx
        if (alfa < omx) fin_qq=fin_qq-(2._dp*log((2._dp-x)/omx)+gamq)/omx
        return
      endif
      
      fin_qq=-(two*log(omx)+gamq)/omx
      
      return
      end

***************************** Gluon-Gluon *****************************
      function fin_gg(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: fin_gg
      
      integer:: vorz
      real(dp):: x,L,omx,gamg
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--MSbar
c-- Id,aqg=(-2/3*(epinv-L)-10/9)*[delta(1-x)]
c--  +0
c--  +2/3/[1-xp]
c-- Id,agg=
c--  (2*epinv*(epinv-L)+L^2+(epinv-L)*11/3+67/9-[pi]^2)
c--  *[delta(1-x)]
c--  +4*[ln(2-x)]/[1-x]
c--  +2*(-2*[ln(1-x)/(1-xp)]-11/6/[1-xp])

      gamg=11._dp/3._dp-2._dp/3._dp*real(nf,dp)/xn
      if (vorz == 1) then
        fin_gg=two*epinv*(epinv2-L)+L**2+gamg*(epinv-L)
     &         +67._dp/9._dp-real(nf,dp)/xn*10._dp/9._dp-pisq
        if (scheme == 'tH-V') then
        return
        elseif (scheme == 'dred') then
        fin_gg=fin_gg-real(nf,dp)/xn/3._dp
        return
        endif
      endif
      
      omx=one-x
      
      if (vorz == 2) then
        fin_gg=four*log(two-x)/omx
        if (omx > alfa) 
     &  fin_gg=fin_gg-(four*log((2._dp-x)/omx)+gamg)/omx
        return
      endif
      
      fin_gg=-(four*log(omx)+gamg)/omx
      return
      end




***************************** Quark-Gluon *****************************
c      function fin_qg(x,L,vorz)
c      implicit none
      include 'types.f'
      real(dp):: fin_qg
c      
c      integer:: vorz
c      real(dp):: x,L,omx
c      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
c      include 'epinv.f'
c      include 'epinv2.f'
c      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--MSbar
c--Id,aqg=(-2/3*(epinv-L)-10/9)*[delta(1-x)]
c-- +0
c-- +2/3/[1-xp]

      
c      if (vorz == 1) then
c       fin_qg=2._dp/3._dp*(-epinv+L)-10._dp/9._dp
c       if (scheme == 'tH-V') then
c          return
c       elseif (scheme == 'dred') then
c          fin_qg=fin_qg-1._dp/3._dp
c          return
c       endif
c      elseif (vorz == 2) then
c        fin_qg=0._dp
c      elseif (vorz == 3) then
c        fin_qg=2._dp/3._dp/(one-x)
c      endif
c      return
c      end


***********************************************************************
***************************** FINAL-FINAL *****************************
***********************************************************************

***************************** Quark-Quark *****************************
      function ffn_qq(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: ffn_qq
      
      integer:: vorz
      real(dp):: x,L
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- final-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C --MSbar
c Id,aqq=epinv*(epinv-L)+1/2*L^2+3/2*(epinv-L)+5-[pi]^2/2
      
      ffn_qq=0._dp
      if (vorz == 1) then
        ffn_qq=epinv*(epinv2-L)+half*L**2+1.5_dp*(epinv-L)+5._dp-half*pisq
        if (alfa<1._dp) 
     &  ffn_qq=ffn_qq-log(alfa)**2+1.5_dp*(alfa-1._dp-log(alfa))
        if (scheme == 'tH-V') then
        return
        elseif (scheme == 'dred') then
        ffn_qq=ffn_qq-half
        return
        endif
      endif
      return
      end

***************************** Quark-Gluon *****************************
c      function ffn_qg(x,L,vorz)
c      implicit none
      include 'types.f'
      real(dp):: ffn_qg
c      
c      integer:: vorz
c      real(dp):: x,L
c      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
c      include 'epinv.f'
c      include 'epinv2.f'
c      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)
C --MSbar
c Id,aqg=-2/3*(epinv-L)-16/9
c      
c      ffn_qg=0._dp
c      if (vorz == 1) then
c        ffn_qg=-2._dp/3._dp*(epinv-L)-16._dp/9._dp
c        if (scheme == 'tH-V') then
c          return
c        elseif (scheme == 'dred') then
c          ffn_qg=ffn_qg-1._dp/3._dp
c          return
c       endif
c      endif
c      return
c      end

***************************** Gluon-Gluon *****************************
      function ffn_gg(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: ffn_gg
      
      integer:: vorz
      real(dp):: x,L
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C --MSbar
c Id,aqg=-2/3*(epinv-L)-16/9
c Id,agg=2*epinv*(epinv-L)+L^2+11/3*(epinv-L)+100/9-[pi]^2
      
      ffn_gg=0._dp
      if (vorz == 1) then
        ffn_gg=two*epinv*(epinv2-L)+L**2+100._dp/9._dp-pisq
     &        +11._dp/3._dp*(epinv-L)
        ffn_gg=ffn_gg+real(nf,dp)/xn*(-2._dp/3._dp*(epinv-L)-16._dp/9._dp)
        if (alfa<1._dp) 
     &  ffn_gg=ffn_gg-log(alfa)**2
     &  +(11._dp/3._dp-2._dp/3._dp*real(nf,dp)/xn)*(alfa-1._dp-log(alfa))
        if (scheme == 'tH-V') then
          return
        elseif (scheme == 'dred') then
          ffn_gg=ffn_gg-real(nf,dp)/xn/3._dp
          return
        endif
      endif
      return
      end

