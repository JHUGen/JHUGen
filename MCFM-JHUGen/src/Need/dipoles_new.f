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
      double precision function iin_qq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,Pqqreg,alfax
      include 'constants.f'
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

      
      if (vorz .eq. 1) then
        iin_qq=epinv*(epinv2-L)+0.5d0*L**2+1.5d0*epinv-pisqo6
     .   -epinv*1.5d0
        if (scheme .eq. 'tH-V') then
           return
        elseif (scheme .eq. 'dred') then
           iin_qq=iin_qq-half
           return
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        Pqqreg=-one-x
        iin_qq=omx+Pqqreg*(two*lomx+L-epinv)-(one+x**2)/omx*lx
        alfax=alfa/omx
        if (alfax .lt. 1d0) iin_qq=iin_qq+(two/omx+Pqqreg)*log(alfax)
        return
      endif
      
      iin_qq=two/omx*(two*lomx+L-epinv)
      
      return
      end

***************************** Quark-Gluon *****************************
      double precision function iin_qg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,Pqgreg,alfax
      include 'constants.f'
      include 'epinv.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--TH-V
c-- Id,aqg=(1-2*x*(1-x))*(-[ln(x)]+L+2*[ln(1-x)])+2*x*(1-x)
      iin_qg=0d0
      if ((vorz .eq. 1) .or. (vorz .eq. 3)) return
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        Pqgreg=one-two*x*omx
        iin_qg=Pqgreg*(two*lomx-lx+L-epinv)+two*x*omx
        alfax=alfa/omx
        if (alfax .lt. 1d0) iin_qg=iin_qg+Pqgreg*log(alfax)
      endif
      return
      end
      
***************************** Gluon-Quark *****************************
      double precision function iin_gq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,Pgqreg,alfax
      include 'constants.f'
      include 'epinv.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial quark-quark (--> gluon) antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--TH-V
c-- Id,agq=(1+(1-x)^2)/x*(-[ln(x)]+L+2*[ln(1-x)])+x

      
      iin_gq=0d0
      if ((vorz .eq. 1) .or. (vorz .eq. 3)) return
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        Pgqreg=(one+omx**2)/x
        iin_gq=Pgqreg*(two*lomx-lx+L-epinv)+x
        alfax=alfa/omx
        if (alfax .lt. 1d0) iin_gq=iin_gq+Pgqreg*log(alfax)
        return
      endif

      return
      end

***************************** Gluon-Gluon *****************************
      double precision function iin_gg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,Pggreg,alfax
      include 'constants.f'
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
      
      if (vorz .eq. 1) then
        iin_gg=epinv*(epinv2-L)+half*L**2-pisqo6
     .       +(11d0-two*nf/xn)/6d0*epinv
     .       -(11d0-two*nf/xn)/6d0*epinv
        if (scheme .eq. 'tH-V') then
        return
        elseif (scheme .eq. 'dred') then
        iin_gg=iin_gg-1d0/6d0
        return
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      
      if (vorz .eq. 2) then
        Pggreg=omx/x+x*omx-one
        lx=dlog(x)
        iin_gg=two*Pggreg*(two*lomx-lx+L-epinv)-two*lx/omx
        alfax=alfa/omx
        if (alfax .lt. 1d0) iin_gg=iin_gg+Pggreg*log(alfax)
        return
      endif
      
      iin_gg=two*(two*lomx+L-epinv)/omx
      
      return
      end

***********************************************************************
**************************** INITIAL-FINAL ****************************
***********************************************************************

***************************** Quark-Quark *****************************
      double precision function ifn_qq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,ltmx,Pqqpr,Pqqreg
      include 'constants.f'
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
      
      if (vorz .eq. 1) then
        ifn_qq=epinv*(epinv2-L)+half*L**2+pisqo6
        if (scheme .eq. 'tH-V') then
        return
        elseif (scheme .eq. 'dred') then
        ifn_qq=ifn_qq-half
        return
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      
      if (vorz .eq. 2) then
        Pqqreg=-one-x
        Pqqpr=omx
        ltmx=dlog(two-x)
        lx=dlog(x)
        ifn_qq=Pqqpr+Pqqreg*(lomx+log(alfa)+L-epinv-lx)
     .    -two/omx*(lx+log((omx+alfa)/alfa))
        return
      endif
      
      ifn_qq=two/omx*(two*lomx+L-epinv)
      
      return
      end

***************************** Gluon-Gluon *****************************
      double precision function ifn_gg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,ltmx,Pggreg
      include 'constants.f'
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
      
      if (vorz .eq. 1) then
        ifn_gg=epinv*(epinv2-L)+half*L**2+pisqo6
c     .       +(11d0-two*nf/xn)*epinv/6d0
        if (scheme .eq. 'tH-V') then
        return
        elseif (scheme .eq. 'dred') then 
        ifn_gg=ifn_gg-1d0/6d0
        return
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      
      if (vorz .eq. 2) then
        Pggreg=two*(omx/x+x*omx-one)
        ltmx=dlog(two-x)
        lx=dlog(x)
        ifn_gg=Pggreg*(lomx+log(alfa)+L-epinv-lx)
     .    -two/omx*(lx+log((omx+alfa)/alfa))
        return
      endif
      
       ifn_gg=two/omx*(two*lomx+L-epinv)
      
      return
      end

C----Not necessary because non-soft singular initial states
C----are included as initial-initial dipoles
***************************** Quark-Gluon *****************************
c      double precision function ifn_qg(x,L,vorz)
c      implicit none
c      integer vorz
c      double precision x,L,omx,lx,lomx
c      include 'constants.f'
c      include 'epinv.f'
c--- returns the integral of the subtraction term for an
c--- initial-final gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c-- TH-V
c-- Id,aqg=(1-2*x*(1-x))*(L-[ln(x)]+[ln(1-x)])+2*x*(1-x)
c      
c      ifn_qg=0d0
c      if ((vorz .eq. 1).or.(vorz .eq. 3)) return
c      
c      omx=one-x
c      lomx=dlog(omx)
c      lx=dlog(x)
c      
c      if (vorz .eq. 2) then
c        ifn_qg=(one-two*x*omx)*(lomx-lx+L-epinv)+two*x*omx
c      endif
      
c      return
c      end
      
***************************** Gluon-Quark *****************************
c      double precision function ifn_gq(x,L,vorz)
c      implicit none
c      integer vorz
c      double precision x,L,omx,lx,lomx
c      include 'constants.f'
c      include 'epinv.f'
c--- returns the integral of the subtraction term for an
c--- initial-final gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c-- TH-V
c-- Id,agq=(1+(1-x)^2)/x*(L-[ln(x)]+[ln(1-x)])+x
c      
c      ifn_gq=0d0
c      if ((vorz .eq. 1).or.(vorz .eq. 3)) return
      
c      omx=one-x
c      lomx=dlog(omx)
c      lx=dlog(x)
      
c      if (vorz .eq. 2) then
c        ifn_gq=(one+omx**2)/x*(lomx-lx+L-epinv)+x
c      endif
      
c      return
c      end

***********************************************************************
**************************** FINAL-INITIAL ****************************
***********************************************************************

***************************** Quark-Quark *****************************
      double precision function fin_qq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,gamq
      include 'constants.f'
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
      
      gamq=1.5d0
      if (vorz .eq. 1) then
         fin_qq=epinv*(epinv2-L)+half*L**2+gamq*(epinv-L)
     .   +3.5d0-half*pisq

         if (scheme .eq. 'tH-V') then
           return
         elseif (scheme .eq. 'dred') then
           fin_qq=fin_qq-half
         return
         endif
      endif
      
      omx=one-x
      
      if (vorz .eq. 2) then
        fin_qq=two*dlog(two-x)/omx
        if (alfa .lt. omx) fin_qq=fin_qq-(2d0*log((2d0-x)/omx)+gamq)/omx
        return
      endif
      
      fin_qq=-(two*dlog(omx)+gamq)/omx
      
      return
      end

***************************** Gluon-Gluon *****************************
      double precision function fin_gg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,gamg
      include 'constants.f'
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

      gamg=11d0/3d0-2d0/3d0*dfloat(nf)/xn
      if (vorz .eq. 1) then
        fin_gg=two*epinv*(epinv2-L)+L**2+gamg*(epinv-L)
     .         +67d0/9d0-dfloat(nf)/xn*10d0/9d0-pisq
        if (scheme .eq. 'tH-V') then
        return
        elseif (scheme .eq. 'dred') then
        fin_gg=fin_gg-dfloat(nf)/xn/3d0
        return
        endif
      endif
      
      omx=one-x
      
      if (vorz .eq. 2) then
        fin_gg=four*dlog(two-x)/omx
        if (omx .gt. alfa) 
     .  fin_gg=fin_gg-(four*dlog((2d0-x)/omx)+gamg)/omx
        return
      endif
      
      fin_gg=-(four*dlog(omx)+gamg)/omx
      return
      end




***************************** Quark-Gluon *****************************
c      double precision function fin_qg(x,L,vorz)
c      implicit none
c      integer vorz
c      double precision x,L,omx
c      include 'constants.f'
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

      
c      if (vorz .eq. 1) then
c       fin_qg=2d0/3d0*(-epinv+L)-10d0/9d0
c       if (scheme .eq. 'tH-V') then
c          return
c       elseif (scheme .eq. 'dred') then
c          fin_qg=fin_qg-1d0/3d0
c          return
c       endif
c      elseif (vorz .eq. 2) then
c        fin_qg=0d0
c      elseif (vorz .eq. 3) then
c        fin_qg=2d0/3d0/(one-x)
c      endif
c      return
c      end


***********************************************************************
***************************** FINAL-FINAL *****************************
***********************************************************************

***************************** Quark-Quark *****************************
      double precision function ffn_qq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- final-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C --MSbar
c Id,aqq=epinv*(epinv-L)+1/2*L^2+3/2*(epinv-L)+5-[pi]^2/2
      
      ffn_qq=0d0
      if (vorz .eq. 1) then
        ffn_qq=epinv*(epinv2-L)+half*L**2+1.5d0*(epinv-L)+5d0-half*pisq
        if (alfa.lt.1d0) 
     .  ffn_qq=ffn_qq-log(alfa)**2+1.5d0*(alfa-1d0-log(alfa))
        if (scheme .eq. 'tH-V') then
        return
        elseif (scheme .eq. 'dred') then
        ffn_qq=ffn_qq-half
        return
        endif
      endif
      return
      end

***************************** Quark-Gluon *****************************
c      double precision function ffn_qg(x,L,vorz)
c      implicit none
c      integer vorz
c      double precision x,L
c      include 'constants.f'
c      include 'epinv.f'
c      include 'epinv2.f'
c      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)
C --MSbar
c Id,aqg=-2/3*(epinv-L)-16/9
c      
c      ffn_qg=0d0
c      if (vorz .eq. 1) then
c        ffn_qg=-2d0/3d0*(epinv-L)-16d0/9d0
c        if (scheme .eq. 'tH-V') then
c          return
c        elseif (scheme .eq. 'dred') then
c          ffn_qg=ffn_qg-1d0/3d0
c          return
c       endif
c      endif
c      return
c      end

***************************** Gluon-Gluon *****************************
      double precision function ffn_gg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L
      include 'constants.f'
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
      
      ffn_gg=0d0
      if (vorz .eq. 1) then
        ffn_gg=two*epinv*(epinv2-L)+L**2+100d0/9d0-pisq
     .        +11d0/3d0*(epinv-L)
        ffn_gg=ffn_gg+dfloat(nf)/xn*(-2d0/3d0*(epinv-L)-16d0/9d0)
        if (alfa.lt.1d0) 
     .  ffn_gg=ffn_gg-log(alfa)**2
     .  +(11d0/3d0-2d0/3d0*dfloat(nf)/xn)*(alfa-1d0-log(alfa))
        if (scheme .eq. 'tH-V') then
          return
        elseif (scheme .eq. 'dred') then
          ffn_gg=ffn_gg-dfloat(nf)/xn/3d0
          return
        endif
      endif
      return
      end

