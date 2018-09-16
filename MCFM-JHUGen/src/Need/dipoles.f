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
      function ii_qq(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: ii_qq
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx
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

      
cIIqq = 
c  + 1 - x - 2*[ln(x)]*[1-x]^-1 + [ln(x)]*[1+x] + [ln(al(x))]*
c [(1+x^2)/(1-x)] - [1+x]*L - 2*[1+x]*[ln(1-x)] + [1+x]*epinv
c 
c + [delta(1-x)]
c  * ( 1/2*L^2 - 1/6*pisq - epinv*L + epinv^2 )
c 
c + [1/(1-x)_(0)]
c  * ( 2*L + 4*[ln(1-x)] - 2*epinv )
 
      if (vorz == 1) then
        ii_qq=epinv*(epinv2-L)+0.5_dp*L**2-pisqo6
        if (scheme == 'tH-V') then
          return
        elseif (scheme == 'dred') then
          ii_qq=ii_qq-half
          return
      else
        write(6,*) 'Value of scheme not implemented properly ',scheme
        stop
        endif
      endif
      
      omx=one-x
      lomx=log(omx)
      lx=log(x)
      
      if (vorz == 2) then
        ii_qq=omx-(one+x)*(two*lomx+L-epinv)-(one+x**2)/omx*lx
        if (omx > aii) ii_qq=ii_qq+(one+x**2)/omx*log(aii/omx)
        return
      endif
      
      ii_qq=two/omx*(two*lomx+L-epinv)
      
      return
      end

***************************** Quark-Gluon *****************************
      function ii_qg(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: ii_qg
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx
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

c IIqg = 
c  + 2*x - 2*x^2 - [ln(x)]*[x^2+(1-x)^2] + [ln(al(x))]*
c [x^2+(1-x)^2] + [x^2+(1-x)^2]*L + 2*[x^2+(1-x)^2]*[ln(1-x)]
c  - [x^2+(1-x)^2]*epinv
 

      ii_qg=0._dp
      if ((vorz == 1) .or. (vorz == 3)) return
      
      omx=one-x
      lomx=log(omx)
      lx=log(x)
      
      if (vorz == 2) then
        ii_qg=(one-two*x*omx)*(two*lomx-lx+L-epinv)+two*x*omx
        if (omx > aii) ii_qg=ii_qg+(one-two*x*omx)*log(aii/omx)
      endif
      return
      end
      
***************************** Gluon-Quark *****************************
      function ii_gq(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: ii_gq
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx
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

      
c IIgq =  + x - [ln(x)]*[(1+(1-x)^2)/x] + [ln(al(x))]*
c  [(1+(1-x)^2)/x] + [(1+(1-x)^2)/x]*L + 2*[(1+(1-x)^2)/x]*
c  [ln(1-x)] - [(1+(1-x)^2)/x]*epinv
 
      ii_gq=0._dp
      if ((vorz == 1) .or. (vorz == 3)) return
      
      omx=one-x
      lomx=log(omx)
      lx=log(x)
      
      if (vorz == 2) then
        ii_gq=(one+omx**2)/x*(two*lomx-lx+L-epinv)+x
        if (omx > aii) ii_gq=ii_gq+(one+omx**2)/x*log(aii/omx)
        return
      endif

      return
      end

***************************** Gluon-Gluon *****************************
      function ii_gg(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: ii_gg
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx
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
c--  -nflav/3/xn*epinv)*[delta(1-x)]
c--  -2*[ln(x)]/[1-x]
c--  +2*(-1+x*(1-x)+(1-x)/x)*(-[ln(x)]+L+2*[ln(1-x)])
c--  +(4*[ln(1-x)/(1-xp)]+2*L/[1-xp])
      
c  IIgg = 
c   - 2*[ln(x)]*[1-x]^-1 - 2*[ln(x)]*[(1-x)/x-1+x*(1-x)] + 2*
c   [1-x]^-1*[ln(al(x))] + 2*[ln(al(x))]*[(1-x)/x-1+x*(1-x)]
c   + 2*[(1-x)/x-1+x*(1-x)]*L + 4*[(1-x)/x-1+x*(1-x)]*[ln(1-x)]
c   - 2*[(1-x)/x-1+x*(1-x)]*epinv
c 
c   + [delta(1-x)]
c    * ( 1/2*L^2 - 1/6*pisq - epinv*L + epinv^2 )
c 
c   + [1/(1-x)_(0)]
c    * ( 2*L + 4*[ln(1-x)] - 2*epinv )
 
      if (vorz == 1) then
        ii_gg=epinv*(epinv2-L)+half*L**2-pisqo6
        if (scheme == 'tH-V') then
          return
        elseif (scheme == 'dred') then
          ii_gg=ii_gg-1._dp/6._dp
          return
      else
        write(6,*) 'Value of scheme not implemented properly ',scheme
        stop
        endif
      endif
      
      omx=one-x
      lomx=log(omx)
      
      if (vorz == 2) then
        lx=log(x)
        ii_gg=two*(omx/x+x*omx-one)*(two*lomx-lx+L-epinv)-two*lx/omx
        if (omx > aii) ii_gg=ii_gg
     &   +two*(one/omx+omx/x+x*omx-one)*log(aii/omx)
        return
      endif
      
      ii_gg=two*(two*lomx+L-epinv)/omx
      
      return
      end

***********************************************************************
**************************** INITIAL-FINAL ****************************
***********************************************************************

***************************** Quark-Quark *****************************
 
 
      function if_qq(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: if_qq
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx,ltmx
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
        if_qq=epinv*(epinv2-L)+half*L**2+pisqo6
        if (scheme == 'tH-V') then
          return
        elseif (scheme == 'dred') then
          if_qq=if_qq-half
          return
      else
        write(6,*) 'Value of scheme not implemented properly ',scheme
        stop
        endif
      endif
      
      omx=one-x
      lomx=log(omx)
      
cIFqq = 
c  + 1 - x - 2*[ln(x)]*[1-x]^-1 + [ln(x)]*[1+x] - 2*[1-x]^-1*
c [ln((1+al-x)/al)] - [ln(al)]*[1+x] - [1+x]*L - [1+x]*[ln(1-x)]
c  + [1+x]*epinv
c 
c + [delta(1-x)]
c  * ( 1/2*L^2 + 1/6*pisq - epinv*L + epinv^2 )
c 
c + [1/(1-x)_(0)]
c  * ( 2*L + 4*[ln(1-x)] - 2*epinv )

      if (vorz == 2) then
        ltmx=log((omx+aif)/aif)
        lx=log(x)
        if_qq=omx-two/omx*ltmx-(one+x)*(lomx+L-epinv)-(one+x**2)/omx*lx
        if_qq=if_qq-log(aif)*(one+x)
        return
      endif
      
      if_qq=two/omx*(two*lomx+L-epinv)
      
      return
      end

***************************** Gluon-Gluon *****************************
      function if_gg(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: if_gg
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx,ltmx
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
c--  epinv*(epinv-L)+1/2*L^2+11/6*epinv+[pi]^2/6-1/3*epinv*nflav/xn)
c--  +2*(-1+(1-x)/x+x*(1-x))*(L-[ln(x)]+[ln(1-x)])
c--  -2*[ln(2-x)]/[1-x]-2*[ln(x)]/[1-x]
c--  +4*[ln(1-x)/(1-xp)]+2*L/[1-xp]
      
      if (vorz == 1) then
        if_gg=epinv*(epinv2-L)+half*L**2+pisq/6._dp
        if (scheme == 'tH-V') then
          return
        elseif (scheme == 'dred') then 
          if_gg=if_gg-1._dp/6._dp
          return
      else
        write(6,*) 'Value of scheme not implemented properly ',scheme
        stop
        endif
      endif
      
      omx=one-x
      lomx=log(omx)
      
cIFgg = 
c  - 2*[ln(x)]*[1-x]^-1 - 2*[ln(x)]*[(1-x)/x-1+x*(1-x)] - 2*
c [1-x]^-1*[ln((1+al-x)/al)] + 2*[ln(al)]*[(1-x)/x-1+x*(1-x)]
c  + 2*[(1-x)/x-1+x*(1-x)]*L + 2*[(1-x)/x-1+x*(1-x)]*[ln(1-x)]
c  - 2*[(1-x)/x-1+x*(1-x)]*epinv
c 
c + [delta(1-x)]
c  * ( 1/2*L^2 + 1/6*pisq - epinv*L + epinv^2 )
c 
c + [1/(1-x)_(0)]
c  * ( 2*L + 4*[ln(1-x)] - 2*epinv )
 


      if (vorz == 2) then
        ltmx=log((omx+aif)/aif)
        lx=log(x)
        if_gg=two*((lomx-lx+L-epinv+log(aif))*(omx/x+x*omx-one)
     &  -(ltmx+lx)/omx)
        return
      endif
      
      if_gg=two/omx*(two*lomx+L-epinv)
      
      return
      end

***************************** Quark-Gluon *****************************
      function if_qg(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: if_qg
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-final gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c-- TH-V
c-- Id,aqg=(1-2*x*(1-x))*(L-[ln(x)]+[ln(1-x)])+2*x*(1-x)
cIFqg = 
c  + 2*x - 2*x^2 - [ln(x)]*[x^2+(1-x)^2] + [ln(al)]*
c [x^2+(1-x)^2] + [x^2+(1-x)^2]*L + [x^2+(1-x)^2]*[ln(1-x)]
c  - [x^2+(1-x)^2]*epinv

      
      if_qg=0._dp
      if ((vorz == 1).or.(vorz == 3)) return
      
      omx=one-x
      lomx=log(omx)
      lx=log(x)
      
      if (vorz == 2) then
        if_qg=(one-two*x*omx)*(lomx-lx+L-epinv+log(aif))+two*x*omx
      endif
     
      return
      end
      
***************************** Gluon-Quark *****************************
      function if_gq(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: if_gq
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-final gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c-- TH-V
c-- Id,agq=(1+(1-x)^2)/x*(L-[ln(x)]+[ln(1-x)])+x
cIFgq =  + x - [ln(x)]*[(1+(1-x)^2)/x] + [ln(al)]*
c [(1+(1-x)^2)/x] + [(1+(1-x)^2)/x]*L + [(1+(1-x)^2)/x]*
c [ln(1-x)] - [(1+(1-x)^2)/x]*epinv
 
      
      if_gq=0._dp
      if ((vorz == 1).or.(vorz == 3)) return
      
      omx=one-x
      lomx=log(omx)
      lx=log(x)
      
      if (vorz == 2) then
        if_gq=(one+omx**2)/x*(lomx-lx+L-epinv+log(aif))+x
      endif
      
      return
      end

***********************************************************************
**************************** FINAL-INITIAL ****************************
***********************************************************************

***************************** Quark-Quark *****************************
      function fi_qq(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: fi_qq
      
      integer:: vorz
      real(dp):: x,L,omx,theta
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
      
cFIq = 
c  + 2*[1-x]^-1*[ln(2-x)]*[Theta(x-1+al)]
c 
c + 2*[ln(1/(1-x))/(1-x)_(1-al)]
c 
c + [delta(1-x)]
c  * ( 7/2 - 3/2*L + 1/2*L^2 - 1/2*pisq - 3/2*[ln(al)] - 
c [ln(al)]^2 + 3/2*epinv - epinv*L + epinv^2 )
c 
c - 3/2*[1/(1-x)_(1-al)]
      theta=0._dp 
      if (x > 1._dp-afi) theta=1._dp       
      if (vorz == 1) then
         fi_qq=epinv*(epinv2-L)+half*L**2+1.5_dp*(epinv-L)
     &   +3.5_dp-half*pisq-log(afi)*(1.5_dp+log(afi))
         if (scheme == 'tH-V') then
           return
         elseif (scheme == 'dred') then
           fi_qq=fi_qq-half
           return
       else
         write(6,*) 'Value of scheme not implemented properly ',scheme
         stop
         endif
      endif
      
      omx=one-x
      
      if (vorz == 2) then
        fi_qq=two*log(two-x)/omx*theta
        return
      endif
      
      fi_qq=-(two*log(omx)+1.5_dp)/omx*theta
      
      return
      end

***************************** Gluon-Gluon *****************************
      function fi_gg(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: fi_gg
      
      integer:: vorz
      real(dp):: x,L,omx,theta
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
      include 'nflav.f'
      include 'b0.f'
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

cFIg = 
c  + 4*[1-x]^-1*[ln(2-x)]*[Theta(x-1+al)]
c 
c + 4*[ln(1/(1-x))/(1-x)_(1-al)]
c 
c + [delta(1-x)]
c  * ( 67/9 + L^2 - pisq - 2*[ln(al)]^2 - 2*epinv*L + 2*epinv^2 )
c 
c + [delta(1-x)]*CA^-1
c  * ( - 20/9*Tr*nflav - 2*[ln(al)]*b0 - 2*b0*L + 2*b0*epinv )
c 
c + [1/(1-x)_(1-al)]*CA^-1
c  * ( - 2*b0 )

      theta=0._dp 
      if (x > 1._dp-afi) theta=1._dp       
      if (vorz == 1) then
        fi_gg=two*epinv*(epinv2-L)+L**2
     &  +67._dp/9._dp-10._dp/9._dp*real(nflav,dp)/xn
     &  -pisq+2._dp*b0/xn*(epinv-L)-2._dp*log(afi)*(b0/xn+log(afi))
        if (scheme == 'tH-V') then
          return
        elseif (scheme == 'dred') then
          fi_gg=fi_gg-1._dp/3._dp
          return
      else
        write(6,*) 'Value of scheme not implemented properly ',scheme
        stop
        endif
      endif
      
      omx=one-x
      
      if (vorz == 2) then
        fi_gg=four*log(two-x)/omx*theta
        return
      endif
      
      fi_gg=-(four*log(omx)
     &       +11._dp/3._dp-real(nflav,dp)/xn*2._dp/3._dp)/omx*theta
      return
      end




***************************** Quark-Gluon *****************************
c      function fi_qg(x,L,vorz)
c      implicit none
c      include 'types.f'
c      real(dp):: fi_qg
c      
c      integer:: vorz
c      real(dp):: x,L,omx
c      include 'constants.f'
c      include 'nf.f'
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
c       fi_qg=2._dp/3._dp*(-epinv+L)-10._dp/9._dp
c       if (scheme == 'tH-V') then
c          return
c       elseif (scheme == 'dred') then
c          fi_qg=fi_qg-1._dp/3._dp
c          return
c       endif
c      elseif (vorz == 2) then
c        fi_qg=0._dp
c      elseif (vorz == 3) then
c        fi_qg=2._dp/3._dp/(one-x)
c      endif
c      return
c      end


***********************************************************************
***************************** FINAL-FINAL *****************************
***********************************************************************

***************************** Quark-Quark *****************************
      function ff_qq(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: ff_qq
      
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
      
cFFq = 
c  + 7/2 - 3/2*L + 1/2*L^2 + 3/2*al - 1/2*pisq - 3/2*[ln(al)]
c  - [ln(al)]^2 + 3/2*epinv - epinv*L + epinv^2

      ff_qq=0._dp
      if (vorz == 1) then
        ff_qq=epinv*(epinv2-L)+half*L**2+1.5_dp*(epinv-L)+5._dp-half*pisq
        ff_qq=ff_qq+1.5_dp*(aff-1._dp-log(aff))-log(aff)**2
        if (scheme == 'tH-V') then
          return
        elseif (scheme == 'dred') then
          ff_qq=ff_qq-half
          return
      else
        write(6,*) 'Value of scheme not implemented properly ',scheme
        stop
        endif
      endif
      return
      end

***************************** Quark-Gluon *****************************
c      function ff_qg(x,L,vorz)
c      implicit none
c      include 'types.f'
c      real(dp):: ff_qg
c      
c      integer:: vorz
c      real(dp):: x,L
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'epinv.f'
c      include 'epinv2.f'
c      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)
C --MSbar
c Id,aqg=-2/3*(epinv-L)-16/9
c      
c      ff_qg=0._dp
c      if (vorz == 1) then
c        ff_qg=-2._dp/3._dp*(epinv-L)-16._dp/9._dp
c        if (scheme == 'tH-V') then
c          return
c        elseif (scheme == 'dred') then
c          ff_qg=ff_qg-1._dp/3._dp
c          return
c       endif
c      endif
c      return
c      end

***************************** Gluon-Gluon *****************************
      function ff_gg(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: ff_gg
      
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
      include 'nflav.f'
      include 'colstruc.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c--- 26/11/09: modified to enable separation of CA and TR pieces
c---           (used in checks of single top + b process)

C --MSbar
c Id,aqg=-2/3*(epinv-L)-16/9
c Id,agg=2*epinv*(epinv-L)+L^2+11/3*(epinv-L)+100/9-[pi]^2
      
cFFg = 
c  + 67/9 + L^2 - pisq - 2*[ln(al)]^2 - 2*epinv*L + 2*epinv^2
c + CA^-1
c  * ( 2*al*b0 - 20/9*Tr*nflav - 2*[ln(al)]*b0 - 2*b0*L + 2*b0*epinv ) + 0.

      ff_gg=0._dp
      if (vorz == 1) then
        if (nfonly) then
          ff_gg=0._dp
        else
          ff_gg=two*epinv*(epinv2-L)+L**2+100._dp/9._dp-pisq
     &          +two*11._dp/6._dp*(epinv-L)
          ff_gg=ff_gg-two*log(aff)**2+two*11._dp/6._dp*(aff-1._dp-log(aff))
          if (scheme == 'tH-V') then
            continue ! the above is the CT in this scheme
          elseif (scheme == 'dred') then
            ff_gg=ff_gg-1._dp/3._dp
                     ! no return yet, need to include nflav piece
          else
            write(6,*)'Value of scheme not implemented properly ',scheme
            stop
          endif
        endif
        if (caonly) then
          continue ! nothing more to do
        else
          ff_gg=ff_gg-4._dp/3._dp*tr*real(nflav,dp)/ca*(epinv-L)
     &          -real(nflav,dp)/ca*16._dp/9._dp
          ff_gg=ff_gg-4._dp/3._dp*tr*real(nflav,dp)/ca*(aff-1._dp-log(aff))
        endif
      endif
      
      return
      end


*********** Gluon -> Quark-Antiquark final state splitting ***********
*********** (just the nf pieces of ff_gg)
      function ff_gq(x,L,vorz)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nflav.f'
      include 'epinv.f'
      include 'alfacut.f'
      real(dp):: ff_gq
      integer:: vorz
      real(dp):: x,L
      
      ff_gq=0._dp
      if (vorz == 1) then
        ff_gq=ff_gq-4._dp/3._dp*tr*real(nflav,dp)/ca*(epinv-L)
     &        -real(nflav,dp)/ca*16._dp/9._dp
        ff_gq=ff_gq-4._dp/3._dp*tr*real(nflav,dp)/ca*(aff-1._dp-log(aff))
      endif
      
      return
      end

*********** Gluon -> Quark-Antiquark final state splitting ***********
*********** (just the nf pieces of fi_gg)
      function fi_gq(x,L,vorz)
      implicit none
      include 'types.f'
      real(dp):: fi_gq
      integer:: vorz
      real(dp):: x,L,omx,theta
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'alfacut.f'
      include 'nflav.f'
      include 'b0.f'

      theta=0._dp 
      if (x > 1._dp-afi) theta=1._dp
      if (vorz == 1) then
        fi_gq=
     &  -2._dp*real(nflav,dp)/3._dp/xn*(epinv-L)
     &  -10._dp/9._dp*real(nflav,dp)/xn
     &  +2._dp*log(afi)*real(nflav,dp)/3._dp/xn
        return
      endif
      
      omx=one-x
      
      if (vorz == 2) then
        fi_gq=zip
        return
      endif
      
      fi_gq=+real(nflav,dp)/xn*2._dp/3._dp/omx*theta
      
      return
      end


