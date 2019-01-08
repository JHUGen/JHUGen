***************************** Gluon-Gluon *****************************
      double precision function ii_gg_fac(x,L,Lfac,vorz)
      implicit none
      integer vorz
      double precision x,L,Lfac,omx,lx,lomx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
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
        ii_gg_fac=epinv*(epinv2-L)+half*L**2-pisqo6
     .       +(11/6d0-dfloat(nf)/3d0/xn)*(epinv-L)
     .       -(11/6d0-dfloat(nf)/3d0/xn)*(epinv-Lfac)
        if (scheme .eq. 'tH-V') then
           return
        elseif (scheme .eq. 'dred') then
           ii_gg_fac=ii_gg_fac-1d0/6d0
           return
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      
      if (vorz .eq. 2) then
        lx=dlog(x)
        ii_gg_fac=two*(omx/x+x*omx-one)*(two*lomx-lx+Lfac-epinv)
     .           -two*lx/omx
        return
      endif
      
      ii_gg_fac=two*(two*lomx+Lfac-epinv)/omx
      
      return
      end


***************************** Gluon-Gluon *****************************
      double precision function if_gg_fac(x,L,Lfac,vorz)
      implicit none
      integer vorz
      double precision x,L,Lfac,omx,lx,lomx,ltmx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
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
        if_gg_fac=epinv*(epinv2-L)+half*L**2+pisq/6d0
     .       +(11/6d0-dfloat(nf)/3d0/xn)*(epinv-L)
     .       -(11/6d0-dfloat(nf)/3d0/xn)*(epinv-Lfac)
        if (scheme .eq. 'tH-V') then
        return
        elseif (scheme .eq. 'dred') then 
        if_gg_fac=if_gg_fac-1d0/6d0
        return
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      
      if (vorz .eq. 2) then
        ltmx=dlog(two-x)
        lx=dlog(x)
        if_gg_fac=two*((lomx-lx+Lfac-epinv)*(omx/x+x*omx-one)
     .           -(ltmx+lx)/omx)
        return
      endif
      
      if_gg_fac=two/omx*(two*lomx+Lfac-epinv)
      
      return
      end

***************************** Quark-Quark *****************************
      double precision function ii_qq_fac(x,L,Lfac,vorz)
      implicit none
      integer vorz
      double precision x,L,Lfac,omx,lx,lomx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--TH-V
c-- Id,aqq=
c--  [delta(1-x)]*(epinv*(epinv-L)+1/2*L^2+3/2*epinv-[pi]^2/6)
c--  +(1-x)-(1+x)*(L+2*[ln(1-x)])-(1+x^2)*[ln(x)]/[1-x]
c--  +4*[ln(1-x)/(1-xp)]+2*L/[1-xp]

      
      if (vorz .eq. 1) then
        ii_qq_fac=epinv*(epinv2-L)+0.5d0*L**2-pisqo6
     .   +1.5d0*(epinv-L)
     .   -1.5d0*(epinv-Lfac)
        if (scheme .eq. 'tH-V') then
           return
        elseif (scheme .eq. 'dred') then
           ii_qq_fac=ii_qq_fac-half
           return
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        ii_qq_fac=omx-(one+x)*(two*lomx+Lfac-epinv)-(one+x**2)/omx*lx
        return
      endif
      
      ii_qq_fac=two/omx*(two*lomx+Lfac-epinv)
      
      return
      end

***************************** Quark-Quark *****************************
      double precision function if_qq_fac(x,L,Lfac,vorz)
      implicit none
      integer vorz
      double precision x,L,Lfac,omx,lx,lomx,ltmx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- initial-final quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)
c-- TH-V
c-- Id,aqq=(epinv*(epinv-L)+1/2*L^2+3/2*epinv+[pi]^2/6)*[delta(1-x)]
c--  +(1-x-2/[1-x]*[ln(2-x)]
c--  -(1+x)*(L+[ln(1-x)])-(1+x^2)*[ln(x)]/[1-x]
c--  +4*[ln(1-x)/(1-xp)]+2*L/[1-xp]
      
      if (vorz .eq. 1) then
        if_qq_fac=epinv*(epinv2-L)+half*L**2+pisqo6
     .   +1.5d0*(epinv-L)
     .   -1.5d0*(epinv-Lfac)
        if (scheme .eq. 'tH-V') then
          return
        elseif (scheme .eq. 'dred') then
          if_qq_fac=if_qq_fac-half
          return
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      
      if (vorz .eq. 2) then
        ltmx=dlog(two-x)
        lx=dlog(x)
        if_qq_fac=omx-two*ltmx/omx-(one+x)*(lomx+Lfac-epinv)
     .           -(one+x**2)/omx*lx
        return
      endif
      
      if_qq_fac=two/omx*(two*lomx+Lfac-epinv)
      
      return
      end

