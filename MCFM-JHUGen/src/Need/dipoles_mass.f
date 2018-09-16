************************************************************************
*     Author: J. M. Campbell                                           *
*     April 15th, 2004                                                 *
*                                                                      *
*     Routines which return various pieces of the integrated           *
*     MASSIVE subtraction terms, used in both _v and _z routines       *
*                                                                      *
*     The formulae implemented here are derived in the FORM programs   *
*     testif.frm, testif_gg.frm and testfi.frm                         *
*                                                                      *
*     Other final-initial and all final-final dipoles are untested     *
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
*************************** Initial-INITIAL ***************************
***********************************************************************

***************************** Quark-Quark *****************************
      function ii_mqq(x,L,mbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: ii_mqq
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx,mbar,Pqqreg,alfax
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
      
      if (vorz == 1) then
        ii_mqq=epinv*(epinv2-L)+0.5_dp*L**2-pisqo6
        if (scheme == 'tH-V') then
           return
        elseif (scheme == 'dred') then
           ii_mqq=ii_mqq-half
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
        Pqqreg=-one-x
        ii_mqq=omx+Pqqreg*(two*lomx+L-epinv)-(one+x**2)/omx*lx
        alfax=aii/omx
        if (alfax < 1._dp) ii_mqq=ii_mqq+(two/omx+Pqqreg)*log(alfax)
        return
      endif
      
      ii_mqq=two/omx*(two*lomx+L-epinv)
      
      return
      end

***************************** Quark-Gluon *****************************
      function ii_mqg(x,L,mbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: ii_mqg
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx,Pqgreg,alfax,mbar
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial gluon-quark antenna, either

      ii_mqg=0._dp
      if ((vorz == 1) .or. (vorz == 3)) return
      
      omx=one-x
      lomx=log(omx)
      lx=log(x)
      
      if (vorz == 2) then
        Pqgreg=one-two*x*omx
        ii_mqg=Pqgreg*(two*lomx-lx+L-epinv)+two*x*omx
        alfax=aii/omx
        if (alfax < 1._dp) ii_mqg=ii_mqg+Pqgreg*log(alfax)
      endif
      return
      end
      
***************************** Gluon-Quark *****************************
      function ii_mgq(x,L,mbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: ii_mgq
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx,Pgqreg,alfax,mbar
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial quark-quark (--> gluon) antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)  
      
      ii_mgq=0._dp
      if ((vorz == 1) .or. (vorz == 3)) return
      
      omx=one-x
      lomx=log(omx)
      lx=log(x)
      
      if (vorz == 2) then
        Pgqreg=(one+omx**2)/x
        ii_mgq=Pgqreg*(two*lomx-lx+L-epinv)+x
        alfax=aii/omx
        if (alfax < 1._dp) ii_mgq=ii_mgq+Pgqreg*log(alfax)
        return
      endif

      return
      end

***************************** Gluon-Gluon *****************************
      function ii_mgg(x,L,mbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: ii_mgg
      
      integer:: vorz
      real(dp):: x,L,omx,lx,lomx,Pggreg,alfax,mbar
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
      
      if (vorz == 1) then
        ii_mgg=epinv*(epinv2-L)+half*L**2-pisqo6
        if (scheme == 'tH-V') then
          return
        elseif (scheme == 'dred') then
          ii_mgg=ii_mgg-1._dp/6._dp
          return
      else
        write(6,*) 'Value of scheme not implemented properly ',scheme
        stop
        endif
      endif
      
      omx=one-x
      lomx=log(omx)
      
      if (vorz == 2) then
        Pggreg=omx/x+x*omx-one
        lx=log(x)
        ii_mgg=two*Pggreg*(two*lomx-lx+L-epinv)-two*lx/omx
        alfax=aii/omx
        if (alfax < 1._dp) ii_mgg=ii_mgg
     &   +two*(one/omx+Pggreg)*log(alfax)
        return
      endif
      
      ii_mgg=two*(two*lomx+L-epinv)/omx
      
      return
      end

***********************************************************************
**************************** INITIAL-FINAL ****************************
***********************************************************************

***************************** Quark-Quark *****************************
      function if_mqq(x,L,mbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: if_mqq
      
      integer:: vorz
      real(dp):: x,L,mbar,mbarsq,omx,lx,lomx,Pqqreg,ddilog,zp
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'alfacut.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- initial-final quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     

      mbarsq=mbar**2
      if_mqq=0._dp
      if (vorz == 1) then
        if_mqq=(epinv+log(1._dp+mbarsq))*(epinv-L)+half*L**2
     &   -half*log(1._dp+mbarsq)**2+2._dp*log(mbarsq)*log(1._dp+mbarsq)
     &   +2._dp*ddilog(-mbarsq)+pisqo6
c--- correct form for double pole (normally zero)
        if_mqq=if_mqq-epinv**2+epinv*epinv2     
        if (scheme == 'tH-V') then
          return
        elseif (scheme == 'dred') then
          if_mqq=if_mqq-half
          return
      else
        write(6,*) 'Value of scheme not implemented properly ',scheme
        stop
        endif
      endif
      omx=one-x
      lomx=log(omx)
      zp=omx/(omx+x*mbarsq)
      if (vorz == 2) then
         Pqqreg=-(1._dp+x)
         lx=log(x)
         if_mqq=Pqqreg*(-(epinv-L)+2._dp*lomx-lx-log(x*mbarsq+omx))
     &    +omx-2._dp/omx*(lx+log((1._dp+x*mbarsq+omx)/(1._dp+mbarsq)))
         if (aif < zp) then
         if_mqq=if_mqq-(two/omx*(log(zp*(omx+aif)/(aif*(omx+zp))))
     &    +Pqqreg*log(zp/aif))
         endif
      elseif (vorz == 3) then
         if_mqq=2._dp/omx*(-(epinv-L)+2._dp*lomx-log(1._dp+mbarsq))
      endif
      return
      end

***************************** Gluon-Gluon *****************************
      function if_mgg(x,L,mbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: if_mgg
      
      integer:: vorz
      real(dp):: x,L,omx,lx,mbar,mbarsq,Pggreg,ddilog,zp,omzp
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'alfacut.f'
      include 'scheme.f'

      mbarsq=mbar**2
      if_mgg=0._dp

      if (vorz == 1) then
          if_mgg=(epinv+log(1._dp+mbarsq))*(epinv-L)+half*L**2
     &     -half*log(1._dp+mbarsq)**2+2._dp*log(mbarsq)*log(1._dp+mbarsq)
     &     +2._dp*ddilog(-mbarsq)+pisqo6
c--- correct form for double pole (normally zero)
        if_mgg=if_mgg-epinv**2+epinv*epinv2     
        if (scheme == 'tH-V') then
          return
        elseif (scheme == 'dred') then
          if_mgg=if_mgg-1._dp/6._dp
          return
       else
        write(6,*) 'Value of scheme not implemented properly ',scheme
        stop
        endif
        return
      endif
      omx=one-x
      zp=omx/(omx+x*mbarsq)
      if (vorz == 2) then
         Pggreg=2._dp*(omx/x-1._dp+x*omx)
         lx=log(x)
         if_mgg=Pggreg*(-(epinv-L)+2._dp*log(omx)-lx-log(x*mbarsq+omx))
     &    +2._dp*mbarsq*log(x*mbarsq/(x*mbarsq+omx))
     &    -2._dp/omx*(lx+log((1._dp+x*mbarsq+omx)/(1._dp+mbarsq)))
         if (aif < zp) then
           if (aif == 1._dp) then
             write(6,*) 'zp > 1 in dipoles_mass.f - this is forbidden'
             stop
           endif  
         omzp=x*mbarsq/(omx+x*mbarsq)
           if_mgg=if_mgg-(two/omx*(log(zp*(omx+aif)/(aif*(omx+zp))))
     &     +Pggreg*log(zp/aif)+2._dp*mbarsq*log(omzp/(1._dp-aif)))
         endif
         return
      elseif (vorz == 3) then
         if_mgg=2._dp/omx*(-(epinv-L)+2._dp*log(omx)-log(1._dp+mbarsq))
         return 
      endif
      return
      end

***************************** Quark-Gluon *****************************
c--- Not necessary because for off-diagonal (no soft singularity)
c--- we always choose to use the initial spectator
      
***************************** Gluon-Quark *****************************
c--- Not necessary because for off-diagonal (no soft singularity)
c--- we always choose to use the initial spectator


***********************************************************************
**************************** FINAL-INITIAL ****************************
***********************************************************************

***************************** Quark-Quark *****************************
      function fi_mqq(x,L,mbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: fi_mqq
      
      integer:: vorz
      real(dp):: x,L,mbar,mbarsq,omx,ddilog
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- final-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     

      mbarsq=mbar**2
      if (vorz == 1) then
        fi_mqq=(1._dp+log(mbarsq/(1._dp+mbarsq)))*(epinv-L)
     &        +log(mbarsq)+half*log(mbarsq)**2
     &        +half*log(1._dp+mbarsq)**2
     &        -2._dp*log(mbarsq)*log(1._dp+mbarsq)
     &        -2._dp*ddilog(-mbarsq)
     &        +2._dp-pisq/3._dp
     &        +2._dp*log(afi)*(log((1._dp+mbarsq)/mbarsq)-1._dp)
        return
      endif
      
      omx=one-x
      
      if (vorz == 2) then
        if (x > 1._dp-afi) then
          fi_mqq=+omx/2._dp/(x*mbarsq+omx)**2
     &     +2._dp/omx*(log((1._dp+x*mbarsq+omx)*mbarsq/
     &                    ((1._dp+mbarsq)*(omx+x*mbarsq))))
        else
          fi_mqq=0._dp
        endif        
        return
      endif
      
      if (vorz == 3) then
        if (x > 1._dp-afi) then
          fi_mqq=2._dp/omx*(log((1._dp+mbarsq)/(mbarsq))-1._dp)
        else
          fi_mqq=0._dp
        endif        
      endif
      
      return
      end


***********************************************************************
***************************** FINAL-FINAL *****************************
***********************************************************************

***************************** Quark-Quark *****************************
*           (case when both emitter and spectator are massive)        *
***********************************************************************
      function ff_mqq(x,L,mbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: ff_mqq
      
      integer:: vorz
      real(dp):: x,L,mbar,mbarsq,Lro,ro,vtijk,arg,ddilog,ff_alpha

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
C     mbarsq=mass**2/Qsq
C     L=log(Qsq/musq)
c--- returns the integral of the subtraction term for an
c--- final-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C --MSbar
c Id,aqq=epinv*(epinv-L)+1/2*L^2+3/2*(epinv-L)+5-[pi]^2/2

c g zipffqq=colfac*del(omx)*(
c  1/vtijk*((epinv-L)*Lro-1/2*Lro**2-2*Lro*log(1-4*mbarsq)
c           +4*dilog(-ro)-6*dilog(1-ro)+pisq/3)
c +epinv-L+3-2*log(1-2*mbar)+log(1-mbar)+1/2*log(mbarsq)
c  -2*mbar/(1-2*mbarsq)*(1-2*mbar)-mbar/(1-mbar)
c  -2*mbarsq/(1-2*mbarsq)*log(mbar/(1-mu))
c  )-ffqq;

c      if (aff < 1._dp) then
c        write(6,*) 'Integrated dipole routine ff_mqq does not yet'
c      write(6,*) ' support values of alpha < 1.'
c      stop
c      endif

      ff_mqq=0._dp
      if (vorz == 1) then
        mbarsq=mbar**2
        arg=1._dp-4._dp*mbarsq
        if (arg < 0._dp) then
            write(6,*) 'Threshold problem in ff_mqq'
            stop
        endif
        vtijk=sqrt(arg)/(1._dp-2._dp*mbarsq)
        ro=sqrt((1._dp-vtijk)/(1._dp+vtijk))
        Lro=log(ro)
        ff_mqq=1._dp/vtijk*((epinv-L)*Lro-half*Lro**2-2._dp*Lro*log(arg)
     &   +4._dp*ddilog(-ro)-6._dp*ddilog(1._dp-ro)+pisq/3._dp)
     &   +epinv-L+3._dp
     &   -2._dp*log(1._dp-2._dp*mbar)+log(1._dp-mbar)+half*log(mbarsq)
     &   -2._dp*mbar/(1._dp-2._dp*mbarsq)*(1._dp-2._dp*mbar)-mbar/(1._dp-mbar)
     &   -2._dp*mbarsq/(1._dp-2._dp*mbarsq)*log(mbar/(1._dp-mbar))
c        Ieik=1._dp/vtijk*(half*(epinv-L)*Lro-Lro*log(arg)
c     &  -log(roj)**2+pisq/6._dp
c     &  +2._dp*ddilog(-ro)-2._dp*ddilog(1._dp-ro)-ddilog(1._dp-roj**2))
c        Icoll=(epinv-L)+half*Lmbarsq-2._dp-2._dp*log((1._dp-mbar)**2-mbarsq)
c     & +log(1._dp-mbar)-2._dp*mbarsq/(1._dp-2._dp*mbarsq)*log(mbar/(1._dp-mbar))
c     & +5._dp-mbar/(1._dp-mbar)-2._dp*mbar*(1._dp-2._dp*mbar)/(1._dp-2._dp*mbarsq)
c        ff_mqq=2._dp*Ieik+Icoll

c--- now add on the contribution for aff < 1
c        if (abs(aff-1._dp) > 1.e-8_dp) then
      ff_mqq=ff_mqq+ff_alpha(mbar,mbar)
c      endif 
        return
      endif
      return
      end


***************************** Quark-Quark *****************************
*   (case when emitter and spectator have different non-zero masses)  *
***********************************************************************
      function ff_2mqq(x,L,mbar,sbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: ff_2mqq
      
      integer:: vorz
      real(dp):: x,L,mbar,mbarsq,sbar,sbarsq,
     & Lro,ro,roj,rok,vtijk,arg,arg2,den,ddilog

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'alfacut.f'

      if (aff < 1._dp) then
        write(6,*) 'Integrated dipole routine ff_2mqq does not yet'
      write(6,*) ' support values of alpha < 1.'
      stop
      endif

      ff_2mqq=0._dp
      if (vorz == 1) then
        mbarsq=mbar**2
      sbarsq=sbar**2
        arg=1._dp+mbarsq**2+sbarsq**2-2._dp*(mbarsq+sbarsq+mbarsq*sbarsq)
      arg2=1._dp-2._dp*sbar+sbarsq-mbarsq
      den=1._dp-mbarsq-sbarsq
        if (arg < 0._dp) then
            write(6,*) 'Threshold problem in ff_2mqq'
            stop
        endif
        vtijk=sqrt(arg)/den
c--- some redundancy here: ro=roj*rok
      roj=sqrt((1._dp-vtijk+2._dp*mbarsq/den)/(1._dp+vtijk+2._dp*mbarsq/den))
      rok=sqrt((1._dp-vtijk+2._dp*sbarsq/den)/(1._dp+vtijk+2._dp*sbarsq/den))
        ro=sqrt((1._dp-vtijk)/(1._dp+vtijk))
        Lro=log(ro)
      ff_2mqq=1._dp/vtijk*((epinv-L)*Lro-log(roj)**2-log(rok)**2
     &   -2._dp*Lro*log(1._dp-(mbar+sbar)**2)
     &   +4._dp*ddilog(-ro)-4._dp*ddilog(1._dp-ro)
     &   -ddilog(1._dp-roj**2)-ddilog(1._dp-rok**2)+pisq/3._dp)
     &   +epinv-L+3._dp
     &   +log(1._dp-sbar)+half*log(mbarsq)-2._dp*log(arg2)
     &   +(4._dp*sbarsq-2._dp*sbar-2._dp*mbarsq*log(mbar/(1._dp-sbar)))/den
     &   -sbar/(1._dp-sbar)
        return
      endif
      
      return
      end


***************************** Quark-Quark *****************************
*    (case when the emitter is massless and spectator is massive,     *
*    plus the case when emitter is massive and spectator is massless) *
***********************************************************************
      function ff_1mqq(x,L,mbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: ff_1mqq
      
      integer:: vorz
      real(dp):: x,L,mbar,mbarsq,ddilog,afftmp
     & ,Icolla,Ieika,Icollb,Ieikb,arg,ommsq,logm,logomm,xp,yp
     & ,arg1,arg2,arg3,ypp,ypm


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
      include 'phi.f'

      ff_1mqq=0._dp 
      phi=1._dp


      if (vorz == 1) then
        if (mbar > 1._dp) then
            write(6,*) 'Problem with mbar in ff_1mqq, mbar=',mbar
            stop
        endif


      mbarsq=mbar**2
      ommsq=1._dp-mbarsq
      logm=log(mbarsq)
      logomm=log(ommsq)

c--- Note: this routine uses the exact form of the eikonal integrals
c---       given in the appendix, not the symmetric combination

C----radiation from massive line with massless spectator
      afftmp=aff       
      arg=afftmp+(1._dp-afftmp)*mbarsq
      Ieika=
     & half*logm*(epinv-L)-2._dp*ddilog(ommsq)
     & -logm*logomm-0.25_dp*logm**2
     & -log(afftmp)*logm-ddilog(-ommsq/mbarsq)
     & +ddilog(-afftmp*ommsq/mbarsq)
      Icolla=
     &  epinv-L+phi+2._dp+(1._dp+half*phi)*logm
     & +(phi-2._dp)*logm/(ommsq)-2._dp*logomm
     & +half*phi*(3._dp*afftmp-2._dp-(3._dp-mbarsq)/ommsq*log(arg)
     & -afftmp/arg)
     & -2._dp*log(afftmp)+2._dp*log(arg)/ommsq


C----radiation from massless line with massive spectator
      yp=(1._dp-mbar)/(1._dp+mbar)
      afftmp=aff*yp       
      xp=(yp-afftmp)+sqrt((yp-afftmp)*(1._dp/yp-afftmp))
      arg1=0.25_dp*(1._dp-yp**2+2*xp*yp)
      arg2=half*(1._dp+YP-xp)
      arg3=half*(1._dp-YP+xp)
      ypp=half*(1._dp+yp)
      ypm=half*(1._dp-yp)

      Ieikb=0._dp
     & +half*epinv*epinv2-half*epinv*L+0.25_dp*L**2
     & -logomm*(epinv-L)
     & +ddilog(ommsq)-2.5_dp*pisqo6+logomm**2
     & +half*log(arg1/(arg2*arg3))**2-log(arg2/ypp)**2
     & +2._dp*(log(ypp)*log(arg3/ypm)+log(ypp/yp)*log(arg1/(ypp*ypm))
     & +DDILOG(ypm/ypp)-DDILOG(arg1/ypp**2)
     & +DDILOG(arg3)-DDILOG(ypm))

      Icollb=1.5_dp*(epinv-L)-3._dp*log(1._dp-mbar)+5._dp-mbar/(1._dp-mbar)
     & -2._dp*mbar*(1._dp-2._dp*mbar)/ommsq
     & +1.5_dp*(log(YP/afftmp)-YP+afftmp)

c--- Note: extra factor of half because we include this term once for each
c---  leg, but this is the sum of both legs
      ff_1mqq=half*(2._dp*(Ieika+Ieikb)+Icolla+Icollb)
    

        if (scheme == 'tH-V') then
           return
        elseif (scheme == 'dred') then
c--- Note: extra factor of half because we include this term once for each
c---  leg, but this is the sum of both legs (as above)
           ff_1mqq=ff_1mqq-half*half
           return
        endif
      endif
      
      return
      end
      
***************************** Quark-Quark *****************************
*     (case when the emitter is massive and spectator is massless)    *
***********************************************************************
      function ff_mqq0(x,L,mbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: ff_mqq0
      
      integer:: vorz
      real(dp):: x,L,mbar,mbarsq,ddilog
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for a
c--- final-final quark-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C --MSbar
c g zipffqq0=colfac*del(omx)*(
c epinv^2/2+epinv*(1-L/2+log(mbarsq)/2-log(1-mbarsq))-L+3+L**2/4-pisq*5/12
c  -dilog(1-mbarsq)-log(mbarsq)*(log(mbarsq)/2-1+L)/2
c  -log(1-mbarsq)*(log(mbarsq,1-mbarsq)-L+2)
c  -mbarsq/(1-mbarsq)*log(mbarsq)
c  )-ffqq0;

c--- Note: this uses the symmetric form of the eikonal integrals, so
c---       that it can match with ff_mgg

      if (aff < 1._dp) then
        write(6,*) 'Integrated dipole routine ff_mqq0 does not yet'
      write(6,*) ' support values of alpha < 1.'
      stop
      endif

      ff_mqq0=0._dp 
      if (vorz == 1) then
        if (mbar > 1._dp) then
            write(6,*) 'Problem with mbar in ff_mqq0, mbar=',mbar
            stop
        endif

        mbarsq=mbar**2
        ff_mqq0=epinv*epinv2/2._dp
     &         +epinv*(1._dp-L/2._dp+log(mbarsq)/2._dp-log(1._dp-mbarsq))
     &         -L+3._dp+L**2/4._dp-pisq*5._dp/12._dp
     &         -ddilog(1._dp-mbarsq)
     &         -log(mbarsq)*(log(mbarsq)/2._dp-1._dp+L)/2._dp
     &         -log(1._dp-mbarsq)*(log(mbarsq/(1._dp-mbarsq))-L+2._dp)
     &         -mbarsq/(1._dp-mbarsq)*log(mbarsq)
      
      endif
      
      return
      end

***************************** Gluon-Gluon *****************************
*                (in this case the spectator is massive)              *
***********************************************************************
      function ff_mgg(x,L,mbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: ff_mgg
      
      integer:: vorz
      real(dp):: x,L,mbar,mbarsq,ddilog

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'scheme.f'
      include 'nflav.f'
      include 'alfacut.f'
      include 'colstruc.f'
c--- returns the integral of the subtraction term for a
c--- final-final gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C --MSbar
c g zipffgg=colfac*del(omx)*(
c  (epinv-L)*(epinv+log(mbarsq)-2*log(1-mbarsq))+L**2/2+100/9-pisq*5/6
c  +2*b0/xn*(epinv-L)
c  -2*dilog(1-mbarsq)-log(mbarsq)^2/2
c  -2*log(1-mbarsq)*(log(mbarsq,1-mbarsq))
c  -22/3*(log(1-mbar)+mbar/(1+mbar))
c  +4/3*mbarsq/(1-mbarsq)*log(2*mbar/(1+mu))
c  -nfoN*4/3*(4/3-log(1-mbar)-mbar/(1+mbar))
c  )-ffgg;

      if (aff < 1._dp) then
        write(6,*) 'Integrated dipole routine ff_mgg does not yet'
      write(6,*) ' support values of alpha < 1.'
      stop
      endif

      ff_mgg=0._dp
      if (vorz == 1) then
        mbarsq=mbar**2
c      ff_mgg=(epinv-L)*(epinv+log(mbarsq)-2._dp*log(1._dp-mbarsq))
c     &   +L**2/2._dp+100._dp/9._dp-pisq*5._dp/6._dp
c     &   +2._dp*b0/xn*(epinv-L)
c     &   -2._dp*ddilog(1._dp-mbarsq)-log(mbarsq)**2/2._dp
c     &   -2._dp*log(1._dp-mbarsq)*log(mbarsq/(1._dp-mbarsq))
c     &   -22._dp/3._dp*(log(1._dp-mbar)+mbar/(1._dp+mbar))
c     &   +4._dp/3._dp*mbarsq/(1-mbarsq)*log(2._dp*mbar/(1._dp+mbar))
c     &   -real(nflav,dp)/xn*4._dp/3._dp*(
c     &       4._dp/3._dp-log(1._dp-mbar)-mbar/(1._dp+mbar)
c     &      +mbarsq/(1._dp-mbarsq)*log(2._dp*mbar/(1._dp+mbar)))
c--- NB: last term added by JC, 6/10/08 
        if (nfonly) then
        ff_mgg=0._dp
      else
         ff_mgg=(epinv-L)*(epinv+log(mbarsq)-2._dp*log(1._dp-mbarsq))
     &    +L**2/2._dp+100._dp/9._dp-pisq*5._dp/6._dp
     &    +2._dp*11._dp/6._dp*(epinv-L)
     &    -2._dp*ddilog(1._dp-mbarsq)-log(mbarsq)**2/2._dp
     &    -2._dp*log(1._dp-mbarsq)*log(mbarsq/(1._dp-mbarsq))
     &    -22._dp/3._dp*(log(1._dp-mbar)+mbar/(1._dp+mbar))
     &    +4._dp/3._dp*mbarsq/(1._dp-mbarsq)*log(2._dp*mbar/(1._dp+mbar))
          if (scheme == 'tH-V') then
            continue ! the above is the CT in this scheme
          elseif (scheme == 'dred') then
            ff_mgg=ff_mgg-1._dp/3._dp
            return
        else
          write(6,*)'Value of scheme not implemented properly ',scheme
          stop
          endif
        endif
        if (caonly) then
          continue ! nothing more to do
      else
        ff_mgg=ff_mgg
     &   -tr*4._dp/3._dp*real(nflav,dp)/ca*(epinv-L)
     &   -real(nflav,dp)/ca*2._dp*tr*4._dp/3._dp*(
     &       4._dp/3._dp-log(1._dp-mbar)-mbar/(1._dp+mbar)
     &      +mbarsq/(1._dp-mbarsq)*log(2._dp*mbar/(1._dp+mbar)))
         endif
      endif
      
      return
      end

      
      
      
      
      
      

************************************************************************
************************************************************************
*     BELOW HERE functionS HAVE NOT BEEN CHECKED (AND ARE NOT USED)    *
************************************************************************
************************************************************************


***************************** Gluon-Gluon *****************************
      function fi_mgg(x,L,mbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: fi_mgg
      
      integer:: vorz
      real(dp):: x,L,mbar,omx
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
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


      if (vorz == 1) then
        fi_mgg=two*epinv*(epinv2-L)+L**2+67._dp/9._dp-pisq
     &    +11._dp*(epinv-L)/3._dp
      endif
      
      omx=one-x
      
      if (vorz == 2) then
        fi_mgg=four*log(two-x)/omx
        return
      endif
      
      fi_mgg=-(four*log(omx)+11._dp/3._dp)/omx
      return
      end




***************************** Quark-Gluon *****************************
      function fi_mqg(x,L,mbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: fi_mqg
      
      integer:: vorz
      real(dp):: x,L,mbar,mbarsq,omx,rt,JaS,JaNS
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--MSbar

      mbarsq=mbar**2
      if (1._dp-4._dp*mbarsq < 0._dp) then
      write(6,*) 'error in fi_mqg,(1._dp-4._dp*mbarsq < 0._dp)'
      stop
      else
      rt=sqrt(1._dp-4._dp*mbarsq)
      endif
      if (vorz == 1) then
CDTS 5.63
      JaS=-2._dp/3._dp*log(mbarsq)-10._dp/9._dp
CDTS 5.64
      JaNS=10._dp/9._dp*(1._dp-rt)-8._dp/9._dp*mbarsq*rt
     & +4._dp/3._dp*log(half*(1._dp+rt))
      fi_mqg=JaS+JaNS
       
      elseif (vorz == 2) then
C regular
        fi_mqg=0._dp
      elseif (vorz == 3) then
C plus at x=x+
      omx=1._dp-x
CDTS 5.62
      fi_mqg=2._dp/3._dp*(omx+2._dp*mbarsq)/omx**2*sqrt(1._dp-4._dp*mbarsq/omx)
      endif
      return
      end









***************************** Quark-Gluon *****************************
      function ff_mqg(x,L,mbar,vorz)
      implicit none
      include 'types.f'
      real(dp):: ff_mqg
      
      integer:: vorz
      real(dp):: x,L,mbar,mbarsq,ro,arg
C-----L=Log(Qsq/musq)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
     
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)
C --MSbar
c Id,aqg=-2/3*(epinv-L)-16/9
c      
      ff_mqg=0._dp
      if (vorz == 1) then
        mbarsq=mbar**2
        arg=1._dp-4._dp*mbarsq
        if (arg < 0._dp) then
        write(6,*) 'Threshold problem in ff_mqg'
        stop
        else
        ro=sqrt(arg)  
        ff_mqg=-2._dp/3._dp*(2._dp*log(mbarsq)
     &   -2._dp*log(half*(1._dp+ro))+2._dp/3._dp*ro*(3._dp+ro**2))
        return
        endif
      endif
      return
      end

***************************** Gluon-Gluon *****************************
c      function ff_mgg(x,L,mbar,vorz)
c      implicit none
c      include 'types.f'
c      real(dp):: ff_mgg
c      
c      integer:: vorz
c      real(dp):: x,L,mbar,Ieik,Icoll
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'epinv.f'
c      include 'epinv2.f'
cc--- returns the integral of the subtraction term for an
cc--- final-initial gluon-gluon antenna, either
cc--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
cC --MSbar
cc Id,aqg=-2/3*(epinv-L)-16/9
cc Id,agg=2*epinv*(epinv-L)+L^2+11/3*(epinv-L)+100/9-[pi]^2
c      
c      ff_mgg=0._dp
cCDST Eq.5.32
c      if (vorz == 1) then
c      Ieik=half*epinv*(epinv2-L)+half*L**2-pisq/4._dp
c      Icoll=11._dp/6._dp*(epinv-L)+50._dp/9._dp
cCDST Eq.5.36 needs to be added too - need to fix this
c      Icoll=Icoll+half*real(nf,dp)/xn*(-2._dp/3._dp*(epinv-L)-16._dp/9._dp)
c      ff_mgg=2._dp*(2._dp*Ieik+Icoll)
c      return
c      endif
c      return
c      end

c
c * Now write these expressions in a neater form
c g zipffqq=colfac*del(omx)*(
c  1/vtijk*((epinv-L)*Lro-1/2*Lro**2-2*Lro*log(1-4*mbarsq)
c           +4*dilog(-ro)-6*dilog(1-ro)+pisq/3)
c +epinv-L+3-2*log(1-2*mbar)+log(1-mbar)+1/2*log(mbarsq)
c  -2*mbar/(1-2*mbarsq)*(1-2*mbar)-mbar/(1-mbar)
c  -2*mbarsq/(1-2*mbarsq)*log(mbar/(1-mu))
c  )-ffqq;

c g zipifqq=colfac*(
c  del(omx)*((epinv+log(1+mbarsq))*(epinv-L)+1/2*L**2
c            +1/2*log(1+mbarsq)**2+2*dilog(1/(1+musq))-pisq/6)
c +Preg(q,q)/colfac*(-(epinv-L)+2*log(omx)-log(x)-log(x*mbarsq+omx))
c  +omx-2/omx*(log(x)+log(1+x*mbarsq+omx,1+mbarsq))
c +2/omxp*(-[epinv-L]+2*log(omx)-log(1+mbarsq))
c  )-ifqq;

c g zipifgg=colfac*(
c  del(omx)*((epinv+log(1+mbarsq))*(epinv-L)+1/2*L**2
c            +1/2*log(1+mbarsq)**2+2*dilog(1/(1+musq))-pisq/6)
c +Preg(g,g)/colfac*(-(epinv-L)+2*log(omx)-log(x)-log(x*mbarsq+omx))
c  +2*mbarsq*log(x*mbarsq,x*mbarsq+omx)
c  -2/omx*(log(x)+log(1+x*mbarsq+omx,1+mbarsq))
c +2/omxp*(-[epinv-L]+2*log(omx)-log(1+mbarsq))
c  )-ifgg;

c g zipfiqq=colfac*(
c  del(omx)*((1+log(mbarsq,1+mbarsq))*(epinv-L)
c            +1/2*log(mbarsq)-1/2*log(mbarsq)**2
c            +1/2*log(1+mbarsq)+1/2*log(1+mbarsq)**2
c            -2*log(mbarsq)*log(1+mbarsq)
c            -4*dilog(-mbarsq)+mbarsq/2/(1+mbarsq)
c            +3/2-2/3*pisq)
c +2/omx*(log(1+x*mbarsq+omx,1+mbarsq))
c +omx/2/(x*mbarsq+omx)**2
c +2/omxp*(log(1+mbarsq,omx+x*mbarsq)-1)
c  )-fiqq;

