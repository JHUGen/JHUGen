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
      double precision function ii_mqq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,mbar,Pqqreg,alfax
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
      
      if (vorz .eq. 1) then
        ii_mqq=epinv*(epinv2-L)+0.5d0*L**2-pisqo6
        if (scheme .eq. 'tH-V') then
           return
        elseif (scheme .eq. 'dred') then
           ii_mqq=ii_mqq-half
           return
      else
        write(6,*) 'Value of scheme not implemented properly ',scheme
        stop
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        Pqqreg=-one-x
        ii_mqq=omx+Pqqreg*(two*lomx+L-epinv)-(one+x**2)/omx*lx
        alfax=aii/omx
        if (alfax .lt. 1d0) ii_mqq=ii_mqq+(two/omx+Pqqreg)*dlog(alfax)
        return
      endif
      
      ii_mqq=two/omx*(two*lomx+L-epinv)
      
      return
      end

***************************** Quark-Gluon *****************************
      double precision function ii_mqg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,Pqgreg,alfax,mbar
      include 'constants.f'
      include 'epinv.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial gluon-quark antenna, either

      ii_mqg=0d0
      if ((vorz .eq. 1) .or. (vorz .eq. 3)) return
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        Pqgreg=one-two*x*omx
        ii_mqg=Pqgreg*(two*lomx-lx+L-epinv)+two*x*omx
        alfax=aii/omx
        if (alfax .lt. 1d0) ii_mqg=ii_mqg+Pqgreg*dlog(alfax)
      endif
      return
      end
      
***************************** Gluon-Quark *****************************
      double precision function ii_mgq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,Pgqreg,alfax,mbar
      include 'constants.f'
      include 'epinv.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial quark-quark (--> gluon) antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)  
      
      ii_mgq=0d0
      if ((vorz .eq. 1) .or. (vorz .eq. 3)) return
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        Pgqreg=(one+omx**2)/x
        ii_mgq=Pgqreg*(two*lomx-lx+L-epinv)+x
        alfax=aii/omx
        if (alfax .lt. 1d0) ii_mgq=ii_mgq+Pgqreg*dlog(alfax)
        return
      endif

      return
      end

***************************** Gluon-Gluon *****************************
      double precision function ii_mgg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,Pggreg,alfax,mbar
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
      
      if (vorz .eq. 1) then
        ii_mgg=epinv*(epinv2-L)+half*L**2-pisqo6
        if (scheme .eq. 'tH-V') then
          return
        elseif (scheme .eq. 'dred') then
          ii_mgg=ii_mgg-1d0/6d0
          return
      else
        write(6,*) 'Value of scheme not implemented properly ',scheme
        stop
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      
      if (vorz .eq. 2) then
        Pggreg=omx/x+x*omx-one
        lx=dlog(x)
        ii_mgg=two*Pggreg*(two*lomx-lx+L-epinv)-two*lx/omx
        alfax=aii/omx
        if (alfax .lt. 1d0) ii_mgg=ii_mgg
     &   +two*(one/omx+Pggreg)*dlog(alfax)
        return
      endif
      
      ii_mgg=two*(two*lomx+L-epinv)/omx
      
      return
      end

***********************************************************************
**************************** INITIAL-FINAL ****************************
***********************************************************************

***************************** Quark-Quark *****************************
      double precision function if_mqq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,omx,lx,lomx,Pqqreg,ddilog,zp
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'alfacut.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- initial-final quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     

      mbarsq=mbar**2
      if_mqq=0d0
      if (vorz .eq. 1) then
        if_mqq=(epinv+dlog(1d0+mbarsq))*(epinv-L)+half*L**2
     &   -half*dlog(1d0+mbarsq)**2+2d0*dlog(mbarsq)*dlog(1d0+mbarsq)
     &   +2d0*ddilog(-mbarsq)+pisqo6
c--- correct form for double pole (normally zero)
        if_mqq=if_mqq-epinv**2+epinv*epinv2     
        if (scheme .eq. 'tH-V') then
          return
        elseif (scheme .eq. 'dred') then
          if_mqq=if_mqq-half
          return
      else
        write(6,*) 'Value of scheme not implemented properly ',scheme
        stop
        endif
      endif
      omx=one-x
      lomx=dlog(omx)
      zp=omx/(omx+x*mbarsq)
      if (vorz .eq. 2) then
         Pqqreg=-(1d0+x)
         lx=dlog(x)
         if_mqq=Pqqreg*(-(epinv-L)+2d0*lomx-lx-dlog(x*mbarsq+omx))
     &    +omx-2d0/omx*(lx+dlog((1d0+x*mbarsq+omx)/(1d0+mbarsq)))
         if (aif .lt. zp) then
         if_mqq=if_mqq-(two/omx*(dlog(zp*(omx+aif)/(aif*(omx+zp))))
     &    +Pqqreg*dlog(zp/aif))
         endif
      elseif (vorz .eq. 3) then
         if_mqq=2d0/omx*(-(epinv-L)+2d0*lomx-dlog(1d0+mbarsq))
      endif
      return
      end

***************************** Gluon-Gluon *****************************
      double precision function if_mgg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,mbar,mbarsq,Pggreg,ddilog,zp,omzp
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'alfacut.f'
      include 'scheme.f'

      mbarsq=mbar**2
      if_mgg=0d0

      if (vorz .eq. 1) then
          if_mgg=(epinv+dlog(1d0+mbarsq))*(epinv-L)+half*L**2
     &     -half*dlog(1d0+mbarsq)**2+2d0*dlog(mbarsq)*dlog(1d0+mbarsq)
     &     +2d0*ddilog(-mbarsq)+pisqo6
c--- correct form for double pole (normally zero)
        if_mgg=if_mgg-epinv**2+epinv*epinv2     
        if (scheme .eq. 'tH-V') then
          return
        elseif (scheme .eq. 'dred') then
          if_mgg=if_mgg-1d0/6d0
          return
       else
        write(6,*) 'Value of scheme not implemented properly ',scheme
        stop
        endif
        return
      endif
      omx=one-x
      zp=omx/(omx+x*mbarsq)
      if (vorz .eq. 2) then
         Pggreg=2d0*(omx/x-1d0+x*omx)
         lx=dlog(x)
         if_mgg=Pggreg*(-(epinv-L)+2d0*dlog(omx)-lx-dlog(x*mbarsq+omx))
     &    +2d0*mbarsq*dlog(x*mbarsq/(x*mbarsq+omx))
     &    -2d0/omx*(lx+dlog((1d0+x*mbarsq+omx)/(1d0+mbarsq)))
         if (aif .lt. zp) then
           if (aif .eq. 1d0) then
             write(6,*) 'zp > 1 in dipoles_mass.f - this is forbidden'
             stop
           endif  
         omzp=x*mbarsq/(omx+x*mbarsq)
           if_mgg=if_mgg-(two/omx*(dlog(zp*(omx+aif)/(aif*(omx+zp))))
     &     +Pggreg*dlog(zp/aif)+2d0*mbarsq*dlog(omzp/(1d0-aif)))
         endif
         return
      elseif (vorz .eq. 3) then
         if_mgg=2d0/omx*(-(epinv-L)+2d0*dlog(omx)-dlog(1d0+mbarsq))
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
      double precision function fi_mqq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,omx,ddilog
      include 'constants.f'
      include 'epinv.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- final-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     

      mbarsq=mbar**2
      if (vorz .eq. 1) then
        fi_mqq=(1d0+dlog(mbarsq/(1d0+mbarsq)))*(epinv-L)
     &        +dlog(mbarsq)+half*dlog(mbarsq)**2
     &        +half*dlog(1d0+mbarsq)**2
     &        -2d0*dlog(mbarsq)*dlog(1d0+mbarsq)
     &        -2d0*ddilog(-mbarsq)
     &        +2d0-pisq/3d0
     &        +2d0*dlog(afi)*(dlog((1d0+mbarsq)/mbarsq)-1d0)
        return
      endif
      
      omx=one-x
      
      if (vorz .eq. 2) then
        if (x .gt. 1d0-afi) then
          fi_mqq=+omx/2d0/(x*mbarsq+omx)**2
     &     +2d0/omx*(dlog((1d0+x*mbarsq+omx)*mbarsq/
     &                    ((1d0+mbarsq)*(omx+x*mbarsq))))
        else
          fi_mqq=0d0
        endif        
        return
      endif
      
      if (vorz .eq. 3) then
        if (x .gt. 1d0-afi) then
          fi_mqq=2d0/omx*(dlog((1d0+mbarsq)/(mbarsq))-1d0)
        else
          fi_mqq=0d0
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
      double precision function ff_mqq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,Lro,ro,vtijk,arg,ddilog,ff_alpha

      include 'constants.f'
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

c      if (aff .lt. 1d0) then
c        write(6,*) 'Integrated dipole routine ff_mqq does not yet'
c      write(6,*) ' support values of alpha < 1.'
c      stop
c      endif

      ff_mqq=0d0
      if (vorz .eq. 1) then
        mbarsq=mbar**2
        arg=1d0-4d0*mbarsq
        if (arg .lt. 0d0) then
            write(6,*) 'Threshold problem in ff_mqq'
            stop
        endif
        vtijk=dsqrt(arg)/(1d0-2d0*mbarsq)
        ro=dsqrt((1d0-vtijk)/(1d0+vtijk))
        Lro=dlog(ro)
        ff_mqq=1d0/vtijk*((epinv-L)*Lro-half*Lro**2-2d0*Lro*dlog(arg)
     &   +4d0*ddilog(-ro)-6d0*ddilog(1d0-ro)+pisq/3d0)
     &   +epinv-L+3d0
     &   -2d0*dlog(1d0-2d0*mbar)+dlog(1d0-mbar)+half*dlog(mbarsq)
     &   -2d0*mbar/(1d0-2d0*mbarsq)*(1d0-2d0*mbar)-mbar/(1d0-mbar)
     &   -2d0*mbarsq/(1d0-2d0*mbarsq)*dlog(mbar/(1d0-mbar))
c        Ieik=1d0/vtijk*(half*(epinv-L)*Lro-Lro*dlog(arg)
c     &  -dlog(roj)**2+pisq/6d0
c     &  +2d0*ddilog(-ro)-2d0*ddilog(1d0-ro)-ddilog(1d0-roj**2))
c        Icoll=(epinv-L)+half*Lmbarsq-2d0-2d0*dlog((1d0-mbar)**2-mbarsq)
c     & +dlog(1d0-mbar)-2d0*mbarsq/(1d0-2d0*mbarsq)*dlog(mbar/(1d0-mbar))
c     & +5d0-mbar/(1d0-mbar)-2d0*mbar*(1d0-2d0*mbar)/(1d0-2d0*mbarsq)
c        ff_mqq=2d0*Ieik+Icoll

c--- now add on the contribution for aff < 1
c        if (abs(aff-1d0) .gt. 1d-8) then
      ff_mqq=ff_mqq+ff_alpha(mbar,mbar)
c      endif 
        return
      endif
      return
      end


***************************** Quark-Quark *****************************
*   (case when emitter and spectator have different non-zero masses)  *
***********************************************************************
      double precision function ff_2mqq(x,L,mbar,sbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,sbar,sbarsq,
     & Lro,ro,roj,rok,vtijk,arg,arg2,den,ddilog

      include 'constants.f'
      include 'epinv.f'
      include 'alfacut.f'

      if (aff .lt. 1d0) then
        write(6,*) 'Integrated dipole routine ff_2mqq does not yet'
      write(6,*) ' support values of alpha < 1.'
      stop
      endif

      ff_2mqq=0d0
      if (vorz .eq. 1) then
        mbarsq=mbar**2
      sbarsq=sbar**2
        arg=1d0+mbarsq**2+sbarsq**2-2d0*(mbarsq+sbarsq+mbarsq*sbarsq)
      arg2=1d0-2d0*sbar+sbarsq-mbarsq
      den=1d0-mbarsq-sbarsq
        if (arg .lt. 0d0) then
            write(6,*) 'Threshold problem in ff_2mqq'
            stop
        endif
        vtijk=dsqrt(arg)/den
c--- some redundancy here: ro=roj*rok
      roj=dsqrt((1d0-vtijk+2d0*mbarsq/den)/(1d0+vtijk+2d0*mbarsq/den))
      rok=dsqrt((1d0-vtijk+2d0*sbarsq/den)/(1d0+vtijk+2d0*sbarsq/den))
        ro=dsqrt((1d0-vtijk)/(1d0+vtijk))
        Lro=dlog(ro)
      ff_2mqq=1d0/vtijk*((epinv-L)*Lro-dlog(roj)**2-dlog(rok)**2
     &   -2d0*Lro*dlog(1d0-(mbar+sbar)**2)
     &   +4d0*ddilog(-ro)-4d0*ddilog(1d0-ro)
     &   -ddilog(1d0-roj**2)-ddilog(1d0-rok**2)+pisq/3d0)
     &   +epinv-L+3d0
     &   +dlog(1d0-sbar)+half*dlog(mbarsq)-2d0*dlog(arg2)
     &   +(4d0*sbarsq-2d0*sbar-2d0*mbarsq*dlog(mbar/(1d0-sbar)))/den
     &   -sbar/(1d0-sbar)
        return
      endif
      
      return
      end


***************************** Quark-Quark *****************************
*    (case when the emitter is massless and spectator is massive,     *
*    plus the case when emitter is massive and spectator is massless) *
***********************************************************************
      double precision function ff_1mqq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,ddilog,afftmp
     & ,Icolla,Ieika,Icollb,Ieikb,arg,ommsq,logm,logomm,xp,yp
     & ,arg1,arg2,arg3,ypp,ypm


      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
      include 'phi.f'

      ff_1mqq=0d0 
      phi=1d0


      if (vorz .eq. 1) then
        if (mbar .gt. 1d0) then
            write(6,*) 'Problem with mbar in ff_1mqq, mbar=',mbar
            stop
        endif


      mbarsq=mbar**2
      ommsq=1d0-mbarsq
      logm=dlog(mbarsq)
      logomm=dlog(ommsq)

c--- Note: this routine uses the exact form of the eikonal integrals
c---       given in the appendix, not the symmetric combination

C----radiation from massive line with massless spectator
      afftmp=aff       
      arg=afftmp+(1d0-afftmp)*mbarsq
      Ieika=
     & half*logm*(epinv-L)-2d0*ddilog(ommsq)
     & -logm*logomm-0.25d0*logm**2
     & -dlog(afftmp)*logm-ddilog(-ommsq/mbarsq)
     & +ddilog(-afftmp*ommsq/mbarsq)
      Icolla=
     &  epinv-L+phi+2d0+(1d0+half*phi)*logm
     & +(phi-2d0)*logm/(ommsq)-2d0*logomm
     & +half*phi*(3d0*afftmp-2d0-(3d0-mbarsq)/ommsq*dlog(arg)
     & -afftmp/arg)
     & -2d0*dlog(afftmp)+2d0*dlog(arg)/ommsq


C----radiation from massless line with massive spectator
      yp=(1d0-mbar)/(1d0+mbar)
      afftmp=aff*yp       
      xp=(yp-afftmp)+dsqrt((yp-afftmp)*(1d0/yp-afftmp))
      arg1=0.25d0*(1d0-yp**2+2*xp*yp)
      arg2=half*(1d0+YP-xp)
      arg3=half*(1d0-YP+xp)
      ypp=half*(1d0+yp)
      ypm=half*(1d0-yp)

      Ieikb=0d0
     & +half*epinv*epinv2-half*epinv*L+0.25d0*L**2
     & -logomm*(epinv-L)
     & +ddilog(ommsq)-2.5d0*pisqo6+logomm**2
     & +half*dlog(arg1/(arg2*arg3))**2-dLOG(arg2/ypp)**2
     & +2d0*(dLOG(ypp)*dLOG(arg3/ypm)+dLOG(ypp/yp)*dLOG(arg1/(ypp*ypm))
     & +DDILOG(ypm/ypp)-DDILOG(arg1/ypp**2)
     & +DDILOG(arg3)-DDILOG(ypm))

      Icollb=1.5d0*(epinv-L)-3d0*dlog(1d0-mbar)+5d0-mbar/(1d0-mbar)
     & -2d0*mbar*(1d0-2d0*mbar)/ommsq
     & +1.5d0*(dLOG(YP/afftmp)-YP+afftmp)

c--- Note: extra factor of half because we include this term once for each
c---  leg, but this is the sum of both legs
      ff_1mqq=half*(2d0*(Ieika+Ieikb)+Icolla+Icollb)
    

        if (scheme .eq. 'tH-V') then
           return
        elseif (scheme .eq. 'dred') then
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
      double precision function ff_mqq0(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,ddilog
      include 'constants.f'
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

      if (aff .lt. 1d0) then
        write(6,*) 'Integrated dipole routine ff_mqq0 does not yet'
      write(6,*) ' support values of alpha < 1.'
      stop
      endif

      ff_mqq0=0d0 
      if (vorz .eq. 1) then
        if (mbar .gt. 1d0) then
            write(6,*) 'Problem with mbar in ff_mqq0, mbar=',mbar
            stop
        endif

        mbarsq=mbar**2
        ff_mqq0=epinv*epinv2/2d0
     &         +epinv*(1d0-L/2d0+dlog(mbarsq)/2d0-dlog(1d0-mbarsq))
     &         -L+3d0+L**2/4d0-pisq*5d0/12d0
     &         -ddilog(1d0-mbarsq)
     &         -dlog(mbarsq)*(dlog(mbarsq)/2d0-1d0+L)/2d0
     &         -dlog(1d0-mbarsq)*(dlog(mbarsq/(1d0-mbarsq))-L+2d0)
     &         -mbarsq/(1d0-mbarsq)*dlog(mbarsq)
      
      endif
      
      return
      end

***************************** Gluon-Gluon *****************************
*                (in this case the spectator is massive)              *
***********************************************************************
      double precision function ff_mgg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,ddilog

      include 'constants.f'
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

      if (aff .lt. 1d0) then
        write(6,*) 'Integrated dipole routine ff_mgg does not yet'
      write(6,*) ' support values of alpha < 1.'
      stop
      endif

      ff_mgg=0d0
      if (vorz .eq. 1) then
        mbarsq=mbar**2
c      ff_mgg=(epinv-L)*(epinv+dlog(mbarsq)-2d0*dlog(1d0-mbarsq))
c     &   +L**2/2d0+100d0/9d0-pisq*5d0/6d0
c     &   +2d0*b0/xn*(epinv-L)
c     &   -2d0*ddilog(1d0-mbarsq)-dlog(mbarsq)**2/2d0
c     &   -2d0*dlog(1d0-mbarsq)*dlog(mbarsq/(1d0-mbarsq))
c     &   -22d0/3d0*(dlog(1d0-mbar)+mbar/(1d0+mbar))
c     &   +4d0/3d0*mbarsq/(1-mbarsq)*dlog(2d0*mbar/(1d0+mbar))
c     &   -dfloat(nflav)/xn*4d0/3d0*(
c     &       4d0/3d0-dlog(1d0-mbar)-mbar/(1d0+mbar)
c     &      +mbarsq/(1d0-mbarsq)*dlog(2d0*mbar/(1d0+mbar)))
c--- NB: last term added by JC, 6/10/08 
        if (nfonly) then
        ff_mgg=0d0
      else
         ff_mgg=(epinv-L)*(epinv+dlog(mbarsq)-2d0*dlog(1d0-mbarsq))
     &    +L**2/2d0+100d0/9d0-pisq*5d0/6d0
     &    +2d0*11d0/6d0*(epinv-L)
     &    -2d0*ddilog(1d0-mbarsq)-dlog(mbarsq)**2/2d0
     &    -2d0*dlog(1d0-mbarsq)*dlog(mbarsq/(1d0-mbarsq))
     &    -22d0/3d0*(dlog(1d0-mbar)+mbar/(1d0+mbar))
     &    +4d0/3d0*mbarsq/(1d0-mbarsq)*dlog(2d0*mbar/(1d0+mbar))
          if (scheme .eq. 'tH-V') then
            continue ! the above is the CT in this scheme
          elseif (scheme .eq. 'dred') then
            ff_mgg=ff_mgg-1d0/3d0
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
     &   -tr*4d0/3d0*dfloat(nflav)/ca*(epinv-L)
     &   -dfloat(nflav)/ca*2d0*tr*4d0/3d0*(
     &       4d0/3d0-dlog(1d0-mbar)-mbar/(1d0+mbar)
     &      +mbarsq/(1d0-mbarsq)*dlog(2d0*mbar/(1d0+mbar)))
         endif
      endif
      
      return
      end

      
      
      
      
      
      

************************************************************************
************************************************************************
*     BELOW HERE FUNCTIONS HAVE NOT BEEN CHECKED (AND ARE NOT USED)    *
************************************************************************
************************************************************************


***************************** Gluon-Gluon *****************************
      double precision function fi_mgg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,omx
      include 'constants.f'
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


      if (vorz .eq. 1) then
        fi_mgg=two*epinv*(epinv2-L)+L**2+67d0/9d0-pisq
     &    +11d0*(epinv-L)/3d0
      endif
      
      omx=one-x
      
      if (vorz .eq. 2) then
        fi_mgg=four*dlog(two-x)/omx
        return
      endif
      
      fi_mgg=-(four*dlog(omx)+11d0/3d0)/omx
      return
      end




***************************** Quark-Gluon *****************************
      double precision function fi_mqg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,omx,rt,JaS,JaNS
      include 'constants.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--MSbar

      mbarsq=mbar**2
      if (1d0-4d0*mbarsq .lt. 0d0) then
      write(6,*) 'error in fi_mqg,(1d0-4d0*mbarsq .lt. 0d0)'
      stop
      else
      rt=dsqrt(1d0-4d0*mbarsq)
      endif
      if (vorz .eq. 1) then
CDTS 5.63
      JaS=-2d0/3d0*dlog(mbarsq)-10d0/9d0
CDTS 5.64
      JaNS=10d0/9d0*(1d0-rt)-8d0/9d0*mbarsq*rt
     & +4d0/3d0*dlog(half*(1d0+rt))
      fi_mqg=JaS+JaNS
       
      elseif (vorz .eq. 2) then
C regular
        fi_mqg=0d0
      elseif (vorz .eq. 3) then
C plus at x=x+
      omx=1d0-x
CDTS 5.62
      fi_mqg=2d0/3d0*(omx+2d0*mbarsq)/omx**2*sqrt(1d0-4d0*mbarsq/omx)
      endif
      return
      end









***************************** Quark-Gluon *****************************
      double precision function ff_mqg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,ro,arg
C-----L=Log(Qsq/musq)
      include 'constants.f'
     
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)
C --MSbar
c Id,aqg=-2/3*(epinv-L)-16/9
c      
      ff_mqg=0d0
      if (vorz .eq. 1) then
        mbarsq=mbar**2
        arg=1d0-4d0*mbarsq
        if (arg .lt. 0d0) then
        write(6,*) 'Threshold problem in ff_mqg'
        stop
        else
        ro=dsqrt(arg)  
        ff_mqg=-2d0/3d0*(2d0*dlog(mbarsq)
     &   -2d0*dlog(half*(1d0+ro))+2d0/3d0*ro*(3d0+ro**2))
        return
        endif
      endif
      return
      end

***************************** Gluon-Gluon *****************************
c      double precision function ff_mgg(x,L,mbar,vorz)
c      implicit none
c      integer vorz
c      double precision x,L,mbar,Ieik,Icoll
c      include 'constants.f'
c      include 'epinv.f'
c      include 'epinv2.f'
cc--- returns the integral of the subtraction term for an
cc--- final-initial gluon-gluon antenna, either
cc--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
cC --MSbar
cc Id,aqg=-2/3*(epinv-L)-16/9
cc Id,agg=2*epinv*(epinv-L)+L^2+11/3*(epinv-L)+100/9-[pi]^2
c      
c      ff_mgg=0d0
cCDST Eq.5.32
c      if (vorz .eq. 1) then
c      Ieik=half*epinv*(epinv2-L)+half*L**2-pisq/4d0
c      Icoll=11d0/6d0*(epinv-L)+50d0/9d0
cCDST Eq.5.36 needs to be added too - need to fix this
c      Icoll=Icoll+half*dfloat(nf)/xn*(-2d0/3d0*(epinv-L)-16d0/9d0)
c      ff_mgg=2d0*(2d0*Ieik+Icoll)
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

