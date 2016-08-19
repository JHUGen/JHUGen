      subroutine qq_tbg_z(p,z)
************************************************************************
*     Virtual ct's for s-channel single top + jet                      *
*                                                                      *
*     q(p1) + g(p2) -> t(p3) + b(p4) + q'(p5)                          *      
*                                                                      *
*     Author: J. Campbell, June 20, 2008                               *
*                                                                      *
************************************************************************
*                                                                      *
*     IMPORTANT NOTE!                                                  *
*                                                                      *
*     For now, we only include radiation from heavy quark line         *
*      and virtual corrections on that line too                        *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'agq.f'
      include 'PR_new.f'
      include 'colstruc.f'
      include 'stopscales.f'
      include 'breit.f'
      double precision z,p(mxpart,4),metric,Q34sq,Q35sq,Q45sq,
     & xl34,xl35,xl45,
     & mbar35,mbar45,mbar34a,mbar34b,
     . ff_2mqq,ff_mqq0,ff_mgg,
     . ason4pi_H
      integer is,nu
          
      Q34sq=0d0
      Q35sq=0d0
      Q45sq=0d0
      metric=1d0
      do nu=4,1,-1
      Q34sq=Q34sq+metric*(p(3,nu)+p(4,nu))**2
      Q35sq=Q35sq+metric*(p(3,nu)+p(5,nu))**2
      Q45sq=Q45sq+metric*(p(4,nu)+p(5,nu))**2
      metric=-1d0
      enddo

CDST (5.5)
      mbar34a=mass2/dSqrt(Q34sq)
      mbar34b=mass3/dSqrt(Q34sq)
      mbar35 =mass2/dSqrt(Q35sq)
      mbar45 =mass3/dSqrt(Q45sq)
      
      xl34=dlog(Q34sq/renscale_H**2)
      xl35=dlog(Q35sq/renscale_H**2)
      xl45=dlog(Q45sq/renscale_H**2)

      ason4pi_H=as_H/fourpi
C      ason4pi_L=as_L/fourpi

c--- The variables R and P provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (Q1(a,b,c,is) and leg 2 (Q2(a,b,c,is)
c--- In each case the parton labelling is Using the normal QM notation of putting
c--- everything backward
c---       emitted line after emission =   a
c---       emitter before emission     =   b
c---       spectator                   =   c
c--- There is no label for he or she who is emitted.
c--- Note that in general each piece will be composed of many different
c--- dipole contributions

      do is=1,3

c--- counterterms for extra gluon emission on the heavy line
c      Q1(q,q,a,is)=ason4pi_L*two*cf*(ii_qq(z,xl12,is)+ii_qq(z,xl12,is))

c--- leading colour (Ca) goes in (0,0)
c--- subleading colour (2*Cf-Ca) goes in (0,1)
c--- Tr terms go in (1,0)
      caonly=.true.
      nfonly=.false.
      Q1(g,g,g,is)=ason4pi_H*(
     .       ca*(ff_mgg(z,xl35,mbar35,is)/2d0+ff_mqq0(z,xl35,mbar35,is))
     .      +ca*(ff_mgg(z,xl45,mbar45,is)/2d0+ff_mqq0(z,xl45,mbar45,is))
     .                       )
      Q1(g,g,q,is)=ason4pi_H*(
     .      +(2d0*cf-ca)*
     .          (ff_2mqq(z,xl34,mbar34a,mbar34b,is)
     .          +ff_2mqq(z,xl34,mbar34b,mbar34a,is))
     .                       )
      
      caonly=.false.
      nfonly=.true.
      Q1(q,q,g,is)=ason4pi_H*(
     .       ca*(ff_mgg(z,xl35,mbar35,is)/2d0)
     .      +ca*(ff_mgg(z,xl45,mbar45,is)/2d0)
     .                       )
      enddo
      
      return
      end

