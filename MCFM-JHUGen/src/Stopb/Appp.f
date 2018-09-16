      subroutine Aamp_ppp(q,mc,ms,Appp)
      implicit none
      include 'types.f'
c--- u + g  ->  c + s + d  (t-channel single-charm)
************************************************************************
*                                                                      *
* AUTHORS: R. FREDERIX AND F. TRAMONTANO                               *
* DATE  : 12/17/2008                                                    *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'stopf1inc.f'
      real(dp):: q(mxpart,4),dot,cDs,gDs,cDg,mc,ms,
     & mc2,ms2,qsq,s,t,u,xsn,xsd,xs
      complex(dp):: trc,trg,trs,trsgc,zppp,Appp

      mc2=mc**2
      ms2=ms**2

      cDs=dot(q,3,4)+mc2*dot(q,4,2)/2._dp/dot(q,3,2)
     & +ms2*dot(q,3,2)/2._dp/dot(q,4,2)
      cDg=dot(q,3,2)
      gDs=dot(q,4,2)
      qsq=mc2+ms2+2._dp*cDs+2._dp*cDg+2._dp*gDs
      s=ms2+2._dp*gDs
      t=mc2+2._dp*cDg
      u=mc2+ms2+2._dp*cDs

      xsn=(1._dp-sqrt(1._dp-4._dp*ms*mc/(u-(ms-mc)**2)))
      xsd=(1._dp+sqrt(1._dp-4._dp*ms*mc/(u-(ms-mc)**2)))
      xs=-xsn/xsd

      trg=2._dp*za(5,2)*zb(2,1)
      trs=2._dp*za(5,4)*zb(4,1)+ms**2*za(5,2)*zb(2,1)/gDs
      trc=2._dp*za(5,3)*zb(3,1)+mc**2*za(5,2)*zb(2,1)/cDg
      trsgc=2._dp*zb(1,4)*za(4,2)*zb(2,3)*za(3,5)
      zppp=za(2,4)*za(2,3)

      Appp = mc*ms2*(cDg*(ms2*trg-(2*trs+trg)*gDs)+trsgc*gDs)*tr3Xs/gDs
     &   **2

      Appp = mc*mc2*(cDg*(ms2*trg-(2*trs+trg)*gDs)+trsgc*gDs)*tr3Xc/cDg
     &   **2+Appp

      Appp = mc*(cDg*(ms2*trg-2*(trs+trg)*gDs)+trsgc*gDs)*tr1Xfs/cDg+A
     &   ppp

      Appp = mc*(cDg*(ms2*trg-2*(trs+trg)*gDs)+trsgc*gDs)*tr1Xfc/gDs+A
     &   ppp

      Appp = (-mc2*ms2*(trs+trg+trc)*gDs**3-mc2*(2*trsgc+ms2*trs)*cDg*g
     &   Ds**2+ms2*trc*cDg**3*(ms2-2*(gDs+cDs))+cDg**2*gDs*(4*mc2*(trs+t
     &   rg)*gDs-ms2*(trsgc+mc2*(3*sck-2)*trg))-2*ms2*trc*cDg**4)*lVs/(m
     &   c*cDg**2*gDs**2)/4.0_dp+Appp

      Appp = mc*(mc2*(trs+trg+trc)*gDs**3-2*(trg+trc)*cDg**2*gDs**2+(mc
     &   2*(trs+trg)-trsgc)*cDg*gDs**2+trc*cDg**3*(2*(gDs+cDs)-ms2)+2*tr
     &   c*cDg**4)*lVc/(cDg**2*gDs**2)/4.0_dp+Appp

      Appp = mc*(cDg*(ms2*trg-2*(trs+trg)*gDs)+trsgc*gDs)*lRs2/(cDg*gDs
     &   )/2.0_dp+Appp

      Appp = epinv*mc*(cDg*(ms2*trg*(2*lRs1+2*lRc1-2*sck+1)-2*(trs+trg)
     &   *gDs*(2*lRs1+2*lRc1-1))+trsgc*gDs*(2*lRs1+2*lRc1-1))/(cDg*gDs)/
     &   4.0_dp+Appp

      Appp = mc*(cDg*(ms2*trg-2*(trs+trg)*gDs)+trsgc*gDs)*lRc2/(cDg*gDs
     &   )/2.0_dp+Appp

      Appp = (mc2*cDg**2*(4*(10*trs-trg+8*trc)*gDs**3-2*(18*trg*cDs+33*
     &   trsgc-mc2*(52*trs+15*trg+3*trc)+ms2*(-3*trs+26*trg+8*trc))*gDs*
     &   *2-mc2*ms2*(10*trs+15*sck*trg-5*trg+8*trc)*gDs+9*ms2*(mc2*trg*(
     &   cDs+ms2)+trsgc*(cDs+mc2)))+cDg**3*(4*mc2*(29*trs+15*trg+3*trc)*
     &   gDs**2-2*mc2*(ms2*(10*trs+15*sck*trg-5*trg+6*trc)+3*mc2*trc)*gD
     &   s+6*(3*ms2*trsgc+mc2*ms2*(3*trg+trc)-mc2**2*trc)*cDs+mc2*ms2*(2
     &   7*trsgc+mc2*(-26*trs-17*trg+3*trc)+18*ms2*trg+5*ms2*trc))+mc2*c
     &   Dg*gDs**2*(2*(-9*trsgc+mc2*(17*trs-5*trg+13*trc)+3*ms2*(trs+trg
     &   +trc))*gDs-mc2*(36*trg*cDs+54*trsgc+ms2*(-3*trs+44*trg+8*trc))-
     &   18*trsgc*cDs+23*mc2**2*trs)+mc2**2*gDs**2*((-9*trsgc+mc2*(7*trs
     &   -4*trg+5*trc)+3*ms2*(trs+trg+trc))*gDs-9*mc2*trg*(cDs+ms2)-9*tr
     &   sgc*(cDs+mc2))-2*cDg**4*(-4*ms2*trc*gDs+6*mc2*trc*gDs+6*(mc-ms)
     &   *(ms+mc)*trc*cDs-9*ms2*trsgc+26*mc2*ms2*trs+17*mc2*ms2*trg-5*ms
     &   2**2*trc-6*mc2*ms2*trc+3*mc2**2*trc)-12*(mc-ms)*(ms+mc)*trc*cDg
     &   **5)/(mc*cDg**2*(2*cDg+mc2)*gDs**2)/1.2e+1_dp+Appp

      Appp = 2*LsA*mc*(-mc2*trsgc*gDs**3-cDg**2*gDs*(cDs*(2*(trs+trg)*g
     &   Ds-ms2*trg)+ms2*trsgc)+cDg*gDs**2*(mc2*(2*trs+trg)*gDs+trsgc*cD
     &   s-mc2*ms2*trg)+ms2*cDg**3*((2*trs+trg)*gDs-ms2*trg))/(cDg**2*gD
     &   s**2)+Appp

      Appp = B0cgsf*(cDg**2*(mc2*(2*(trs+trc)*gDs**2-ms2*(trs+trc)*gDs+
     &   ms2*trg*(cDs+ms2))+ms2*trsgc*(gDs+cDs+mc2))-mc2*gDs**2*((trsgc-
     &   ms2*(trs+trg+trc)+mc2*trg)*gDs+mc2*trg*(cDs+ms2)+trsgc*(cDs+mc2
     &   ))+mc2*cDg*gDs**2*(2*(trs+trc)*gDs-2*trg*cDs-trsgc+ms2*(trs-3*t
     &   rg-trc)+2*mc2*trs-mc2*trg)+cDg**3*(2*(ms**2-mc**2)*trc*gDs+2*(m
     &   s**2-mc**2)*trc*cDs+ms2*(trsgc+mc2*(-3*trs-2*trg+trc)))-2*(mc-m
     &   s)*(ms+mc)*trc*cDg**4)/(mc*cDg**2*gDs**2)/4.0_dp+Appp

      Appp = BfunX*(-cDg**3*(4*trc*gDs**2+2*trc*cDs*(4*gDs+mc2)+2*mc2*t
     &   rc*gDs+4*trc*cDs**2+mc2*ms2*(trs+trg-trc))-mc2*cDg*gDs**2*(2*tr
     &   s*(gDs+cDs)-ms2*(-trs+trg+trc))-2*trc*cDg**4*(4*(gDs+cDs)+mc2)-
     &   mc2*(trs+trg+trc)*gDs**3*(2*(gDs+cDs)+ms2)-mc2*cDg**2*gDs*(ms2*
     &   (trs+trg+trc)-2*(trg+trc)*gDs)-4*trc*cDg**5)/(mc*cDg**2*gDs**2)
     &   /4.0_dp+Appp

      Appp = Appp-ms2*tr3c00fs*(mc2*(trg*(2*gDs+cDs+ms2)+trs*gDs)+cDg*
     &   (4*trc*gDs+trsgc-mc2*(2*trs+trg)+3*ms2*trc)+trsgc*(cDs+mc2))/(m
     &   c*gDs**2)

      Appp = mc*tr3s00ft*(mc2*(trg*(gDs+cDs+ms2)+trs*gDs)+cDg*(2*trg*(g
     &   Ds+cDs)+2*trs*(gDs-mc2)+3*trsgc+ms2*(5*trg+3*trc))+trsgc*(gDs+c
     &   Ds+mc2)-2*trs*cDg**2)/cDg**2+Appp

      Appp = epinv**2*mc*(cDg*(ms2*trg-2*(trs+trg)*gDs)+trsgc*gDs)/(cDg
     &   *gDs)/2.0_dp+Appp

      Appp = Appp-3*ms2*tr3c001fs*(trc*cDg*(2*gDs+ms2)+mc2*(trs+trg)*g
     &   Ds)/(mc*gDs**2)

      Appp = 3*mc*tr3s001ft*(trs*(2*cDg+mc2)*gDs+ms2*(trg+trc)*cDg)/cDg
     &   **2+Appp

      Appp = 3*mc*tr3s002ft*(2*cDg+mc2)*((trg+trc)*gDs+trs*cDg)/cDg**2+
     &   Appp

      Appp = Appp-3*mc*ms2*tr3c002fs*(trc*gDs+(trs+trg)*cDg)/gDs**2

      Appp = ls*(ms2*trsgc+mc2*ms2*trg)/(4*mc*gDs+2*mc*ms2)+Appp

      Appp = lc*mc*(mc2*trsgc/(2*cDg+mc2)**2+trg)/2.0_dp+Appp

      Appp=Appp/zppp

      return
      end
