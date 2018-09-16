      subroutine Aamp_ppm(q,mc,ms,Appm)
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
      complex(dp):: trc,trg,trs,trsgc,zppm,Appm

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
      zppm=za(2,3)**2*zb(4,3)

      Appm = mc*ms*(cDg**2*(ms2*trg*(gDs+ms2)-2*trs*gDs**2)+mc2*ms2*trg
     &   *gDs**2+cDg*gDs*(trsgc*gDs-2*ms2*trg*cDs))*tr3Xs/gDs**3

      Appm = mc*mc2*ms*trg*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)*tr3Xc/
     &   (cDg**2*gDs)+Appm

      Appm = mc*ms*trg*(mc2*gDs/cDg+ms2*cDg/gDs-2*cDs)*tr1Xfs+Appm

      Appm = mc*ms*trg*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)*tr1Xfc/gDs
     &   **2+Appm

      Appm = Appm-ms*(2*mc2*ms2*(trs+trg+trc)*cDg*gDs**3+mc2*gDs**4*(t
     &   rs*(mc2-2*(gDs+cDs))+ms2*(-trs+trg+trc))+cDg**3*gDs*(ms2*(trsgc
     &   +mc2*trg)-2*(-trsgc+mc2*trs+ms2*trc)*gDs)+mc2*cDg**2*gDs**2*(2*
     &   (trs+trg)*gDs+trsgc)-ms2*trc*cDg**4*(ms2-2*cDs)+2*ms2*trc*cDg**
     &   5)*lVs/(mc*cDg**2*gDs**3)/4.0_dp

      Appm = mc*(mc2*(ms2*(2*trs+trg+2*trc)-trsgc)*cDg*gDs**3-2*ms2*(tr
     &   g+trc)*cDg**3*gDs**2+mc2*gDs**4*(trs*(mc2-2*(gDs+cDs))+ms2*(-tr
     &   s+trg+trc))+cDg**2*gDs**2*(2*(-trsgc+mc2*trs+ms2*trc)*gDs-ms2*(
     &   trsgc+2*mc2*trg))-ms2*trc*cDg**4*(ms2-2*cDs)+2*ms2*trc*cDg**5)*
     &   lVc/(ms*cDg**2*gDs**3)/4.0_dp+Appm

      Appm = mc*ms*trg*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)*lRs2/(cDg*
     &   gDs**2)/2.0_dp+Appm

      Appm = epinv*mc*ms*trg*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)*(2*l
     &   Rs1+2*lRc1-1)/(cDg*gDs**2)/4.0_dp+Appm

      Appm = mc*ms*trg*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)*lRc2/(cDg*
     &   gDs**2)/2.0_dp+Appm

      Appm = (mc2*cDg**3*(-8*trs*gDs**4+4*ms2*(-13*trs+trg+6*trc)*gDs**
     &   3-2*ms2*(33*trsgc+mc2*(-27*trs+22*trg+13*trc)+29*ms2*trg-18*ms2
     &   *trc)*gDs**2+ms2*(-18*(ms2*(trg+trc)-mc2*(trs+2*trg))*cDs+mc2*(
     &   9*trsgc+16*ms2*trs+7*ms2*trg-16*ms2*trc)-18*ms2*trsgc)*gDs+9*ms
     &   2**2*((trsgc+mc2*trg)*cDs+mc2*(trsgc+ms2*trg)))-mc2*cDg**2*gDs*
     &   (4*(-9*trsgc+3*mc2*trs+8*ms2*trg+17*ms2*trc)*gDs**3-9*cDs*(4*(t
     &   rsgc-mc2*trs+ms2*(2*trg+trc))*gDs**2-2*ms2*(trsgc+mc2*trg)*gDs-
     &   mc2*ms2*(trsgc+ms2*trg))-2*ms2*(9*trsgc+6*ms2*(trs+trg+trc)+mc2
     &   *(-35*trs+23*trg+19*trc))*gDs**2+mc2*ms2*(54*trsgc-18*mc2*trs+4
     &   7*ms2*trg+23*mc2*trg-18*ms2*trc)*gDs+9*mc2*ms2**2*(trsgc+mc2*tr
     &   g))+mc2*cDg*gDs**2*(12*(mc**2-ms**2)*trs*gDs**3+3*cDs*(4*(mc**2
     &   -ms**2)*trs*gDs**2+6*mc2*(2*trsgc-mc2*trs+ms2*(3*trg+trc))*gDs-
     &   3*mc2*ms2*(trsgc+mc2*trg))+2*(mc2*(18*trsgc+ms2*(6*trs-28*trg-4
     &   6*trc))+3*mc2**2*trs+3*ms2**2*(-trs+trg+trc))*gDs**2+mc2*ms2*(2
     &   7*trsgc+6*ms2*(trs+trg+trc)+mc2*(-22*trs+37*trg+10*trc))*gDs-9*
     &   mc2**2*ms2*(trsgc+ms2*trg))+ms2*cDg**4*(4*mc2*(9*trs+trg-14*trc
     &   )*gDs**2+(2*mc2*(18*trsgc+16*ms2*trs-17*mc2*trs+25*ms2*trg-8*mc
     &   2*trg-18*ms2*trc)-36*(-trsgc-mc2*trs-2*mc2*trg+ms2*trc)*cDs)*gD
     &   s+6*(3*ms2*trsgc+mc2*ms2*(3*trg+trc)-mc2**2*trc)*cDs+mc2*ms2*(2
     &   7*trsgc+mc2*(-26*trs-17*trg+3*trc)+18*ms2*trg+5*ms2*trc))+mc2**
     &   2*gDs**3*(6*(mc**2-ms**2)*trs*gDs**2+(mc2*(6*trs*cDs+9*trsgc+ms
     &   2*(6*trs-20*trg-29*trc))+3*ms2*(ms2*(-trs+trg+trc)-2*trs*cDs)+5
     &   *mc2**2*trs)*gDs+9*mc2*((trsgc+ms2*trg)*cDs+ms2*(trsgc+mc2*trg)
     &   ))-2*ms2*cDg**5*(4*trc*gDs**2+2*(-9*trsgc+17*mc2*trs+8*mc2*trg+
     &   2*ms2*trc)*gDs+6*(mc**2-ms**2)*trc*cDs-9*ms2*trsgc+26*mc2*ms2*t
     &   rs+17*mc2*ms2*trg-5*ms2**2*trc-6*mc2*ms2*trc+3*mc2**2*trc)-12*(
     &   mc2-ms2)*ms2*trc*cDg**6)/(mc*ms*cDg**2*(2*cDg+mc2)*gDs**3)/1.2d
     &   +1+Appm

      Appm = Appm-2*LsA*mc*ms*(mc2**2*trg*gDs**4-3*mc2*trg*cDg*cDs*gDs
     &   **3+cDg**4*(ms2*trg*(gDs+ms2)-2*trs*gDs**2)+2*trg*cDg**2*(cDs**
     &   2+mc2*ms2)*gDs**2+cDg**3*gDs*(trsgc*gDs-3*ms2*trg*cDs))/(cDg**2
     &   *gDs**3)

      Appm = Appm-ms*tr3c00fs*(cDg*(mc2*((6*trs+8*trg)*gDs**2+(trsgc+6
     &   *ms2*trs+8*ms2*trg)*gDs+ms2*(trsgc+ms2*trg))+cDs*(2*(trsgc+mc2*
     &   trs+2*mc2*trg-ms2*trc)*gDs+ms2*(trsgc+mc2*trg)))+cDg**2*(8*trc*
     &   gDs**2+2*(trsgc-mc2*trs+4*ms2*trc)*gDs+ms2*(trsgc-mc2*(2*trs+tr
     &   g)+3*ms2*trc))+mc2*gDs*(6*trc*gDs**2+(-trsgc-2*mc2*trg+5*ms2*tr
     &   c)*gDs-ms2*trg*(cDs+mc2)+trsgc*(-cDs-ms2)))/(mc*gDs**3)

      Appm = 3*mc*ms*tr3s001ft*((trg+trc)*((2*cDg+mc2)*gDs**2+ms2*cDg**
     &   2)+2*trs*cDg*(2*cDg+mc2)*gDs)/(cDg**2*gDs)+Appm

      Appm = B0cgsf*(-mc2*cDg**2*gDs*(-2*(trsgc-mc2*trs+ms2*(trg-trs))*
     &   gDs**2+ms2*(2*trsgc-2*mc2*trs+3*ms2*trg+mc2*trg-2*ms2*trc)*gDs+
     &   ms2*trsgc*(cDs+ms2)+ms2**2*trg*(cDs+mc2))+mc2*cDg*gDs**2*(2*(tr
     &   sgc-ms2*(trg+2*trc))*gDs**2+cDs*(2*(trsgc-mc2*trs+ms2*(2*trg+tr
     &   c))*gDs-ms2*(trsgc+mc2*trg))+(ms2+mc2)*trsgc*gDs+2*ms2*(ms2*(tr
     &   s+trg+trc)+mc2*(trg-2*trs))*gDs-mc2*ms2*(trsgc+ms2*trg))+ms2*cD
     &   g**3*(-2*(-trsgc+mc2*(trc-trg)+ms2*trc)*gDs**2+(2*(trsgc+mc2*(t
     &   rs+2*trg)-ms2*trc)*cDs+mc2*(trsgc+2*ms2*(trs+2*trg-trc))+ms2*tr
     &   sgc)*gDs+ms2*(mc2*trg*(cDs+ms2)+trsgc*(cDs+mc2)))+mc2*gDs**3*(2
     &   *(mc-ms)*(ms+mc)*trs*gDs**2+(mc2*(2*trs*cDs+trsgc+ms2*(2*trs-3*
     &   trg-4*trc))+ms2*(ms2*(-trs+trg+trc)-2*trs*cDs))*gDs+mc2*(trsgc*
     &   (cDs+ms2)+ms2*trg*(cDs+mc2)))+ms2*cDg**4*(mc2*(ms2*(-3*trs-2*tr
     &   g+trc)-2*(2*trs+trg)*gDs)+trsgc*(2*gDs+ms2)+2*(ms**2-mc**2)*trc
     &   *cDs)+2*ms2*(ms2-mc2)*trc*cDg**5)/(mc*ms*cDg**2*gDs**3)/4.0_dp+
     &   Appm

      Appm = mc*tr3s00ft*(cDg*((-2*trsgc+4*mc2*trs+6*ms2*trg+8*ms2*trc)
     &   *gDs**2+cDs*(ms2*(trsgc+mc2*trg)-2*(trsgc-mc2*trs+2*ms2*trg+ms2
     &   *trc)*gDs)-ms2*(trsgc-6*mc2*trs+2*mc2*trg)*gDs+mc2*ms2*(trsgc+m
     &   s2*trg))+cDg**2*(4*trs*gDs**2+4*ms2*(3*trs+trg+trc)*gDs+ms2*(2*
     &   trg*cDs+3*trsgc-2*mc2*trs+5*ms2*trg+3*ms2*trc))+mc2*gDs*(ms2*(5
     &   *trc*gDs-trg*(-4*gDs+cDs+mc2))-trsgc*(gDs+cDs+ms2))-2*ms2*trs*c
     &   Dg**3)/(ms*cDg**2*gDs)+Appm

      Appm = Appm-3*mc*ms*tr3c002fs*(mc2*(trs+trg)*gDs**2+2*trc*cDg*gD
     &   s*(2*gDs+ms2)+(trs+trg)*cDg**2*(2*gDs+ms2))/gDs**3

      Appm = 3*mc*tr3s002ft*(2*cDg+mc2)*(mc2*trs*gDs**2+2*cDg*gDs*(trs*
     &   gDs+ms2*(trg+trc))+ms2*trs*cDg**2)/(ms*cDg**2*gDs)+Appm

      Appm = Appm-3*ms*tr3c001fs*(2*gDs+ms2)*(mc2*trc*gDs**2+trc*cDg**
     &   2*(2*gDs+ms2)+2*mc2*(trs+trg)*cDg*gDs)/(mc*gDs**3)

      Appm = epinv**2*mc*ms*trg*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)/(
     &   cDg*gDs**2)/2.0_dp+Appm

      Appm = BfunX*(-ms2*cDg**4*(mc2*(2*(trs+trg)*gDs+ms2*(trs+trg-trc)
     &   )+2*trc*cDs*(2*gDs+mc2)+4*trc*cDs**2)-mc2*gDs**4*(2*(gDs+cDs)+m
     &   s2)*(ms2*(-trs+trg+trc)-2*trs*(gDs+cDs))-2*mc2*cDg*gDs**3*(ms2*
     &   (trs+trg+trc)-trs*gDs)*(2*(gDs+cDs)+ms2)-2*mc2*ms2*cDg**2*gDs**
     &   2*((trs+trc)*gDs-trg*(cDs+ms2))-2*mc2*ms2*cDg**3*gDs*((trs+trc)
     &   *gDs+ms2*(trs+trg+trc))-2*ms2*trc*cDg**5*(2*gDs+4*cDs+mc2)-4*ms
     &   2*trc*cDg**6)/(mc*ms*cDg**2*gDs**3)/4.0_dp+Appm

      Appm = ls*ms*(cDg*(2*(trsgc+mc2*trs-ms2*trc)*gDs+ms2*(trsgc-mc2*t
     &   rg))+mc2*gDs*(2*trg*gDs-trsgc+ms2*trg))/(mc*gDs*(2*gDs+ms2))/2.
     &   0_dp+Appm

      Appm = lc*mc*(4*cDg**2*(trsgc*gDs-mc2*trs*gDs+mc2*ms2*trg)+mc2*cD
     &   g*(ms2*(trsgc+mc2*trg)-2*(-2*trsgc+mc2*trs+ms2*trc)*gDs)+mc2**2
     &   *(trsgc+ms2*trg)*gDs+4*ms2*trg*cDg**3)/(ms*(2*cDg+mc2)**2*gDs)/
     &   2.0_dp+Appm

      Appm=Appm/zppm

      return
      end
