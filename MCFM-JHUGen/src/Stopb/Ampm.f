      subroutine Aamp_mpm(q,mc,ms,Ampm)
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
      complex(dp):: trc,trg,trs,trsgc,zmpm,Ampm

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
      zmpm=za(4,3)*zb(2,4)**2

      Ampm = -mc*ms*ms2*trg*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)*tr3Xs
     &   /(cDg*gDs**2)

      Ampm = Ampm-mc*ms*(mc2**2*trg*gDs**2+cDg**2*(gDs*(trsgc-2*trc*gD
     &   s)+mc2*ms2*trg)+mc2*trg*cDg*gDs*(gDs-2*cDs))*tr3Xc/cDg**3

      Ampm = Ampm-mc*ms*trg*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)*tr1X
     &   fs/cDg**2

      Ampm = mc*ms*trg*(-mc2*gDs/cDg-ms2*cDg/gDs+2*cDs)*tr1Xfc+Ampm

      Ampm = ms*(2*mc2*ms2*(trs+trg+trc)*cDg*gDs**3+mc2*gDs**4*(trs*(mc
     &   2-2*(gDs+cDs))+ms2*(-trs+trg+trc))+cDg**3*gDs*(2*(trsgc+mc2*(2-
     &   3*sck)*trs-ms2*trc)*gDs+ms2*(trsgc+mc2*(3*sck-2)*trg))+mc2*cDg*
     &   *2*gDs**2*(2*(trs+trg)*gDs+(3*sck-2)*trsgc)-ms2*trc*cDg**4*(ms2
     &   -2*cDs)+2*ms2*(trs+trg+trc)*cDg**5)*lVs/(mc*cDg**3*gDs**2)/4.0d
     &   +0+Ampm

      Ampm = mc*(mc2*(trsgc-ms2*(2*trs+trg+2*trc))*cDg*gDs**3+2*ms2*(tr
     &   g+trc)*cDg**3*gDs**2+mc2*gDs**4*(trs*(2*(gDs+cDs)-mc2)-ms2*(-tr
     &   s+trg+trc))+cDg**2*gDs**2*(ms2*(trsgc+2*mc2*trg)-2*(-trsgc+mc2*
     &   trs+ms2*trc)*gDs)+ms2*trc*cDg**4*(ms2-2*cDs)-2*ms2*(trs+trg+trc
     &   )*cDg**5)*lVc/(ms*cDg**3*gDs**2)/4.0_dp+Ampm

      Ampm = Ampm-mc*ms*trg*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)*lRs2
     &   /(cDg**2*gDs)/2.0_dp

      Ampm = Ampm-epinv*mc*ms*(cDg**2*(ms2*trg*(2*lRs1+2*lRc1-2*sck+1)
     &   +4*(sck-1)*trs*gDs)-2*cDg*gDs*(trg*cDs*(2*lRs1+2*lRc1-1)+(sck-1
     &   )*trsgc)+mc2*trg*gDs**2*(2*lRs1+2*lRc1-1))/(cDg**2*gDs)/4.0_dp

      Ampm = Ampm-mc*ms*trg*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)*lRc2
     &   /(cDg**2*gDs)/2.0_dp

      Ampm = 2*LsA*mc*ms*(mc2**2*trg*gDs**4+cDg**2*gDs**2*(gDs*(trsgc-2
     &   *trc*gDs)+2*trg*(cDs**2+mc2*ms2))+mc2*trg*cDg*gDs**3*(gDs-3*cDs
     &   )-3*ms2*trg*cDg**3*cDs*gDs+ms2**2*trg*cDg**4)/(cDg**3*gDs**2)+
     &   Ampm

      Ampm = (mc2*cDg**2*gDs*((2*gDs+ms2)*(4*trs*gDs**3-2*ms2*(-13*trs+
     &   trg+9*trc)*gDs**2+ms2*(-18*mc2*trs+29*ms2*trg+23*mc2*trg-18*ms2
     &   *trc)*gDs+9*ms2**2*trg*(cDs+mc2))+3*ms2*trsgc*(2*(5*sck+6)*gDs*
     &   *2+3*cDs*(2*gDs+ms2)+ms2*(5*sck+13)*gDs+3*ms2**2))-ms2*cDg**3*(
     &   4*mc2*(3*(5*sck-3)*trs+trg-13*trc)*gDs**3+(2*mc2*(9*trsgc+ms2*(
     &   5*(3*sck+2)*trs+(44-15*sck)*trg-29*trc))-36*(-trsgc-mc2*(trs+2*
     &   trg)+ms2*trc)*cDs)*gDs**2+ms2*(mc2*(27*trsgc+ms2*(16*trs-15*sck
     &   *trg+58*trg-16*trc))-18*(-2*trsgc-mc2*(trs+3*trg)+ms2*trc)*cDs)
     &   *gDs+9*ms2**2*(mc2*trg*(cDs+ms2)+trsgc*(cDs+mc2)))-mc2*cDg*gDs*
     &   *2*(2*gDs+ms2)*(-2*(-9*trsgc+2*mc2*trs+8*ms2*trg+17*ms2*trc)*gD
     &   s**2+9*cDs*(2*(trsgc-mc2*trs+ms2*(2*trg+trc))*gDs-ms2*(trsgc+mc
     &   2*trg))+ms2*(9*trsgc+6*ms2*(trs+trg+trc)+mc2*(-22*trs+19*trg+10
     &   *trc))*gDs-9*mc2*ms2*(trsgc+ms2*trg))-mc2*gDs**3*(2*gDs+ms2)*(6
     &   *(mc-ms)*(ms+mc)*trs*gDs**2+(mc2*(6*trs*cDs+9*trsgc+ms2*(6*trs-
     &   20*trg-29*trc))+3*ms2*(ms2*(-trs+trg+trc)-2*trs*cDs)+5*mc2**2*t
     &   rs)*gDs+9*mc2*(trsgc*(cDs+ms2)+ms2*trg*(cDs+mc2)))+ms2*cDg**4*(
     &   2*gDs+ms2)*(4*trc*gDs**2+2*(-9*trsgc+17*mc2*trs+8*mc2*trg+2*ms2
     &   *trc)*gDs+6*(mc-ms)*(ms+mc)*trc*cDs+ms2*(-9*trsgc+mc2*(26*trs+1
     &   7*trg-3*trc)-5*ms2*trc))-2*ms2*cDg**5*(2*gDs+ms2)*(8*(trs+trg)*
     &   (gDs+cDs)+ms2*(3*trc-5*(trs+trg))-3*mc2*(trs+trg+trc))-16*ms2*(
     &   trs+trg)*cDg**6*(2*gDs+ms2))/(mc*ms*cDg**3*gDs**2*(2*gDs+ms2))/
     &   1.2e+1_dp+Ampm

      Ampm = ms*tr3c00fs*(cDg*(mc2*(8*(trs+trg)*gDs**2+(trsgc+6*ms2*trs
     &   +8*ms2*trg)*gDs+ms2*(trsgc+ms2*trg))+cDs*(2*(trsgc+mc2*trs+2*mc
     &   2*trg-ms2*trc)*gDs+ms2*(trsgc+mc2*trg)))+mc2*gDs*(8*trc*gDs**2+
     &   (-2*trg*(cDs+mc2)-3*trsgc+5*ms2*trc)*gDs-ms2*trg*(cDs+mc2)+trsg
     &   c*(-cDs-ms2))+cDg**2*(8*trc*gDs**2+2*(trsgc-mc2*trs+4*ms2*trc)*
     &   gDs+ms2*(trsgc-mc2*(2*trs+trg)+3*ms2*trc))-6*(trs+trg)*cDg**3*(
     &   gDs+ms2))/(mc*cDg*gDs**2)+Ampm

      Ampm = Ampm-ls*ms*(mc2*gDs*(trg*(2*gDs+ms2)**2+ms2*trsgc)-ms2*cD
     &   g*(4*trc*gDs**2+2*(mc2*trs+ms2*trc)*gDs-mc2*ms2*trg)+trsgc*cDg*
     &   (2*gDs+ms2)**2)/(mc*cDg*(2*gDs+ms2)**2)/2.0_dp

      Ampm = tr3s001ft*(-3*mc*ms*(trg+trc)*((2*cDg+mc2)*gDs**2+ms2*cDg*
     &   *2)-6*mc*ms*trs*cDg*(2*cDg+mc2)*gDs)/cDg**3+Ampm

      Ampm = B0cgsf*(mc2*cDg**2*gDs*(-2*(trsgc-mc2*trs+ms2*(trg-trs))*g
     &   Ds**2+ms2*(2*trsgc-2*mc2*trs+3*ms2*trg+mc2*trg-2*ms2*trc)*gDs+m
     &   s2*trsgc*(cDs+ms2)+ms2**2*trg*(cDs+mc2))-mc2*cDg*gDs**2*(2*(trs
     &   gc-ms2*(trg+2*trc))*gDs**2+cDs*(2*(trsgc-mc2*trs+ms2*(2*trg+trc
     &   ))*gDs-ms2*(trsgc+mc2*trg))+(ms2+mc2)*trsgc*gDs+2*ms2*(ms2*(trs
     &   +trg+trc)+mc2*(trg-2*trs))*gDs-mc2*ms2*(trsgc+ms2*trg))-ms2*cDg
     &   **3*(-2*(-trsgc+mc2*(trc-trg)+ms2*trc)*gDs**2+(2*(trsgc+mc2*(tr
     &   s+2*trg)-ms2*trc)*cDs+mc2*(trsgc+2*ms2*(trs+2*trg-trc))+ms2*trs
     &   gc)*gDs+ms2*(mc2*trg*(cDs+ms2)+trsgc*(cDs+mc2)))-mc2*gDs**3*(2*
     &   (mc-ms)*(ms+mc)*trs*gDs**2+(mc2*(2*trs*cDs+trsgc+ms2*(2*trs-3*t
     &   rg-4*trc))+ms2*(ms2*(-trs+trg+trc)-2*trs*cDs))*gDs+mc2*(trsgc*(
     &   cDs+ms2)+ms2*trg*(cDs+mc2)))-ms2*cDg**4*(mc2*(ms2*(-3*trs-2*trg
     &   +trc)-2*(2*trs+trg)*gDs)+trsgc*(2*gDs+ms2)+2*(ms**2-mc**2)*trc*
     &   cDs)-2*ms2*cDg**5*((trs+trg)*(gDs+cDs)-mc2*(trs+trg+trc)+ms2*tr
     &   c)-2*ms2*(trs+trg)*cDg**6)/(mc*ms*cDg**3*gDs**2)/4.0_dp+Ampm

      Ampm = mc*tr3s00ft*(-cDg*((-2*trsgc+4*mc2*trs+6*ms2*trg+8*ms2*trc
     &   )*gDs**2+cDs*(ms2*(trsgc+mc2*trg)-2*(trsgc-mc2*trs+2*ms2*trg+ms
     &   2*trc)*gDs)-ms2*(trsgc-6*mc2*trs+2*mc2*trg)*gDs+mc2*ms2*(trsgc+
     &   ms2*trg))-cDg**2*(4*trs*gDs**2+2*ms2*(6*trs+2*trg+3*trc)*gDs+ms
     &   2*(trsgc-2*mc2*trs+5*ms2*trg+3*ms2*trc))-mc2*gDs*(ms2*(5*trc*gD
     &   s-trg*(-4*gDs+cDs+mc2))-trsgc*(gDs+cDs+ms2)))/(ms*cDg**3)+Ampm

      Ampm = 3*ms*tr3c002fs*(mc2**2*(trs+trg)*gDs**2+2*mc2*trc*cDg*gDs*
     &   (2*gDs+ms2)+mc2*(trs+trg)*cDg**2*(2*gDs+ms2)-2*(trs+trg)*cDg**3
     &   *cDs-2*(trs+trg)*cDg**4)/(mc*cDg*gDs**2)+Ampm

      Ampm = Ampm-3*mc*tr3s002ft*(2*cDg+mc2)*(mc2*trs*gDs**2+2*cDg*gDs
     &   *(trs*gDs+ms2*(trg+trc))+ms2*trs*cDg**2)/(ms*cDg**3)

      Ampm = 3*ms*tr3c001fs*(2*gDs+ms2)*(mc2*trc*gDs**2+trc*cDg**2*(2*g
     &   Ds+ms2)+2*mc2*(trs+trg)*cDg*gDs-2*(trs+trg)*cDg**3)/(mc*cDg*gDs
     &   **2)+Ampm

      Ampm = Ampm-epinv**2*mc*ms*trg*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg
     &   **2)/(cDg**2*gDs)/2.0_dp

      Ampm = BfunX*(ms2*cDg**4*(mc2*(2*(trs+trg)*gDs+ms2*(trs+trg-trc))
     &   +2*trc*cDs*(2*gDs+mc2)+4*trc*cDs**2)+mc2*gDs**4*(2*(gDs+cDs)+ms
     &   2)*(ms2*(-trs+trg+trc)-2*trs*(gDs+cDs))+2*mc2*cDg*gDs**3*(ms2*(
     &   trs+trg+trc)-trs*gDs)*(2*(gDs+cDs)+ms2)+2*ms2*cDg**5*((trs+trg+
     &   2*trc)*gDs+(trs+trg+4*trc)*cDs+mc2*(trs+trg+trc))+2*mc2*ms2*cDg
     &   **2*gDs**2*((trs+trc)*gDs-trg*(cDs+ms2))+2*mc2*ms2*cDg**3*gDs*(
     &   (trs+trc)*gDs+ms2*(trs+trg+trc))+2*ms2*(trs+trg+2*trc)*cDg**6)/
     &   (mc*ms*cDg**3*gDs**2)/4.0_dp+Ampm

      Ampm = Ampm-lc*mc*(cDg*(2*(trsgc-mc2*trs+ms2*trc)*gDs+ms2*(mc2*t
     &   rg-trsgc))+mc2*(trsgc-ms2*trg)*gDs+2*ms2*trg*cDg**2)/(ms*cDg*(2
     &   *cDg+mc2))/2.0_dp

      Ampm=Ampm/zmpm

      return
      end
