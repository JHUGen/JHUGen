      subroutine Bamp_ppm(q,mc,ms,Bppm)
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
     & mc2,ms2,u,xsn,xsd,xs
      complex(dp):: trc,trg,trs,trsgc,zppm,Bppm

      mc2=mc**2
      ms2=ms**2

      cDs=dot(q,3,4)+mc2*dot(q,4,2)/2._dp/dot(q,3,2)
     & +ms2*dot(q,3,2)/2._dp/dot(q,4,2)
      cDg=dot(q,3,2)
      gDs=dot(q,4,2)
      u=mc2+ms2+2._dp*cDs

      xsn=(1._dp-sqrt(1._dp-4._dp*ms*mc/(u-(ms-mc)**2)))
      xsd=(1._dp+sqrt(1._dp-4._dp*ms*mc/(u-(ms-mc)**2)))
      xs=-xsn/xsd

      trg=2._dp*za(5,2)*zb(2,1)
      trs=2._dp*za(5,4)*zb(4,1)+ms**2*za(5,2)*zb(2,1)/gDs
      trc=2._dp*za(5,3)*zb(3,1)+mc**2*za(5,2)*zb(2,1)/cDg
      trsgc=2._dp*zb(1,4)*za(4,2)*zb(2,3)*za(3,5)
      zppm=za(2,3)**2*zb(4,3)

      Bppm = mc*ms*(gDs+cDg)*(cDs*(-4*trc*gDs**4+2*(trsgc+mc2*trg-2*ms2
     &   *trc)*gDs**3+ms2*(-2*(trs+trg)*cDg+trsgc+3*mc2*trg)*gDs**2+ms2*
     &   cDg*(2*(trs+trg)*cDg-trsgc-ms2*trg)*gDs+2*ms2**2*trg*cDg**2)+ms
     &   2*gDs*(cDg*(4*trc*gDs**2-2*(trsgc+mc2*trs-2*ms2*trc)*gDs+ms2*(m
     &   c2*trg-trsgc))+mc2*gDs*(2*trs*gDs-2*trg*gDs+trsgc-3*ms2*trg))+2
     &   *trg*cDs**2*gDs*(2*gDs*(gDs+ms2)-cDg*(2*gDs+3*ms2)))*tr5Xs/((cD
     &   s**2-mc2*ms2)*gDs**3)/2.0_dp

      Bppm = mc*mc2*ms*(gDs+cDg)*(2*mc2*trg*cDs*gDs**2+cDg**2*(cDs*(-2*
     &   trc*gDs+trsgc+3*ms2*trg)+4*mc2*trs*gDs-2*ms2*trc*gDs+4*trg*cDs*
     &   *2+ms2*trsgc-3*mc2*ms2*trg)-cDg*gDs*(cDs*(-2*trc*gDs+trsgc+mc2*
     &   trg)+6*trg*cDs**2+mc2*(trsgc-ms2*trg))+2*cDg**3*(ms2*trc-2*trs*
     &   cDs))*tr5Xc/(cDg**2*(cDs**2-mc2*ms2)*gDs)/2.0_dp+Bppm

      Bppm = mc*ms*(cDs*(4*trc*gDs**4-2*(trsgc+mc2*trg-2*ms2*trc)*gDs**
     &   3+ms2*(2*(trs+trg)*cDg-trsgc-3*mc2*trg)*gDs**2+ms2*cDg*(-2*(trs
     &   +trg)*cDg+trsgc+ms2*trg)*gDs-2*ms2**2*trg*cDg**2)+ms2*gDs*(cDg*
     &   (-4*trc*gDs**2+2*(trsgc+mc2*trs-2*ms2*trc)*gDs+ms2*(trsgc-mc2*t
     &   rg))+mc2*gDs*(2*(trg-trs)*gDs-trsgc+3*ms2*trg))+2*trg*cDs**2*gD
     &   s*(cDg*(2*gDs+3*ms2)-2*gDs*(gDs+ms2)))*tr4Xs/((cDs**2-mc2*ms2)*
     &   gDs**2)/2.0_dp+Bppm

      Bppm = mc*mc2*ms*(-2*mc2*trg*cDs*gDs**2-cDg**2*(cDs*(-2*trc*gDs+t
     &   rsgc+3*ms2*trg)+4*mc2*trs*gDs-2*ms2*trc*gDs+4*trg*cDs**2+ms2*tr
     &   sgc-3*mc2*ms2*trg)+cDg*gDs*(cDs*(-2*trc*gDs+trsgc+mc2*trg)+6*tr
     &   g*cDs**2+mc2*(trsgc-ms2*trg))+cDg**3*(4*trs*cDs-2*ms2*trc))*tr4
     &   Xc/(cDg*(cDs**2-mc2*ms2)*gDs)/2.0_dp+Bppm

      Bppm = Bppm-mc*ms*(2*mc2**2*ms2*trg*gDs**2+cDg**2*(cDs*(4*mc2*tr
     &   s*gDs-2*ms2*trc*gDs+ms2*trsgc+mc2*ms2*trg)+mc2*ms2*(2*trg*gDs-2
     &   *trc*gDs+trsgc+3*ms2*trg)-2*trg*cDs**2*gDs)+cDg*gDs*(2*(trg*cDs
     &   **2+mc2*ms2*(trc-trg))*gDs-mc2*((trsgc+5*ms2*trg)*cDs+ms2*(trsg
     &   c+mc2*trg)))+2*ms2*cDg**3*(trc*cDs-2*mc2*trs))*tr3Xs/(cDg*(cDs*
     &   *2-mc**2*ms**2)*gDs)/2.0_dp

      Bppm = Bppm-mc*ms*(cDg*gDs*(cDs*(4*trc*gDs**2-2*(trsgc+mc2*trs+2
     &   *mc2*trg-2*ms2*trc)*gDs-ms2*(trsgc+5*mc2*trg))-mc2*ms2*(2*(trs+
     &   trg)*gDs+trsgc+ms2*trg))+mc2*gDs**2*(-4*trc*gDs**2+2*((trs+trg)
     &   *cDs+trsgc+mc2*trg-2*ms2*trc)*gDs+trsgc*(cDs+ms2)+ms2*trg*(cDs+
     &   3*mc2))+2*mc2*ms2*cDg**2*((trs+trg)*gDs+ms2*trg))*tr3Xc/((cDs-m
     &   c*ms)*(cDs+mc*ms)*gDs**2)/2.0_dp

      Bppm = ms*(cDg**3*(cDs**2*(8*trc*gDs**3+2*(-trsgc+mc2*trs+3*ms2*t
     &   rc)*gDs**2-ms2*(trsgc+mc2*trg)*gDs-2*mc2*ms2**2*trc)+mc2*ms2*cD
     &   s*(2*(-trs-trg+trc)*gDs**2+ms2**2*trc)+mc2*ms2*gDs*(-4*trc*gDs*
     &   *2-4*(mc2*trs+ms2*trc)*gDs+ms2*(trsgc+mc2*trg))+2*ms2*trc*cDs**
     &   4-ms2**2*trc*cDs**3)+mc2*cDg**2*gDs**2*(2*cDs*gDs*(-trc*(2*gDs+
     &   ms2)+trsgc+mc2*trs)+mc2*ms2*((2*trs+4*trg-2*trc)*gDs+trsgc-ms2*
     &   trs)-cDs**2*(2*trg*gDs+trsgc-ms2*trs))+mc2*(cDs-mc*ms)*(cDs+mc*
     &   ms)*gDs**4*((trg+trc)*(8*gDs-ms2)+trs*(2*gDs+2*cDs+ms2-mc2))+2*
     &   mc2*cDg*(cDs-mc*ms)*(cDs+mc*ms)*gDs**3*(trg*gDs-ms2*trs)+2*ms2*
     &   trc*cDg**4*cDs*(cDs-mc*ms)*(cDs+mc*ms))*lVs/(mc*cDg**2*(cDs-mc*
     &   ms)*(cDs+mc*ms)*gDs**3)/4.0_dp+Bppm

      Bppm = Bppm-mc*(cDs**2*(2*(mc2*(trs+4*trg+4*trc)-trg*cDg)*gDs**5
     &   -(-2*trg*cDg**2+mc2**2*trs+mc2*ms2*(-trs+trg+trc))*gDs**4+cDg*(
     &   -4*trs*cDg**2-2*(-trsgc+2*mc2*trs+ms2*trc)*cDg+mc2*(trsgc-2*ms2
     &   *trs+ms2*trg))*gDs**3+ms2*cDg**2*(3*trg*cDg-3*trc*cDg+trsgc+mc2
     &   *trs+2*mc2*trg)*gDs**2-ms2**2*trc*cDg**3*gDs-2*mc2*ms2**2*trc*c
     &   Dg**3)+ms2*cDs*(2*mc2*trg*cDg*gDs**4-2*mc2**2*trs*gDs**4+cDg**3
     &   *(-4*trc*gDs**3+2*(trsgc+mc2*trs-ms2*trc)*gDs**2+mc2*ms2**2*trc
     &   )+2*mc2*(-trs-2*trg+trc)*cDg**2*gDs**3-2*mc2*ms2*trc*cDg**4)+cD
     &   s**3*(-2*trg*cDg*gDs**4+2*mc2*trs*gDs**4+2*trg*cDg**2*gDs**3+2*
     &   ms2*trc*cDg**4-ms2**2*trc*cDg**3)+mc2*ms2*gDs*(cDg**3*(4*trs*gD
     &   s**2+ms2*(2*trs-trg+trc)*gDs+ms2**2*trc)+cDg*gDs**2*(2*trg*gDs*
     &   *2-mc2*(trsgc+ms2*(trg-2*trs)))+cDg**2*gDs*((4*trc-2*trg)*gDs**
     &   2+2*(-2*trsgc+mc2*trs+2*ms2*trc)*gDs-ms2*(trsgc+mc2*(trs+2*trg)
     &   ))+mc2*gDs**3*(-2*(trs+4*trg+4*trc)*gDs+mc2*trs+ms2*(-trs+trg+t
     &   rc)))+2*ms2*trc*cDg**3*cDs**4)*lVc/(ms*cDg**2*(cDs**2-mc**2*ms*
     &   *2)*gDs**3)/4.0_dp

      Bppm = trg*xs*cDs*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)*lRcs/((xs
     &   **2-1)*cDg*gDs**2)+Bppm

      Bppm = (mc2*cDg**3*(-80*trs*gDs**4-cDs*(72*trs*gDs**3+4*ms2*(5*tr
     &   s+22*trg+30*trc)*gDs**2+2*ms2*(-17*mc2*trs+9*ms2*trg-26*mc2*trg
     &   +34*ms2*trc)*gDs+ms2**2*(-9*trsgc+mc2*(-8*trs-17*trg+3*trc)+5*m
     &   s2*trc))+4*ms2*(-9*trs-4*trg+37*trc)*gDs**3+ms2*(-66*trsgc+ms2*
     &   (88*trc-6*(trs+7*trg))+mc2*(158*trs-13*trg-3*trc))*gDs**2+9*(mc
     &   2-2*ms2)*ms2*trsgc*gDs+mc2*ms2**2*(26*trs+17*trg+13*trc)*gDs+6*
     &   (mc-ms)*(ms+mc)*ms2*trc*cDs**2+mc2*ms2**2*(9*trsgc+ms2*(-8*trs+
     &   trg+8*trc)))+mc2*cDg**2*gDs*(4*(41*trg+32*trc)*gDs**4-12*(-6*tr
     &   g*cDs-3*trsgc+7*mc2*trs+3*ms2*(trc-2*trg))*gDs**3+2*(18*trg*cDs
     &   **2+2*(9*trsgc-27*mc2*trs+26*ms2*trg+17*ms2*trc)*cDs+ms2*(9*trs
     &   gc+6*ms2*trs-31*mc2*trs+8*mc2*trg+58*mc2*trc))*gDs**2-ms2*(2*(9
     &   *trsgc+5*mc2*trs+31*mc2*trg+13*mc2*trc)*cDs+mc2*(54*trsgc+3*ms2
     &   *trs-31*mc2*trs+39*ms2*trg+23*mc2*trg-44*ms2*trc))*gDs-9*mc2*ms
     &   2*((trsgc+ms2*trg)*cDs+ms2*(trsgc+mc2*trg)))+mc2*cDg*gDs**2*(4*
     &   (23*trg+32*trc)*gDs**4+cDs*(8*(7*trg+16*trc)*gDs**3+12*(-ms2*tr
     &   s+mc2*trs+3*mc2*trg)*gDs**2+2*mc2*(18*trsgc-18*mc2*trs+35*ms2*t
     &   rg+17*ms2*trc)*gDs-9*mc2*ms2*(trsgc+mc2*trg))-2*(6*ms2*(trs+4*t
     &   rg+4*trc)+mc2*(-6*trs-trg+8*trc))*gDs**3+(4*mc2*(9*trsgc+3*ms2*
     &   trs+8*ms2*trg-19*ms2*trc)-12*mc2**2*trs+6*ms2**2*(-trs+trg+trc)
     &   )*gDs**2+18*trg*cDs**2*(mc2-2*gDs)*gDs+mc2*ms2*(27*trsgc+6*ms2*
     &   trs+mc2*(-22*trs+27*trg+18*trc))*gDs-9*mc2**2*ms2*(trsgc+ms2*tr
     &   g))+mc2**2*gDs**3*((46*trg+64*trc)*gDs**3+cDs*(4*(7*trg+16*trc)
     &   *gDs**2+6*(mc-ms)*(ms+mc)*trs*gDs+9*mc2*(trsgc+ms2*trg))+2*(mc2
     &   *(3*trs-20*(trg+trc))-3*ms2*(trs+4*(trg+trc)))*gDs**2-18*trg*cD
     &   s**2*gDs+(mc2*(9*trsgc+ms2*(6*trs-2*trg-29*trc))+5*mc2**2*trs+3
     &   *ms2**2*(-trs+trg+trc))*gDs+9*mc2*ms2*(trsgc+mc2*trg))+ms2*cDg*
     &   *4*(mc2*(-24*(-8*trs-3*trg+3*trc)*gDs**2+(36*trsgc+ms2*(52*trs+
     &   70*trg-28*trc)-18*mc2*trs)*gDs+ms2*(27*trsgc-9*mc2*(2*trs+trg)+
     &   2*ms2*(-8*trs+trg+8*trc)))-2*cDs*(68*trc*gDs**2+2*(-9*trsgc-17*
     &   mc2*trs-26*mc2*trg+34*ms2*trc)*gDs+ms2*(5*ms2*trc-9*trsgc)+mc2*
     &   ms2*(-8*trs-17*trg+6*trc)-3*mc2**2*trc)+12*(mc-ms)*(ms+mc)*trc*
     &   cDs**2)-6*ms2*cDg**5*(24*trc*gDs**2+6*(-trsgc+mc2*trs+3*ms2*trc
     &   )*gDs+2*(ms**2-mc**2)*trc*cDs+3*ms2*(-trsgc+2*mc2*trs+mc2*trg))
     &   )/(mc*ms*cDg**2*(2*cDg+mc2)*gDs**3)/1.2e+1_dp+Bppm

      Bppm = tr3s002ft*(3*mc*trs*(2*cDg+mc2)**2*gDs**2-3*mc*(trg+trc)*(
     &   8*(2*cDg+mc2)*gDs**3+ms2*cDg**2*(gDs+cDs)))/(ms*cDg**2*gDs)+Bp
     &   pm

      Bppm = B0cgsf*(cDg**3*(-4*(mc2*trs+2*ms2*trc)*gDs**3+ms2*cDs*(-8*
     &   trc*gDs**2+2*(trsgc+2*mc2*trs+3*mc2*trg-4*ms2*trc)*gDs+ms2*(trs
     &   gc+mc2*(trs+2*trg-trc)))+ms2*(2*trsgc+mc2*(5*trs+7*trg-3*trc)-6
     &   *ms2*trc)*gDs**2+ms2*(ms2+mc2)*trsgc*gDs+mc2*ms2**2*(3*trs+5*tr
     &   g+trc)*gDs+2*(mc-ms)*(ms+mc)*ms2*trc*cDs**2+mc2*ms2**2*(trsgc+m
     &   s2*(trc-trs)))+mc2*cDg**2*gDs*(2*(trg-2*trs)*gDs**3+cDs*(2*(trg
     &   -2*trs)*gDs**2-ms2*(trs+trg+3*trc)*gDs-ms2*(trsgc+ms2*trg))+(2*
     &   trsgc-3*ms2*trs-4*mc2*trs+5*ms2*trg+3*ms2*trc)*gDs**2+ms2*(-2*t
     &   rsgc-ms2*trs+4*mc2*trs+mc2*trg+5*ms2*trc)*gDs-ms2**2*(trsgc+mc2
     &   *trg))+mc2*gDs**3*((6*trg+8*trc)*gDs**3+cDs*(4*(trg+2*trc)*gDs*
     &   *2+2*(mc-ms)*(ms+mc)*trs*gDs+mc2*(trsgc+ms2*trg))+2*(mc2*trs-ms
     &   2*(trs+4*(trg+trc)))*gDs**2-2*trg*cDs**2*gDs+mc2*trsgc*gDs+ms2*
     &   (ms2*(-trs+trg+trc)-mc2*(-2*trs+trg+4*trc))*gDs+mc2*ms2*(trsgc+
     &   mc2*trg))+mc2*cDg*gDs**2*(8*(trg+trc)*gDs**3+2*(trg*cDs+trsgc-m
     &   c2*trs+ms2*(trg-trc))*gDs**2+2*(trg*cDs**2+(-2*mc2*trs+3*ms2*tr
     &   g+2*ms2*trc)*cDs+ms2*(ms2*trs+mc2*(-2*trs+trg+trc)))*gDs+trsgc*
     &   (2*cDs+ms2+mc2)*gDs+ms2*(-mc2*trg*(cDs+ms2)-trsgc*(cDs+mc2)))-m
     &   s2*cDg**4*(8*trc*gDs**2+2*(-trsgc+mc2*trs+3*ms2*trc)*gDs+2*(ms*
     &   *2-mc**2)*trc*cDs+ms2*(-trsgc+2*mc2*trs+mc2*trg)))/(mc*ms*cDg**
     &   2*gDs**3)/4.0_dp+Bppm

      Bppm = Bppm-mc*tr3s00ft*(cDg*((26*trg+24*trc)*gDs**3+2*(2*trg*cD
     &   s+trsgc-3*mc2*trs+2*ms2*trg-ms2*trc)*gDs**2+(2*trg*cDs**2+2*(tr
     &   sgc-2*mc2*trs+5*ms2*trg+4*ms2*trc)*cDs+ms2*(trsgc+2*mc2*(-3*trs
     &   +trg+trc)))*gDs-ms2*((trsgc+mc2*trg)*cDs+mc2*(trsgc+ms2*trg)))+
     &   gDs*((22*trg+24*trc)*gDs**3+cDs*(4*(5*trg+6*trc)*gDs**2+mc2*(tr
     &   sgc+ms2*trg))-2*trg*cDs**2*gDs+mc2*(trsgc-2*ms2*trg-5*ms2*trc)*
     &   gDs+mc2*ms2*(trsgc+mc2*trg))-cDg**2*(8*trs*gDs**2+cDs*(4*trs*gD
     &   s+ms2*(-3*trs+4*trg+2*trc))+ms2*(9*trs-4*trc)*gDs+ms2*(3*trsgc-
     &   4*mc2*trs+2*ms2*trg))+9*ms2*trs*cDg**3)/(ms*cDg**2*gDs)

      Bppm = lc*mc*(cDg*(-4*trg*gDs**3+2*trg*(mc2-2*cDs)*gDs**2+2*mc2*(
     &   trg*cDs+2*trsgc-2*mc2*trs+ms2*(trc-2*trg))*gDs+mc2*ms2*(3*mc2*t
     &   rg-trsgc))+2*cDg**2*(2*trg*gDs**2+2*(trg*cDs+trsgc-3*mc2*trs+2*
     &   ms2*trc)*gDs-ms2*(2*trsgc+mc2*(trc-5*trg)))-mc2*trg*gDs*(2*gDs*
     &   (gDs+cDs)+mc2*ms2)+cDg**3*(6*ms2*(trg-trc)-8*trs*gDs)+mc2**2*tr
     &   sgc*gDs)/(ms*(2*cDg+mc2)**2*gDs)/2.0_dp+Bppm

      Bppm = LsB2*mc*ms*(-2*mc2**2*trg*cDs*gDs**3+cDg**3*(cDs*(8*mc2*tr
     &   s*gDs-2*ms2*trc*gDs+ms2*trsgc-7*mc2*ms2*trg)+4*cDs**2*(-trc*gDs
     &   +trsgc+ms2*trg)+8*trg*cDs**3-mc2*ms2*(3*trsgc+ms2*trg))-2*cDg**
     &   2*gDs*(2*mc2**2*(trs*gDs-ms2*trg)+mc2*cDs*(trsgc-trc*gDs)+cDs**
     &   2*(-2*trc*gDs+trsgc+3*mc2*trg)+4*trg*cDs**3)+mc2*cDg*gDs**2*(cD
     &   s*(-2*trc*gDs+trsgc+mc2*trg)+8*trg*cDs**2+mc2*(trsgc-ms2*trg))+
     &   cDg**4*(-8*trs*cDs**2+2*ms2*trc*cDs+4*mc2*ms2*trs))/(cDg**2*(cD
     &   s**2-mc2*ms2)*gDs)+Bppm

      Bppm = ms*tr3c00fs*(cDg*(cDs*(16*trc*gDs**2+2*(-trsgc-mc2*(trs+2*
     &   trg)+8*ms2*trc)*gDs+ms2*(-trsgc-mc2*trg+3*ms2*trc))-mc2*(2*(trs
     &   +2*trg)*gDs**2+(trsgc+ms2*(trg-trs))*gDs+ms2*(trsgc-3*ms2*trs-2
     &   *ms2*trg)))+cDg**2*(8*trc*gDs**2+2*(-trsgc+mc2*trs+3*ms2*trc)*g
     &   Ds+ms2*(mc2*(2*trs+trg)-trsgc))+mc2*gDs*(-10*trc*gDs**2+(2*(trs
     &   +trg)*cDs+trsgc+2*mc2*trg-7*ms2*trc)*gDs+trsgc*(cDs+ms2)+ms2*tr
     &   g*(cDs+mc2)))/(mc*gDs**3)+Bppm

      Bppm = BfunX*(ms2*cDg**3*(mc2*((trs+trg-trc)*gDs**2+ms2*(trs+trg+
     &   trc)*gDs+ms2**2*(trs+trg+trc))+mc2*cDs*(2*(trs+trg)*gDs+ms2*(tr
     &   s+trg-trc))+2*trc*cDs**2*(2*gDs+mc2)+4*trc*cDs**3)+2*mc2*cDg*gD
     &   s**3*((trg+trc)*(4*gDs**2-ms2*(gDs+cDs))+trs*(gDs-ms2)*(2*(gDs+
     &   cDs)+ms2))+mc2*gDs**4*((trg+trc)*(8*gDs*(gDs+cDs)-2*ms2*(cDs-3*
     &   gDs)-ms2**2)+trs*(2*(gDs+cDs)+ms2)**2)+mc2*ms2*cDg**2*gDs**2*(t
     &   rs*(-gDs+cDs+ms2)-(trc-trg)*(gDs+cDs+ms2))+2*ms2*trc*cDg**4*cDs
     &   *(2*gDs+4*cDs+mc2)+4*ms2*trc*cDg**5*cDs)/(mc*ms*cDg**2*gDs**3)/
     &   4.0_dp+Bppm

      Bppm = LsB1*mc*ms*(-2*cDg*gDs**2*(cDs*(-4*trc*gDs**2+(2*trsgc+mc2
     &   *trs+3*mc2*trg-4*ms2*trc)*gDs+ms2*(trsgc+2*mc2*trg))+mc2*ms2*(2
     &   *trs*gDs+trsgc-ms2*trg)+2*trg*cDs**2*(gDs+ms2)+2*trg*cDs**3)+cD
     &   g**2*gDs*(-4*ms2*trc*gDs**2+2*(2*trg*cDs**2+ms2*((trs+trg)*cDs+
     &   trsgc+2*mc2*trs+mc2*trg)-2*ms2**2*trc)*gDs+ms2*(trg*(8*cDs**2+m
     &   s2*cDs-mc2*ms2)+trsgc*(cDs+ms2)))+gDs**3*(-4*mc2*trc*gDs**2+2*m
     &   c2*((trs+trg)*cDs+trsgc+mc2*trg-2*ms2*trc)*gDs+2*(mc2*trg-trsgc
     &   )*cDs**2+mc2*(trsgc+ms2*trg)*cDs+mc2*ms2*(3*trsgc+mc2*trg))-2*m
     &   s2*cDg**3*cDs*((trs+trg)*gDs+ms2*trg))/((cDs-mc*ms)*(cDs+mc*ms)
     &   *gDs**3)+Bppm

      Bppm = B0csf*mc*ms*(mc2*(2*trc*gDs**2+(-trsgc-mc2*trs+ms2*(trg-tr
     &   s))*gDs+ms2*(ms2+mc2)*trg)+cDs*(2*trc*gDs**2+(-trsgc+mc2*(-2*tr
     &   s-trg+trc)+ms2*trc)*gDs+2*mc2*ms2*trg)+cDg*(ms2*(-trc*(2*gDs+ms
     &   2+mc2)+trsgc+mc2*trg)+cDs*(-2*trc*gDs+trsgc+ms2*(trs+trg-2*trc)
     &   +mc2*trs)+2*trs*cDs**2)+cDs**2*(2*trc*gDs-trg*(2*gDs+ms2+mc2))-
     &   2*trg*cDs**3)/((cDs-mc*ms)*(cDs+mc*ms)*gDs)/2.0_dp+Bppm

      Bppm = ls*ms*(cDg*(-8*trc*gDs**2+2*(trsgc+mc2*trs-3*ms2*trc)*gDs+
     &   ms2*(trsgc-mc2*trg))+mc2*gDs*(-2*trs*gDs+4*trg*gDs-trsgc+3*ms2*
     &   trg))/(mc*gDs*(2*gDs+ms2))/2.0_dp+Bppm

      Bppm = 3*mc*ms*tr3c002fs*(cDg*((trs+trg)*cDs+ms2*trc)*(2*gDs+ms2)
     &   -mc2*(trs+trg)*gDs**2)/gDs**3+Bppm

      Bppm = 3*ms*tr3c001fs*(2*gDs+ms2)*(cDg*(trc*cDs*(2*gDs+ms2)+mc2*m
     &   s2*(trs+trg))-mc2*trc*gDs**2)/(mc*gDs**3)+Bppm

      Bppm = epinv*trg*(4*lp*xs*cDs+mc*(ms-ms*xs**2))*(mc2*gDs**2-2*cDg
     &   *cDs*gDs+ms2*cDg**2)/((xs**2-1)*cDg*gDs**2)/4.0_dp+Bppm

      Bppm = mc*ms*tr2fu*trg*cDs*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)/
     &   (cDg*gDs**2)+Bppm

      Bppm = 3*mc*tr3s001ft*((trg+trc)*gDs*(ms2*(mc2*gDs-2*cDg*cDs)-8*g
     &   Ds**2*(gDs+cDs))+ms2*trs*cDg*(3*cDg*gDs+2*mc2*gDs-cDg*cDs))/(ms
     &   *cDg**2*gDs)+Bppm

      Bppm=Bppm/zppm

      return
      end
