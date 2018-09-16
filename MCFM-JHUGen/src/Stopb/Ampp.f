      subroutine Aamp_mpp(q,mc,ms,Ampp)
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
      complex(dp):: trc,trg,trs,trsgc,zmpp,Ampp

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
      zmpp=za(4,3)**2*zb(2,4)*zb(2,3)

      Ampp = mc*ms2*(-mc2*trsgc*gDs**3-cDg**2*gDs*(4*trs*cDs*gDs-2*ms2*
     &   trg*cDs+ms2*trsgc)+cDg*gDs**2*((trsgc+2*mc2*trs+mc2*trg)*gDs+2*
     &   trsgc*cDs-mc2*ms2*trg)+ms2*cDg**3*(2*trs*gDs-ms2*trg))*tr3Xs/(c
     &   Dg*gDs**3)

      Ampp = mc*(-cDg**3*(4*trc*gDs**3+2*(trsgc+mc2*trg)*gDs**2+2*ms2*(
     &   trsgc-mc2*trs)*gDs+mc2*ms2**2*trg)-mc2**2*trsgc*gDs**3+cDg**2*g
     &   Ds*(mc2*(2*trg*gDs**2-ms2*trsgc)+2*cDs*((trsgc-mc2*(2*trs+trg))
     &   *gDs+mc2*ms2*trg))+mc2*cDg*gDs**2*((mc2*(2*trs+trg)-trsgc)*gDs+
     &   2*trsgc*cDs-mc2*ms2*trg)+2*trg*cDg**5*gDs)*tr3Xc/(cDg**3*gDs)+
     &   Ampp

      Ampp = mc*(2*trs*cDg*gDs-trsgc*gDs-ms2*trg*cDg)*(mc2*gDs**2-2*cDg
     &   *cDs*gDs+ms2*cDg**2)*tr1Xfs/(cDg**2*gDs)+Ampp

      Ampp = mc*(2*trs*cDg*gDs-trsgc*gDs-ms2*trg*cDg)*(mc2*gDs**2-2*cDg
     &   *cDs*gDs+ms2*cDg**2)*tr1Xfc/(cDg*gDs**2)+Ampp

      Ampp = (-cDg**5*(4*(trs+trg)*gDs**3+cDs*(4*trg*gDs**2-2*ms2**2*tr
     &   c)-4*ms2*trc*gDs**2+2*ms2**2*trc*gDs+4*ms2*trc*cDs**2+ms2**3*tr
     &   c)+cDg**4*gDs*(16*trc*gDs**3+4*(trsgc+mc2*trg+2*ms2*trc)*gDs**2
     &   +4*ms2*(trsgc-mc2*trs)*gDs+ms2**2*(trsgc+mc2*trg))+cDg**6*(2*ms
     &   2*(2*(trs+trg)*gDs-4*trc*cDs+ms2*trc)-4*trg*gDs**2)+mc2*cDg**3*
     &   gDs**2*((2*gDs+ms2)*(2*trs*gDs-2*ms2*trs-3*ms2*trg)+8*(trs+trg)
     &   *cDs*gDs+2*ms2*trsgc)+mc2*cDg*gDs**4*(ms2*(2*(-2*trs+trg+trc)*g
     &   Ds+3*trs*(mc2-2*cDs)+3*ms2*(-trs+trg+trc))+2*mc2*trsgc)+mc2*cDg
     &   **2*gDs**3*(2*(2*trsgc-mc2*(2*trs+trg)+ms2*(trs+trg+trc))*gDs-4
     &   *trsgc*cDs-ms2*(trsgc+mc2*trg)+2*ms2**2*(trs+trg+trc))-mc2*ms2*
     &   (trs+trg+trc)*gDs**5*(2*gDs+2*cDs+ms2-mc2)-4*ms2*trc*cDg**7)*lV
     &   s/(mc*cDg**3*gDs**3)/4.0_dp+Ampp

      Ampp = mc*(cDg**3*gDs**2*(4*(trg+2*trc)*gDs**2-4*(trg+trc)*cDs*gD
     &   s+2*(3*trsgc+3*mc2*trs+5*mc2*trg)*gDs+ms2*(3*trsgc+2*mc2*trs+4*
     &   mc2*trg))+cDg**5*(-4*(trg+trc)*gDs**2+2*ms2*trc*gDs+4*trc*cDs**
     &   2-2*ms2*trc*cDs+ms2**2*trc)+2*cDg**2*gDs**3*(mc2*((-trs-trg+trc
     &   )*gDs-ms2*(trs+trg+trc))-trsgc*(gDs+2*cDs))+cDg**6*(-4*(trs+trg
     &   )*gDs+8*trc*cDs-2*ms2*trc)+2*cDg**4*gDs**2*(2*(trs+trc)*gDs+ms2
     &   *(trg+2*trc))+mc2*cDg*gDs**4*(-2*(-2*trs+trg+trc)*gDs+6*trs*cDs
     &   +trsgc-mc2*(3*trs+trg)-3*ms2*(-trs+trg+trc))+mc2*(trs+trg+trc)*
     &   gDs**5*(2*gDs+2*cDs+ms2-mc2)+4*trc*cDg**7)*lVc/(cDg**3*gDs**3)/
     &   4.0_dp+Ampp

      Ampp = mc*(2*trs*cDg*gDs-trsgc*gDs-ms2*trg*cDg)*(mc2*gDs**2-2*cDg
     &   *cDs*gDs+ms2*cDg**2)*lRs2/(cDg**2*gDs**2)/2.0_dp+Ampp

      Ampp = epinv*mc*(2*trs*cDg*gDs-trsgc*gDs-ms2*trg*cDg)*(mc2*gDs**2
     &   -2*cDg*cDs*gDs+ms2*cDg**2)*(2*lRs1+2*lRc1-1)/(cDg**2*gDs**2)/4.
     &   0_dp+Ampp

      Ampp = mc*(2*trs*cDg*gDs-trsgc*gDs-ms2*trg*cDg)*(mc2*gDs**2-2*cDg
     &   *cDs*gDs+ms2*cDg**2)*lRc2/(cDg**2*gDs**2)/2.0_dp+Ampp

      Ampp = 2*LsA*mc*(mc2**2*trsgc*gDs**5+cDg**3*gDs**2*(4*trc*gDs**3-
     &   2*trsgc*gDs**2-(4*trs*cDs**2-ms2*trsgc+mc2*ms2*(4*trs+trg))*gDs
     &   +ms2*(2*trg*(cDs**2+mc2*ms2)-3*trsgc*cDs))+cDg**2*gDs**3*(mc2*(
     &   cDs*(2*(3*trs+trg)*gDs-3*ms2*trg)-2*trg*gDs**2)+2*trsgc*(-cDs*g
     &   Ds+cDs**2+mc2*ms2))+ms2*cDg**4*gDs*(6*trs*cDs*gDs-3*ms2*trg*cDs
     &   +ms2*trsgc)+mc2*cDg*gDs**4*((trsgc-mc2*(2*trs+trg))*gDs-3*trsgc
     &   *cDs+mc2*ms2*trg)+ms2**2*cDg**5*(ms2*trg-2*trs*gDs))/(cDg**3*gD
     &   s**3)+Ampp

      Ampp = B0cgsf*(-cDg**4*(16*trc*gDs**4+4*(trsgc+mc2*(trs-2*trc)+2*
     &   ms2*trc)*gDs**3+cDs*(2*gDs+ms2)*(8*trc*gDs**2+2*(trsgc+2*mc2*tr
     &   s+3*mc2*trg)*gDs+ms2*(trsgc+mc2*trg))+2*ms2*(2*trsgc+mc2*(trg-4
     &   *trc))*gDs**2+ms2*(mc2*(2*trsgc+ms2*(trs+4*trg-3*trc))+ms2*trsg
     &   c)*gDs+mc2*ms2**2*(trsgc+ms2*trg))-cDg**5*(4*(-trs-trg+5*trc)*g
     &   Ds**3+2*cDs*(2*(-trs-2*trg+trc)*gDs**2+ms2*(ms2-mc2)*trc)+4*(tr
     &   sgc+2*mc2*trg+3*ms2*trc-mc2*trc)*gDs**2+4*cDs**2*((mc-ms)*(ms+m
     &   c)*trc-trg*gDs)+2*ms2*(2*trsgc+mc2*(trc-4*trs)+ms2*trc)*gDs+ms2
     &   **2*(trsgc+mc2*(-3*trs-2*trg+trc)))+2*mc2*cDg**3*gDs*(2*(-trs-t
     &   rg+trc)*gDs**3+cDs*(2*(trc-2*(2*trs+trg))*gDs**2+2*trsgc*gDs+ms
     &   2*(trg-2*trs)*gDs+ms2*(trsgc+ms2*trg))+(trsgc+ms2*(2*trs+3*trg+
     &   2*trc)+2*mc2*trs+3*mc2*trg)*gDs**2+ms2*(3*trsgc+ms2*(trs+3*trg-
     &   trc)+4*mc2*trg)*gDs+ms2**2*(trsgc+mc2*trg))-mc2*cDg*gDs**3*(2*(
     &   -trsgc+mc2*(2*trs-trg+trc)+ms2*(-2*trs+trg+trc))*gDs**2-(mc2*(-
     &   2*trs*cDs+4*trg*cDs+trsgc+ms2*(-6*trs+8*trg+6*trc))+2*trsgc*cDs
     &   +6*ms2*trs*cDs-3*ms2**2*(-trs+trg+trc)+mc2**2*(trg-2*trs))*gDs+
     &   2*mc2*(trsgc*(cDs+ms2)+ms2*trg*(cDs+mc2)))-mc2*gDs**4*(2*(mc-ms
     &   )*(ms+mc)*(trs+trg+trc)*gDs**2+cDs*(2*(mc-ms)*(ms+mc)*(trs+trg+
     &   trc)*gDs-mc2*(trsgc+mc2*trg))-(mc2*(trsgc-2*ms2*(trs+trg+trc))+
     &   ms2**2*(trs+trg+trc)+mc2**2*trg)*gDs-mc2**2*(trsgc+ms2*trg))-2*
     &   mc2*cDg**2*gDs**3*(2*(-trs+trg+2*trc)*gDs**2+cDs*(2*(trc-2*trs)
     &   *gDs+trsgc+mc2*(trg-2*trs)+2*ms2*trg)-(mc2*(trs+trg)+ms2*(2*(tr
     &   g+trc)-trs))*gDs+(ms2+mc2)*trsgc+ms2*(ms2*(trs+trg+trc)+mc2*(4*
     &   trg-2*trs)))-2*cDg**6*(2*(trc-trg)*gDs**2-2*((trg-trs)*cDs+mc2*
     &   (trs+trg))*gDs+(mc-ms)*(ms+mc)*trc*(4*cDs-ms2))-4*cDg**7*(trs*g
     &   Ds+(mc-ms)*(ms+mc)*trc))/(mc*cDg**3*gDs**3)/4.0_dp+Ampp

      Ampp = (cDg**5*(-336*trc*gDs**4+4*(-18*trsgc+ms2*(3*(trs+trg)-80*
     &   trc)+6*mc2*(trc-9*trg))*gDs**3+2*cDs*(2*gDs+ms2)*(3*(mc-ms)*(ms
     &   +mc)*ms2*trc-16*(-trs-trg+trc)*gDs**2)-4*ms2*(27*trsgc-34*mc2*t
     &   rs+29*mc2*trg+31*ms2*trc)*gDs**2+12*cDs**2*(2*gDs+ms2)*(3*trg*g
     &   Ds+ms2*trc-mc2*trc)-2*ms2**2*(27*trsgc+3*mc2*(2*trc-5*(4*trs+tr
     &   g))+17*ms2*trc)*gDs+ms2**3*(-9*trsgc+mc2*(26*trs+17*trg-3*trc)-
     &   5*ms2*trc))-cDg**4*(2*gDs+ms2)*(mc2*(-12*(3*trg+7*trc)*gDs**3+m
     &   s2*(34*trs-84*trc)*gDs**2+6*ms2*(3*trsgc+ms2*(trs+4*trg-4*trc))
     &   *gDs+9*ms2**2*(trsgc+ms2*trg))+9*cDs*(2*gDs+ms2)*(8*trc*gDs**2+
     &   2*(trsgc+2*mc2*trs+3*mc2*trg)*gDs+ms2*(trsgc+mc2*trg)))+mc2*cDg
     &   **3*gDs*(2*gDs+ms2)*(108*trc*gDs**3+6*cDs*(2*(-5*trs-2*trg+trc)
     &   *gDs**2+6*trsgc*gDs+3*ms2*(trg-2*trs)*gDs+3*ms2*(trsgc+ms2*trg)
     &   )+6*(12*trsgc+ms2*(3*trs+2*trg+7*trc)+11*mc2*trs+20*mc2*trg)*gD
     &   s**2+ms2*(99*trsgc+ms2*(6*trs+35*trg-16*trc)+10*mc2*trs+96*mc2*
     &   trg)*gDs+18*ms2**2*(trsgc+mc2*trg))-mc2*cDg*gDs**3*(2*gDs+ms2)*
     &   (2*(-9*trsgc+mc2*(8*trs-5*trg+13*trc)+3*ms2*(-2*trs+trg+trc))*g
     &   Ds**2+3*(-mc2*(6*(trs*cDs+2*trg*cDs+2*trsgc)+ms2*(-6*trs+17*trg
     &   +11*trc))+3*(-2*trsgc*cDs-2*ms2*trs*cDs+ms2**2*(-trs+trg+trc))+
     &   11*mc2**2*trs)*gDs+18*mc2*(trsgc*(cDs+ms2)+ms2*trg*(cDs+mc2)))+
     &   mc2*gDs**4*(2*gDs+ms2)*(-6*(mc-ms)*(ms+mc)*(trs+trg+trc)*gDs**2
     &   +cDs*(9*mc2*(trsgc+mc2*trg)-6*(mc-ms)*(ms+mc)*(trs+trg+trc)*gDs
     &   )+(3*mc2*(3*trsgc-2*ms2*(trs+trg+trc))+3*ms2**2*(trs+trg+trc)+m
     &   c2**2*(-7*trs+4*trg-5*trc))*gDs+9*mc2**2*(trsgc+ms2*trg))-mc2*c
     &   Dg**2*gDs**3*(2*gDs+ms2)*(4*(-8*trs+8*trg+17*trc)*gDs**2+18*cDs
     &   *(2*(trc-2*trs)*gDs+5*trsgc+mc2*(trg-2*trs)+2*ms2*trg)-2*(9*trs
     &   gc+3*ms2*(7*(trg+trc)-trs)+mc2*(-14*trs-11*trg+5*trc))*gDs+9*(3
     &   *ms2+2*mc2)*trsgc+ms2*(6*ms2*(trs+trg+trc)+mc2*(-22*trs+89*trg+
     &   10*trc)))-2*cDg**6*(2*gDs+ms2)*(16*trc*gDs**2-2*(2*(5*trg-4*trs
     &   )*cDs+(5*ms2+3*mc2)*(trs+trg))*gDs-3*(mc-ms)*(ms+mc)*trc*(ms2-4
     &   *cDs))+4*cDg**7*(2*gDs+ms2)*((trg-8*trs)*gDs+3*ms2*trc-3*mc2*tr
     &   c))/(mc*cDg**3*gDs**3*(2*gDs+ms2))/1.2e+1_dp+Ampp

      Ampp = tr3c00fs*(cDg**2*(mc2*(4*(6*trs+5*trg)*gDs**3+24*ms2*(trs+
     &   trg)*gDs**2+ms2*(2*trsgc+7*ms2*trs+10*ms2*trg)*gDs+ms2**2*(trsg
     &   c+ms2*trg))+cDs*(2*gDs+ms2)*(8*trc*gDs**2+2*(trsgc+2*mc2*trs+3*
     &   mc2*trg)*gDs+ms2*(trsgc+mc2*trg)))+cDg**3*(20*trc*gDs**3+4*(-3*
     &   (trs+trg)*cDs+trsgc+3*mc2*trg+6*ms2*trc)*gDs**2+4*(-trg*cDs**2+
     &   ms2*(trsgc+mc2*(trg-trs))+4*ms2**2*trc)*gDs+ms2**2*(trsgc-mc2*(
     &   2*trs+trg)+3*ms2*trc))+mc2*cDg*gDs*(12*trc*gDs**3+2*(2*(3*trs+t
     &   rg)*cDs-4*trsgc-2*mc2*trs-5*mc2*trg+12*ms2*trc)*gDs**2+ms2*(4*t
     &   rs*cDs-4*trg*cDs-2*mc2*trs-9*mc2*trg+9*ms2*trc)*gDs-trsgc*(4*cD
     &   s+9*ms2)*gDs-2*ms2*trsgc*(cDs+ms2)-2*ms2**2*trg*(cDs+mc2))+mc2*
     &   gDs**2*(mc2*(trg*(ms2*(6*gDs+cDs)+2*gDs*(3*gDs+cDs)+ms2**2)+ms2
     &   *trs*gDs)+trsgc*(2*(2*cDs+mc2)*gDs+ms2*(cDs+mc2)))-4*cDg**4*gDs
     &   *(3*(trs+trg)*(2*gDs+ms2)+2*trg*cDs)-4*trg*cDg**5*gDs)/(mc*cDg*
     &   gDs**3)+Ampp

      Ampp = BfunX*(-cDg**5*(4*(trs+trg-trc)*gDs**3+cDs*(4*(trs+trg-trc
     &   )*gDs**2-2*mc2*ms2*trc)-4*mc2*trc*gDs**2+4*trc*cDs**2*(2*gDs+mc
     &   2)+2*mc2*ms2*(trc-2*(trs+trg))*gDs+8*trc*cDs**3+mc2*ms2**2*(-tr
     &   s-trg+trc))-2*mc2*cDg**3*gDs**2*((4*(trs+trg)-2*trc)*gDs**2+cDs
     &   *(4*trs*gDs+4*trg*gDs-2*trc*gDs+2*ms2*trs+3*ms2*trg)+ms2*(4*trs
     &   +5*trg-4*trc)*gDs+ms2**2*(trs+2*trg-trc))+mc2*cDg**4*gDs*((8*tr
     &   c-4*(trs+trg))*gDs**2+2*ms2*(2*trs+trg+4*trc)*gDs+3*ms2**2*(trs
     &   +trg+trc))-2*cDg**6*(-2*trc*gDs**2+cDs*((8*trc-2*(trs+trg))*gDs
     &   +4*mc2*trc)-2*mc2*(trs+trg)*gDs+12*trc*cDs**2-mc2*ms2*trc)-3*mc
     &   2*cDg*gDs**4*(2*(gDs+cDs)+ms2)*(2*trs*(gDs+cDs)-ms2*(-trs+trg+t
     &   rc))-mc2*(trs+trg+trc)*gDs**5*(2*(gDs+cDs)+ms2)**2+2*mc2*cDg**2
     &   *gDs**3*((-2*trs+trg+trc)*gDs+ms2*(trs+trg+trc))*(2*(gDs+cDs)+m
     &   s2)-4*cDg**7*((-trs-trg+2*trc)*gDs+6*trc*cDs+mc2*trc)-8*trc*cDg
     &   **8)/(mc*cDg**3*gDs**3)/4.0_dp+Ampp

      Ampp = ls*(4*cDg**3*(2*(trs+trg)*gDs**3+trg*cDs*gDs*(2*gDs+ms2))-
     &   cDg**2*(2*gDs+ms2)**2*(8*trc*gDs**2+2*(trsgc+mc2*trg)*gDs+ms2*(
     &   trsgc-mc2*trg))-2*mc2*cDg*gDs*(2*gDs+ms2)*(gDs*(2*trs*gDs-ms2*t
     &   rg)+2*trg*cDs*(2*gDs+ms2))+mc2*gDs**2*(2*gDs+ms2)*(mc2*trg*(2*g
     &   Ds+ms2)-trsgc*(4*gDs+ms2))+4*trg*cDg**4*gDs*(2*gDs+ms2))/(mc*cD
     &   g*gDs*(2*gDs+ms2)**2)/2.0_dp+Ampp

      Ampp = Ampp-mc*tr3s00ft*(cDg**2*(-4*(trc-3*trs)*gDs**3+cDs*(-4*(
     &   trc-2*trs)*gDs**2-4*(trsgc-mc2*trs+ms2*trg)*gDs+ms2*(trsgc+mc2*
     &   trg))+2*(trsgc+3*mc2*trs+9*ms2*(trg+trc))*gDs**2+ms2*(-3*trsgc+
     &   7*mc2*trs-5*mc2*trg)*gDs+mc2*ms2*(trsgc+ms2*trg))+cDg*gDs*(2*(t
     &   rsgc+2*mc2*(2*trs+trg))*gDs**2+(mc2*(4*(trs+trg)*cDs+3*trsgc+ms
     &   2*(11*trg+9*trc))+2*trsgc*cDs-2*mc2**2*trs)*gDs-2*mc2*((trsgc+m
     &   s2*trg)*cDs+ms2*(trsgc+mc2*trg)))+cDg**3*(12*trs*gDs**2+2*ms2*(
     &   7*trs+3*trg+5*trc)*gDs+ms2*(trsgc-2*mc2*trs+5*ms2*trg+3*ms2*trc
     &   ))+mc2*gDs**2*(mc2*(trg*(gDs+cDs+ms2)+trs*gDs)+trsgc*(gDs+cDs+m
     &   c2)))/(cDg**3*gDs)

      Ampp = 3*tr3c001fs*(cDg**3*(-4*((trs+trg)*cDs-2*ms2*trc)*gDs**2+6
     &   *ms2**2*trc*gDs+ms2**3*trc)-4*(trs+trg)*cDg**4*gDs*(3*gDs+ms2)+
     &   3*mc2*trc*cDg*gDs**2*(2*gDs+ms2)**2+3*mc2*(trs+trg)*cDg**2*gDs*
     &   (2*gDs+ms2)**2+mc2**2*(trs+trg)*gDs**3*(2*gDs+ms2))/(mc*cDg*gDs
     &   **3)+Ampp

      Ampp = 3*tr3c002fs*(cDg**3*(mc2*ms2*(trs+trg)*(4*gDs+ms2)-4*trc*c
     &   Ds*gDs**2)-4*cDg**4*gDs*(trc*gDs+(trs+trg)*cDs)+3*mc2*trc*cDg**
     &   2*gDs*(2*gDs+ms2)**2+mc2**2*trc*gDs**3*(2*gDs+ms2)+3*mc2**2*(tr
     &   s+trg)*cDg*gDs**2*(2*gDs+ms2)-4*(trs+trg)*cDg**5*gDs)/(mc*cDg*g
     &   Ds**3)+Ampp

      Ampp = 3*mc*tr3s002ft*(2*cDg+mc2)*(-trs*cDg*(3*(2*cDg+mc2)*gDs**2
     &   +ms2*cDg**2)-(trg+trc)*gDs*((2*cDg+mc2)*gDs**2+3*ms2*cDg**2))/(
     &   cDg**3*gDs)+Ampp

      Ampp = tr3s001ft*(-3*mc*ms2*(trg+trc)*cDg*(3*(2*cDg+mc2)*gDs**2+m
     &   s2*cDg**2)-3*mc*trs*(2*cDg+mc2)*gDs*((2*cDg+mc2)*gDs**2+3*ms2*c
     &   Dg**2))/(cDg**3*gDs)+Ampp

      Ampp = lc*mc*(cDg**2*(-4*trg*gDs**2+4*trg*cDs*gDs-2*trsgc*gDs+ms2
     &   *trsgc-mc2*ms2*trg)+mc2*(trsgc-mc2*trg)*gDs**2+2*cDg*gDs*((trsg
     &   c-mc2*trg)*gDs+2*mc2*trg*cDs)-2*cDg**3*(2*trs*gDs+ms2*trg))/(cD
     &   g*(2*cDg+mc2)*gDs)/2.0_dp+Ampp

      Ampp = epinv**2*mc*(2*trs*cDg*gDs-trsgc*gDs-ms2*trg*cDg)*(mc2*gDs
     &   **2-2*cDg*cDs*gDs+ms2*cDg**2)/(cDg**2*gDs**2)/2.0_dp+Ampp

      Ampp=Ampp/zmpp

      return
      end
