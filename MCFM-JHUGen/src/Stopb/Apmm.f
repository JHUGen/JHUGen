      subroutine Aamp_pmm(q,mc,ms,Apmm)
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
      complex(dp):: trc,trg,trs,trsgc,zpmm,Apmm

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
      zpmm=za(2,4)*za(2,3)*zb(4,3)**2

      Apmm = -ms*(cDg**3*(4*trs*gDs**3-2*ms2*trg*gDs**2+ms2*(trsgc-ms2*
     &   (trg+2*trc))*gDs+ms2**2*trsgc)+mc2**2*ms2*trg*gDs**3+cDg**2*gDs
     &   *(2*(trsgc+ms2*trg)*gDs**2+2*cDs*((-trsgc+ms2*trg+2*ms2*trc)*gD
     &   s-ms2*trsgc)+mc2*ms2**2*trg)+mc2*cDg*gDs**2*(2*trsgc*gDs-2*ms2*
     &   trc*gDs-2*ms2*trg*cDs+ms2*trsgc))*tr3Xs/(cDg*gDs**3)

      Apmm = mc2*ms*(-mc2**2*trg*gDs**3-cDg**2*gDs*(4*trc*cDs*gDs-2*trs
     &   gc*cDs+mc2*ms2*trg)+cDg**3*((trsgc+ms2*trg+2*ms2*trc)*gDs-ms2*t
     &   rsgc)+mc2*cDg*gDs**2*(2*trc*gDs+2*trg*cDs-trsgc))*tr3Xc/(cDg**3
     &   *gDs)+Apmm

      Apmm = ms*(2*trc*cDg*gDs-mc2*trg*gDs-trsgc*cDg)*(mc2*gDs**2-2*cDg
     &   *cDs*gDs+ms2*cDg**2)*tr1Xfs/(cDg**2*gDs)+Apmm

      Apmm = ms*(2*trc*cDg*gDs-mc2*trg*gDs-trsgc*cDg)*(mc2*gDs**2-2*cDg
     &   *cDs*gDs+ms2*cDg**2)*tr1Xfc/(cDg*gDs**2)+Apmm

      Apmm = Apmm-ms*(mc2*gDs**4*(ms2*(trg+trc)*(2*gDs+2*cDs+ms2-mc2)-
     &   trs*(4*gDs**2+ms2*(4*(gDs+cDs)-mc2)+8*cDs*gDs-2*mc2*gDs+4*cDs**
     &   2-2*mc2*cDs+ms2**2+mc2**2))+mc2*cDg**4*(-4*(2*trs+trg)*gDs**2+2
     &   *(trsgc+ms2*(trs+3*trg))*gDs+ms2*(-trsgc+3*ms2*trs+4*ms2*trg))+
     &   2*mc2*cDg**3*gDs*(-2*(trs+trc)*gDs**2+(2*(trs+trg)*cDs-3*trsgc+
     &   ms2*(trc-2*trg))*gDs+2*trsgc*cDs+ms2**2*(trs+trg+trc))+4*trc*cD
     &   g**5*gDs**2+mc2*cDg*gDs**3*((ms2*(6*trs+4*(trg+trc))-2*mc2*trs)
     &   *gDs+6*ms2*(trs+trg+trc)*cDs+3*ms2*(ms2-mc2)*(trs+trg+trc))-mc2
     &   *cDg**2*gDs**2*(2*(mc2*(2*trs+trg)+ms2*(-trs+2*trg+trc))*gDs+2*
     &   ms2*((trg-2*trs)*cDs+mc2*trs+ms2*(-trs+2*trg+trc))+3*mc2*trsgc)
     &   +4*cDg**6*((trs+trg)*gDs-ms2*trc))*lVs/(mc2*cDg**3*gDs**2)/4.0d
     &   +0

      Apmm = (cDg**4*(16*trs*gDs**3+4*ms2*trc*gDs**2+2*ms2*(2*trsgc+3*m
     &   c2*(trs+trg)-ms2*(trg+2*trc))*gDs+ms2**2*(2*trsgc+3*mc2*(trs+tr
     &   g)))+mc2*gDs**4*(ms2*(trg+trc)*(2*gDs+2*cDs+ms2-mc2)-trs*(4*gDs
     &   **2+ms2*(4*(gDs+cDs)-mc2)+8*cDs*gDs-2*mc2*gDs+4*cDs**2-2*mc2*cD
     &   s+ms2**2+mc2**2))+cDg**3*gDs*(4*(trsgc+2*mc2*trs+ms2*trg)*gDs**
     &   2+4*ms2*cDs*(2*(trg+trc)*gDs-trsgc)+6*mc2*ms2*trc*gDs+mc2*ms2*(
     &   ms2*(2*trs+trg+2*trc)-trsgc))+mc2*cDg**2*gDs**2*(ms2*((2*trs-4*
     &   trg-6*trc)*gDs-(trg-2*trs)*(2*cDs-mc2)-2*ms2*(-trs+2*trg+trc))+
     &   2*trsgc*(2*gDs+ms2))+mc2*cDg*gDs**3*(-2*(mc2-3*ms2)*trs*gDs+4*m
     &   s2*(trg+trc)*gDs+6*ms2*(trs+trg+trc)*cDs+mc2*trsgc+ms2*(3*ms2*(
     &   trs+trg+trc)-mc2*(3*trs+2*trg+3*trc)))-4*ms2*trc*cDg**6)*lVc/(m
     &   s*cDg**3*gDs**2)/4.0_dp+Apmm

      Apmm = ms*(2*trc*cDg*gDs-mc2*trg*gDs-trsgc*cDg)*(mc2*gDs**2-2*cDg
     &   *cDs*gDs+ms2*cDg**2)*lRs2/(cDg**2*gDs**2)/2.0_dp+Apmm

      Apmm = epinv*ms*(2*trc*cDg*gDs-mc2*trg*gDs-trsgc*cDg)*(mc2*gDs**2
     &   -2*cDg*cDs*gDs+ms2*cDg**2)*(2*lRs1+2*lRc1-1)/(cDg**2*gDs**2)/4.
     &   0_dp+Apmm

      Apmm = ms*(2*trc*cDg*gDs-mc2*trg*gDs-trsgc*cDg)*(mc2*gDs**2-2*cDg
     &   *cDs*gDs+ms2*cDg**2)*lRc2/(cDg**2*gDs**2)/2.0_dp+Apmm

      Apmm = 2*LsA*ms*(mc2**3*trg*gDs**5+cDg**5*(4*trs*gDs**3-ms2*gDs*(
     &   2*trg*gDs+ms2*(trg+2*trc))+ms2*trsgc*(gDs+ms2))+cDg**4*gDs*(ms2
     &   *(2*(trg+3*trc)*cDs*gDs+mc2*ms2*trg)-trsgc*(2*gDs*(gDs+cDs)+3*m
     &   s2*cDs))+cDg**3*gDs**2*(mc2*((trsgc-ms2*(trg+4*trc))*gDs+2*ms2*
     &   trsgc)+cDs**2*(2*trsgc-4*trc*gDs)-3*mc2*ms2*trg*cDs)+mc2*cDg**2
     &   *gDs**3*(6*trc*cDs*gDs+2*trg*cDs**2-3*trsgc*cDs+2*mc2*ms2*trg)+
     &   mc2**2*cDg*gDs**4*(-2*trc*gDs-3*trg*cDs+trsgc))/(cDg**3*gDs**3)
     &   +Apmm

      Apmm = (-mc2*cDg**3*gDs*(2*gDs+ms2)*(136*trs*gDs**4+6*cDs*(24*trs
     &   *gDs**3+2*ms2*(-trs+2*trg+5*trc)*gDs**2+3*ms2*(5*trsgc+ms2*trg+
     &   2*mc2*trg-2*ms2*trc)*gDs+3*ms2**2*(trsgc+mc2*trg))-12*ms2*(7*tr
     &   s+3*trg)*gDs**3+6*ms2*(-12*trsgc+mc2*(-7*trs+trg+trc)-23*ms2*tr
     &   g-15*ms2*trc)*gDs**2+ms2*(mc2*(27*trsgc+ms2*(22*trs+101*trg-10*
     &   trc))-6*ms2*(ms2*(trs+trg+trc)-3*trsgc))*gDs+18*mc2*ms2**2*(trs
     &   gc+ms2*trg))-ms2*cDg**5*(-12*(ms2*trc-2*mc2*(7*trs+4*trg))*gDs*
     &   *3+4*mc2*(ms2*(29*trs+2*(trg+trc))-9*trsgc)*gDs**2+16*cDs*gDs*(
     &   2*gDs+ms2)*(mc2*(trs+trg)-2*trc*gDs)+2*mc2*ms2*(ms2*(16*trs-11*
     &   trg+12*trc)-18*trsgc)*gDs+mc2*ms2**2*(ms2*(8*trs-trg+10*trc)-9*
     &   trsgc))-mc2*gDs**4*(2*gDs+ms2)*(-12*(mc-ms)*(ms+mc)*trs*gDs**3-
     &   24*(mc-ms)*(ms+mc)*trs*cDs*gDs**2+6*(mc-ms)*(ms+mc)*(mc2*trs+ms
     &   2*(-2*trs+trg+trc))*gDs**2+9*mc2**2*trsgc*(gDs+cDs+ms2)-12*(mc-
     &   ms)*(ms+mc)*trs*cDs**2*gDs+6*(mc-ms)*(ms+mc)*(mc2*trs+ms2*(-2*t
     &   rs+trg+trc))*cDs*gDs+mc2**2*ms2*(6*trs-20*trg-29*trc)*gDs+5*mc2
     &   **3*trs*gDs-3*ms2**3*(-trs+trg+trc)*gDs+6*mc2*ms2**2*(-trs+trg+
     &   trc)*gDs+9*mc2**2*ms2*trg*cDs+9*mc2**3*ms2*trg)+mc2*ms2*cDg**4*
     &   (2*gDs+ms2)*(108*trs*gDs**3+9*cDs*((8*trc-4*trs)*gDs**2+2*(trsg
     &   c+2*ms2*(trg+trc))*gDs+ms2*(trsgc+ms2*trg))+2*(9*trsgc+ms2*(11*
     &   trs-5*trg-17*trc)+15*mc2*(trs+trg))*gDs**2+3*ms2*(12*trsgc+ms2*
     &   (3*(trs+trg)-14*trc)+5*mc2*trs+11*mc2*trg)*gDs+9*ms2**2*(trsgc+
     &   mc2*trg))-mc2*cDg**2*gDs**3*(2*gDs+ms2)*(4*(9*trsgc+24*mc2*trs+
     &   ms2*(trg-8*trc))*gDs**2+6*cDs*(6*(trsgc+2*mc2*trs+3*ms2*trg+2*m
     &   s2*trc)*gDs+ms2*(-6*trsgc+2*mc2*(trs-2*trg+3*trc)-2*ms2*trs+ms2
     &   *trg))+2*ms2*(3*ms2*(-trs+2*trg+trc)+mc2*(-39*trs-6*trg+14*trc)
     &   )*gDs+ms2*(-mc2*(99*trsgc+2*ms2*(-6*trs+60*trg+11*trc))+mc2**2*
     &   (10*trs-23*trg)+6*ms2**2*(-trs+2*trg+trc)))+mc2*cDg*gDs**3*(2*g
     &   Ds+ms2)*(2*(2*mc2*(ms2*(-6*trs+5*trg+14*trc)-9*trsgc)+3*ms2**2*
     &   (3*trs+2*(trg+trc))-12*mc2**2*trs)*gDs**2-18*cDs*((mc2*(2*trsgc
     &   +ms2*(trs+5*trg+3*trc))-ms2**2*(trs+trg+trc))*gDs-mc2*ms2*(trsg
     &   c+mc2*trg))+3*ms2*(-6*mc2*(trsgc+ms2*(trs+trg+trc))+mc2**2*(11*
     &   trs-5*trg+trc)+3*ms2**2*(trs+trg+trc))*gDs+18*mc2**2*ms2*(trsgc
     &   +ms2*trg))+4*ms2*cDg**6*gDs*(8*(trs+trg-trc)*cDs*(2*gDs+ms2)+(m
     &   s2*(3*(trs+trg)+10*trc)+6*mc2*(trc-4*(trs+trg)))*gDs+ms2*(3*mc2
     &   *(trc-4*(trs+trg))+5*ms2*trc))-32*ms2*(-trs-trg+trc)*cDg**7*gDs
     &   *(2*gDs+ms2))/(mc2*ms*cDg**3*gDs**3*(2*gDs+ms2))/1.2e+1_dp+Apmm

      Apmm = B0cgsf*(-2*mc2*cDg**3*gDs*(8*trs*gDs**4+cDs*(8*trs*gDs**3+
     &   2*ms2*(-trs+2*trg+4*trc)*gDs**2+ms2*(trsgc+ms2*(trg-2*trc)+2*mc
     &   2*trg)*gDs+ms2**2*(trsgc+mc2*trg))+2*(trsgc+2*mc2*trs+ms2*(trc-
     &   2*trs))*gDs**3-ms2*(trsgc-2*mc2*(trc-trs)+6*ms2*(trg+trc))*gDs*
     &   *2+ms2*(ms2*(trsgc+2*mc2*(trs+3*trg))+mc2*trsgc-ms2**2*(trs+trg
     &   +trc))*gDs+mc2*ms2**2*(trsgc+ms2*trg))+mc2*cDg**4*(-16*trs*gDs*
     &   *4-4*ms2*(-trs+trg+trc)*gDs**3+ms2*cDs*((8*trc-4*trs)*gDs**2+2*
     &   (trsgc+2*ms2*(trg+trc))*gDs+ms2*(trsgc+ms2*trg))+2*ms2**2*(2*tr
     &   s+3*trg)*gDs**2+ms2**2*(trsgc+3*ms2*trs+4*ms2*trg+2*mc2*trg-5*m
     &   s2*trc)*gDs+ms2**3*(trsgc+mc2*trg))-mc2*gDs**4*(4*(ms2-mc2)*trs
     &   *gDs**3+cDs*(-8*(mc-ms)*(ms+mc)*trs*gDs**2+2*(mc-ms)*(ms+mc)*(m
     &   c2*trs+ms2*(-2*trs+trg+trc))*gDs+mc2**2*(trsgc+ms2*trg))+2*(mc-
     &   ms)*(ms+mc)*(mc2*trs+ms2*(-2*trs+trg+trc))*gDs**2+4*(ms2-mc2)*t
     &   rs*cDs**2*gDs+(mc2**2*(trsgc+ms2*(2*trs-3*trg-4*trc))-ms2**3*(-
     &   trs+trg+trc)+2*mc2*ms2**2*(-trs+trg+trc))*gDs+mc2**2*ms2*(trsgc
     &   +mc2*trg))-ms2*cDg**5*(-4*trc*gDs**3+2*mc2*(5*trs+3*trg)*gDs**2
     &   +2*cDs*gDs*(mc2*(trs+trg)-2*trc*gDs)-2*mc2*(trsgc+ms2*(trg-trs)
     &   )*gDs+mc2*ms2*(ms2*(trs+trc)-trsgc))-mc2*cDg*gDs**3*(2*(mc2*(2*
     &   trsgc+4*ms2*trs-2*ms2*trc)-ms2**2*(3*trs+2*(trg+trc))+mc2**2*tr
     &   s)*gDs**2+2*cDs*((mc2*(2*trsgc+3*ms2*trs+7*ms2*trg+5*ms2*trc)-3
     &   *ms2**2*(trs+trg+trc))*gDs-mc2*ms2*(trsgc+mc2*trg))+(2*mc2*ms2*
     &   (trsgc+3*ms2*(trs+trg+trc))+mc2**2*(trsgc+ms2*(-6*trs+trg-2*trc
     &   ))-3*ms2**3*(trs+trg+trc))*gDs-2*mc2**2*ms2*(trsgc+ms2*trg))-2*
     &   mc2*cDg**2*gDs**3*(-2*(-trsgc-3*mc2*trs+ms2*trc)*gDs**2+cDs*(2*
     &   (trsgc+2*mc2*trs+3*ms2*trg+2*ms2*trc)*gDs+ms2*(-2*trsgc+2*mc2*(
     &   trs-trg+trc)-2*ms2*trs+ms2*trg))-(mc2*(ms2*(3*trs+trg+trc)-2*tr
     &   sgc)+ms2**2*(trs-2*trg-trc))*gDs-ms2*(mc2*(3*trsgc+2*ms2*(-trs+
     &   4*trg+trc))+ms2**2*(trs-2*trg-trc)+mc2**2*trg))+2*ms2*cDg**6*gD
     &   s*(2*(trs+trg)*gDs+2*(trs+trg-trc)*cDs+mc2*(2*trc-3*(trs+trg)))
     &   -4*ms2*(-trs-trg+trc)*cDg**7*gDs)/(mc2*ms*cDg**3*gDs**3)/4.0_dp
     &   +Apmm

      Apmm = tr3s00ft*(cDg**2*(-4*(-trsgc+2*ms2*trg+3*ms2*trc)*gDs**3+c
     &   Ds*(4*(trsgc+2*mc2*trs+3*ms2*trg+2*ms2*trc)*gDs**2+4*ms2*(-trsg
     &   c-mc2*trg+mc2*trc)*gDs+mc2*ms2*(trsgc+ms2*trg))-12*mc2*ms2*(3*t
     &   rs+trg+trc)*gDs**2-mc2*ms2*(9*trsgc+18*ms2*trg+11*ms2*trc)*gDs+
     &   mc2*ms2**2*(trsgc+mc2*trg))+cDg**3*(8*trs*gDs**3+2*cDs*(8*trs*g
     &   Ds**2+2*ms2*(trg+3*trc)*gDs+ms2*(2*trsgc+ms2*trg))-4*ms2*(9*trs
     &   +4*trg+3*trc)*gDs**2-2*ms2*(4*trsgc+6*mc2*trs+14*ms2*trg+11*ms2
     &   *trc)*gDs+ms2**2*(2*trsgc-3*mc2*trs+3*mc2*trg-2*mc2*trc))-mc2*c
     &   Dg*gDs*(2*(-2*trsgc+mc2*trs+6*ms2*trg+8*ms2*trc)*gDs**2+cDs*(2*
     &   ms2*(trsgc+mc2*trg)-4*(trsgc+2*ms2*trg+ms2*trc)*gDs)+ms2*(-2*tr
     &   sgc+9*mc2*trs-mc2*trg+2*mc2*trc)*gDs+2*mc2*ms2*(trsgc+ms2*trg))
     &   +mc2**2*gDs**2*((trsgc-4*ms2*trg-5*ms2*trc)*gDs+(trsgc+ms2*trg)
     &   *cDs+ms2*(trsgc+mc2*trg))-6*ms2*cDg**4*(4*trs*gDs+ms2*(trs+trc)
     &   ))/(ms*cDg**3*gDs)+Apmm

      Apmm = ms*tr3c00fs*(mc2*cDg**2*(24*trc*gDs**3-cDs*((8*trc-4*trs)*
     &   gDs**2+2*(trsgc+2*ms2*(trg+trc))*gDs+ms2*(trsgc+ms2*trg))+(30*m
     &   s2*trc-2*trsgc)*gDs**2+ms2*(-3*trsgc-2*mc2*trg+11*ms2*trc)*gDs-
     &   ms2**2*(trsgc+mc2*trg))+cDg**3*(mc2*(2*(5*trs+3*trg)*gDs**2-2*t
     &   rsgc*gDs+2*ms2*(3*trs+trg+2*trc)*gDs-ms2*trsgc+ms2**2*(3*trs+2*
     &   (trg+trc)))-12*trc*cDs*gDs**2)+mc2*cDg*gDs*(4*mc2*(2*trs+3*trg+
     &   trc)*gDs**2+(4*(trsgc+mc2*trg-ms2*trc)*cDs+mc2*(3*trsgc+ms2*(9*
     &   trs+14*trg+2*trc)))*gDs+2*ms2*(mc2*trg*(cDs+ms2)+trsgc*(cDs+mc2
     &   )))+mc2**2*gDs**2*(6*trc*gDs**2+(-trsgc-2*mc2*trg+5*ms2*trc)*gD
     &   s-ms2*trg*(cDs+mc2)+trsgc*(-cDs-ms2))-12*cDg**4*gDs*(2*trc*gDs+
     &   (trs+trg)*cDs+ms2*trc)-12*(trs+trg)*cDg**5*gDs)/(mc2*cDg*gDs**3
     &   )+Apmm

      Apmm = BfunX*(ms2*cDg**5*(-4*trc*gDs**3-10*mc2*(trs+trg)*gDs**2-2
     &   *cDs*gDs*(2*trc*gDs+mc2*(trs+trg))+4*mc2*ms2*(-trs-trg+trc)*gDs
     &   +mc2*ms2**2*(trs+trg+trc))-mc2*ms2*cDg**4*gDs*(4*(3*(trs+trg)+t
     &   rc)*gDs**2+2*(6*(trs+trg)*cDs+ms2*(7*(trs+trg)-2*trc))*gDs+3*ms
     &   2*(2*(trs+trg)*cDs+ms2*(trs+trg-trc)))-2*mc2*ms2*cDg**3*gDs**2*
     &   ((6*trs+4*trc)*gDs**2+2*cDs*((3*trs+2*trc)*gDs+ms2*(trs+trg+trc
     &   ))+ms2*(4*trs+trg+4*trc)*gDs+ms2**2*(trs+trg+trc))+mc2*gDs**5*(
     &   2*(gDs+cDs)+ms2)**2*(2*trs*(gDs+cDs)-ms2*(-trs+trg+trc))+mc2*cD
     &   g*gDs**4*(2*trs*gDs-3*ms2*(trs+trg+trc))*(2*(gDs+cDs)+ms2)**2+2
     &   *mc2*ms2*cDg**2*gDs**3*(-(4*trs+trg+2*trc)*gDs+(trg-2*trs)*cDs+
     &   ms2*(-trs+2*trg+trc))*(2*(gDs+cDs)+ms2)-2*ms2*cDg**6*gDs*(2*(tr
     &   s+trg)*gDs+2*(trs+trg-trc)*cDs+mc2*(3*(trs+trg)-2*trc))+4*ms2*(
     &   -trs-trg+trc)*cDg**7*gDs)/(mc2*ms*cDg**3*gDs**3)/4.0_dp+Apmm

      Apmm = ls*ms*(8*trc*cDg**3*gDs**3+mc2*cDg**2*(2*gDs+ms2)*(trsgc*(
     &   2*gDs+ms2)-trg*(4*gDs**2+2*ms2*gDs+ms2**2))-2*mc2*cDg*gDs*(2*gD
     &   s+ms2)*(2*trc*gDs**2+(trsgc-2*trg*cDs)*gDs-2*ms2*trg*cDs)+8*(tr
     &   s+trg)*cDg**4*gDs**2+mc2**2*gDs**2*(2*gDs+ms2)*(trsgc-trg*(2*gD
     &   s+ms2)))/(mc2*cDg*gDs*(2*gDs+ms2)**2)/2.0_dp+Apmm

      Apmm = 3*ms*tr3c002fs*(mc2**3*(trs+trg)*gDs**3+mc2*cDg**3*(ms2**2
     &   *trc-2*((trs+trg)*cDs-2*ms2*trc)*gDs)+3*mc2**2*trc*cDg*gDs**2*(
     &   2*gDs+ms2)+3*mc2**2*(trs+trg)*cDg**2*gDs*(2*gDs+ms2)-2*cDg**4*(
     &   2*trc*cDs+3*mc2*(trs+trg))*gDs-4*trc*cDg**5*gDs)/(mc2*cDg*gDs**
     &   3)+Apmm

      Apmm = 3*ms*tr3c001fs*(cDg**3*(mc2*ms2*(trs+trg)*(2*gDs+ms2)-4*tr
     &   c*cDs*gDs**2)-4*cDg**4*gDs*(3*trc*gDs+(trs+trg)*cDs+ms2*trc)+3*
     &   mc2*trc*cDg**2*gDs*(2*gDs+ms2)**2+mc2**2*trc*gDs**3*(2*gDs+ms2)
     &   +3*mc2**2*(trs+trg)*cDg*gDs**2*(2*gDs+ms2)-4*(trs+trg)*cDg**5*g
     &   Ds)/(mc2*cDg*gDs**3)+Apmm

      Apmm = tr3s002ft*(-3*ms2*(trg+trc)*cDg*(2*cDg+mc2)*(3*(2*cDg+mc2)
     &   *gDs**2+ms2*cDg**2)-3*trs*(2*cDg+mc2)**2*gDs*((2*cDg+mc2)*gDs**
     &   2+3*ms2*cDg**2))/(ms*cDg**3*gDs)+Apmm

      Apmm = 3*ms*tr3s001ft*(2*cDg+mc2)*(-trs*cDg*(3*(2*cDg+mc2)*gDs**2
     &   +ms2*cDg**2)-(trg+trc)*gDs*((2*cDg+mc2)*gDs**2+3*ms2*cDg**2))/(
     &   cDg**3*gDs)+Apmm

      Apmm = lc*(-cDg**2*(4*(trsgc+2*mc2*trs+ms2*trg)*gDs**2+8*ms2*trg*
     &   cDs*gDs-2*mc2*ms2*trg*gDs+mc2*ms2*(trsgc-ms2*trg))-2*cDg**3*(8*
     &   trs*gDs**2+2*ms2*trc*gDs+ms2*(2*trsgc-ms2*trg))-mc2**2*(trsgc-m
     &   s2*trg)*gDs**2-4*mc2*cDg*gDs*(trsgc*gDs+ms2*trg*cDs))/(ms*cDg*(
     &   2*cDg+mc2)*gDs)/2.0_dp+Apmm

      Apmm = epinv**2*ms*(2*trc*cDg*gDs-mc2*trg*gDs-trsgc*cDg)*(mc2*gDs
     &   **2-2*cDg*cDs*gDs+ms2*cDg**2)/(cDg**2*gDs**2)/2.0_dp+Apmm

      Apmm=Apmm/zpmm

      return
      end
