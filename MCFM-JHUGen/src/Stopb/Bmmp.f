      subroutine Bamp_mmp(q,mc,ms,Bmmp)
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
      complex(dp):: trc,trg,trs,trsgc,zmmp,Bmmp

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
      zmmp=za(4,3)*zb(2,3)**2

      Bmmp = ms2*(cDg**3*gDs*(mc2*ms2*(4*trs*gDs-2*trc*gDs+trsgc-ms2*tr
     &   g)+2*cDs**2*(-2*trs*gDs+2*trc*gDs-2*trsgc+ms2*trg)+ms2*cDs*(-6*
     &   trc*gDs+5*trsgc+mc2*trg))+cDg**2*gDs**2*(mc2*ms2*(4*trs*gDs-6*t
     &   rc*gDs+trsgc-4*ms2*trg+mc2*trg)+2*cDs**2*(-2*trs*gDs+6*trc*gDs-
     &   4*trsgc+3*ms2*trg-mc2*trg)+cDs*(-6*ms2*trc*gDs+4*ms2*trsgc+mc2*
     &   trsgc+mc2*ms2*trg))+cDg*gDs**3*(cDs*(-2*ms2*trc*gDs+ms2*trsgc+2
     &   *mc2*trsgc+7*mc2*ms2*trg)+2*cDs**2*(6*trc*gDs-3*trsgc+ms2*trg-3
     &   *mc2*trg)+mc2*ms2*(-6*trc*gDs+trsgc-ms2*trg+4*mc2*trg)-8*trg*cD
     &   s**3)+gDs**4*(cDs**2*(4*trc*gDs-2*(trsgc+mc2*trg))+mc2*ms2*(-2*
     &   trc*gDs+trsgc+mc2*trg)-4*trg*cDs**3+mc2*(trsgc+3*ms2*trg)*cDs)+
     &   2*ms2*cDg**4*cDs*(trsgc-trc*gDs))*tr5Xs/((cDs**2-mc2*ms2)*gDs**
     &   3*(gDs+cDg))/2.0_dp

      Bmmp = mc2*(2*mc2*trsgc*cDs*gDs**4+cDg**3*gDs*(cDs*(8*trs*gDs**2-
     &   2*(trsgc+3*mc2*trs+2*ms2*trg)*gDs+2*ms2*trsgc+mc2*trsgc+7*mc2*m
     &   s2*trg)+ms2*(2*gDs*(-2*trs*gDs+2*trsgc+ms2*trg)+mc2*(-6*trs*gDs
     &   +trsgc+2*ms2*trg)+mc2**2*trg)-2*cDs**2*(-6*trs*gDs+3*trsgc+2*ms
     &   2*trg)-8*trg*cDs**3)+cDg**2*gDs**2*(cDs*(4*trs*gDs**2-2*(2*trsg
     &   c+3*mc2*trs+ms2*trg)*gDs+ms2*trsgc+4*mc2*trsgc+5*mc2*ms2*trg)+m
     &   s2*(mc2*(-2*trs*gDs+3*trsgc+ms2*trg)+2*trsgc*gDs-2*mc2**2*trg)+
     &   2*cDs**2*(2*trs*gDs-5*trsgc-ms2*trg+2*mc2*trg)-4*trg*cDs**3)+cD
     &   g**4*(cDs*(4*trs*gDs**2-2*(mc2*trs+ms2*trg)*gDs+ms2*(trsgc-mc2*
     &   trg))+ms2*(-8*trs*gDs**2+2*(trsgc-mc2*trs+2*ms2*trg)*gDs-mc2*(t
     &   rsgc+ms2*trg))+8*trs*cDs**2*gDs)+cDg*gDs**3*(cDs*(-2*trsgc*gDs-
     &   2*mc2*trs*gDs+5*mc2*trsgc+mc2*ms2*trg)+2*(mc2*trg-2*trsgc)*cDs*
     &   *2+mc2*ms2*(trsgc-mc2*trg))+2*ms2*cDg**5*(trs*(mc2-2*gDs)+ms2*t
     &   rg))*tr5Xc/(cDg**2*(cDs**2-mc2*ms2)*gDs*(gDs+cDg))/2.0_dp+Bmmp

      Bmmp = ms2*(trsgc*(cDs*gDs*(2*cDs*gDs-mc2*gDs+4*cDg*cDs)-ms2*(cDg
     &   *cDs*(gDs+2*cDg)+mc2*gDs*(gDs+cDg)))+trg*gDs*(-ms2*(mc2*(3*cDs+
     &   mc2)*gDs+cDg*cDs*(2*cDs+mc2))+2*cDs**2*(2*cDs+mc2)*gDs+mc2*ms2*
     &   *2*cDg)-2*trc*gDs*(gDs+cDg)*(2*cDs**2*gDs-ms2*(mc2*gDs+cDg*cDs)
     &   ))*tr4Xs/((cDs-mc*ms)*(cDs+mc*ms)*gDs**2)/2.0_dp+Bmmp

      Bmmp = mc2*(trsgc*(-ms2*cDg*(cDg*(2*gDs+cDs)-mc2*(gDs+cDg))-cDs*g
     &   Ds*(mc2*(2*gDs+cDg)-2*cDg*(gDs+cDs)))-ms2*trg*cDg*(-cDg*(cDs*(2
     &   *gDs+mc2)+mc2*ms2)+mc2*(cDs+mc2)*gDs+2*ms2*cDg**2)+2*trs*cDg*(m
     &   s2*cDg-cDs*gDs)*(2*cDg*gDs-mc2*(gDs+cDg)))*tr4Xc/(cDg*(cDs-mc*m
     &   s)*(cDs+mc*ms)*gDs)/2.0_dp+Bmmp

      Bmmp = (cDg**2*(mc2*ms2*(4*trs*gDs**2+2*(trsgc+mc2*trs+2*ms2*trg)
     &   *gDs+ms2*(mc2*trg-trsgc))+ms2*cDs*(mc2*(-2*trs*gDs+trsgc+ms2*tr
     &   g)-2*trsgc*gDs)-2*cDs**2*gDs*(4*trs*gDs+trsgc+ms2*trg))+cDg*gDs
     &   *(-mc2*ms2*(4*trg*gDs**2+2*mc2*(trg-trs)*gDs+mc2*(trsgc+ms2*trg
     &   ))+2*cDs**2*gDs*(2*trg*gDs+trsgc+mc2*trg)+mc2*ms2*(3*trsgc-mc2*
     &   trg)*cDs)-2*mc2**2*ms2*trsgc*gDs**2+cDg**3*(-2*ms2*cDs*(-2*trs*
     &   gDs+mc2*trs+ms2*trg)+4*trs*cDs**2*gDs-4*mc2*ms2*trs*gDs))*tr3Xs
     &   /(cDg*(cDs**2-mc2*ms2)*gDs)/2.0_dp+Bmmp

      Bmmp = mc2*(trg*gDs*(-4*cDs**2*gDs**2+ms2**2*(mc2*gDs-cDg*(cDs+mc
     &   2))+mc2*ms2*gDs*(4*gDs+cDs))-ms2*trsgc*(gDs*(mc2*gDs-cDs*(gDs+3
     &   *cDg))+ms2*cDg*(gDs+2*cDg))-2*ms2*trc*gDs*(gDs+cDg)*(cDs*gDs-ms
     &   2*cDg))*tr3Xc/((cDs-mc*ms)*(cDs+mc*ms)*gDs**2)/2.0_dp+Bmmp

      Bmmp = (-cDg**3*(cDs**2*(4*(3*trs+trc)*gDs**2+(-2*trsgc-5*ms2*trs
     &   -2*ms2*trg+2*ms2*trc)*gDs+ms2*(trsgc-ms2*trg))+mc2*ms2*(-4*(2*t
     &   rs+trc)*gDs**2+(2*trsgc+3*ms2*trs+4*ms2*trg)*gDs+ms2*(ms2*trg-t
     &   rsgc))+2*ms2*(trsgc-mc2*trs+ms2*trc)*cDs*gDs)+cDg**2*gDs*(cDs**
     &   2*(-4*(trc-2*trs)*gDs**2+(6*trsgc+ms2*(3*trs+2*(trg+trc))+2*mc2
     &   *trg)*gDs+ms2*(trsgc+ms2*trs+mc2*trg))-mc2*ms2*(-4*(trc-2*trs)*
     &   gDs**2+(4*trsgc+2*mc2*(trs+trg)+3*ms2*trs+2*ms2*trg)*gDs+ms2*(t
     &   rsgc+ms2*trs+mc2*trg))-mc2*ms2*cDs*((2*trs+6*trg-2*trc)*gDs+4*t
     &   rsgc+ms2*trs)+cDs**3*(8*trg*gDs+4*trsgc+ms2*trs))-2*cDg*(cDs-mc
     &   *ms)*(cDs+mc*ms)*gDs**2*(trg*gDs*(2*gDs+mc2)+trsgc*(gDs+mc2)+ms
     &   2*trs*(-gDs-cDs)+ms2**2*(trg+trc))+ms2*(trs+trg+trc)*(cDs-mc*ms
     &   )*(cDs+mc*ms)*gDs**3*(2*gDs+2*cDs+ms2-mc2)+4*ms2*trs*cDg**4*cDs
     &   *gDs)*lVs/(cDg**2*(cDs-mc*ms)*(cDs+mc*ms)*gDs**2)/4.0_dp+Bmmp

      Bmmp = (cDg**2*gDs*(-mc2*ms2*(8*trs*gDs**3+2*(trsgc+ms2*trg)*gDs*
     &   *2+ms2*(4*trsgc-3*mc2*trs+2*mc2*trg-4*mc2*trc)*gDs+mc2*ms2*(trs
     &   gc-ms2*trs+ms2*trg))+cDs**2*(8*trs*gDs**3+2*(trsgc+ms2*trg)*gDs
     &   **2+ms2*(4*trsgc-5*mc2*trs+4*mc2*trg-2*mc2*trc)*gDs+mc2*ms2*(tr
     &   sgc-ms2*trs+ms2*trg))+cDs**3*(8*trs*gDs**2+2*(trsgc+ms2*trg)*gD
     &   s+4*ms2*trsgc-mc2*ms2*trs)+mc2*ms2*cDs*(-8*trs*gDs**2-2*mc2*trs
     &   *gDs-2*ms2*trg*gDs+2*ms2*trc*gDs-4*ms2*trsgc+mc2*ms2*trs))-2*cD
     &   g**3*(gDs*(mc2*ms2*(-2*trs*gDs**2-2*ms2*(2*trs*gDs+mc2*(trs+trg
     &   ))+ms2**2*(2*trg+trc))+cDs**2*(2*trs*gDs*(gDs+2*ms2)+mc2*ms2*(t
     &   rs+2*trg)-2*ms2**2*trg)+2*trs*cDs**3*gDs+mc2*ms2**2*(-trs+trg+t
     &   rc)*cDs)+ms2**2*trsgc*(mc2*(gDs-ms2)+cDs**2))-cDg*(cDs-mc*ms)*(
     &   cDs+mc*ms)*gDs**2*(mc2*(ms2*(2*trs*(gDs+cDs)+trsgc-2*ms2*(trg+t
     &   rc))+2*trg*gDs*(gDs+cDs))+2*gDs*(gDs+cDs)*(2*trg*gDs+trsgc)-mc2
     &   **2*ms2*trg)+mc2*ms2*(trs+trg+trc)*(mc*ms-cDs)*(cDs+mc*ms)*gDs*
     &   *3*(2*gDs+2*cDs+ms2-mc2)+4*mc2*ms2**2*trs*cDg**4*gDs)*lVc/(ms2*
     &   cDg**2*(cDs-mc*ms)*(cDs+mc*ms)*gDs**2)/4.0_dp+Bmmp

      Bmmp = xs*cDs*(cDg**2*(4*trs*gDs**2-2*ms2*trg*gDs+ms2*trsgc)+mc2*
     &   trsgc*gDs**2-2*trsgc*cDg*gDs*(gDs+cDs))*lRcs/(mc*ms*(xs**2-1)*c
     &   Dg*gDs**2)+Bmmp

      Bmmp = tr3c00fs*(gDs*(-mc2*(2*(2*trs+5*trg)*gDs**2+(2*trsgc+3*ms2
     &   *trs+8*ms2*trg)*gDs+ms2*(trsgc+ms2*trg))+cDs*(2*gDs+ms2)*(2*trc
     &   *gDs-mc2*trg)-trsgc*cDs*(4*gDs+ms2))+cDg*(cDs*(4*(trg+trc)*gDs*
     &   *2+2*(trsgc+ms2*(3*trs+5*trg+trc))*gDs+ms2*(trsgc+ms2*trg))+trs
     &   gc*(2*gDs**2+4*ms2*gDs+ms2**2)-(gDs+ms2)*(2*gDs+ms2)*(2*trc*gDs
     &   -mc2*trg-3*ms2*trc))-cDg**2*(4*(3*trs+2*trg+trc)*gDs**2-2*trsgc
     &   *gDs+2*ms2*(3*trs+trg+3*trc)*gDs-ms2*trsgc+ms2**2*(3*trs+2*(trg
     &   +trc))))/gDs**3+Bmmp

      Bmmp = Bmmp-epinv*(mc*ms*(xs**2-1)-4*lp*xs*cDs)*(cDg**2*(4*trs*g
     &   Ds**2-2*ms2*trg*gDs+ms2*trsgc)+mc2*trsgc*gDs**2-2*trsgc*cDg*gDs
     &   *(gDs+cDs))/(mc*ms*(xs**2-1)*cDg*gDs**2)/4.0_dp

      Bmmp = tr2fu*cDs*(cDg**2*(4*trs*gDs**2-2*ms2*trg*gDs+ms2*trsgc)+m
     &   c2*trsgc*gDs**2-2*trsgc*cDg*gDs*(gDs+cDs))/(cDg*gDs**2)+Bmmp

      Bmmp = Bmmp-3*tr3c001fs*(-ms2*cDg*(4*trc*gDs**2+2*((trs+trg)*cDs
     &   +2*ms2*trc)*gDs+ms2**2*trc)+(trs+trg)*cDg**2*(4*gDs**2+2*ms2*gD
     &   s+ms2**2)+mc2*(trs+trg)*gDs**2*(2*gDs+ms2))/gDs**3

      Bmmp = Bmmp-LsB2*(cDg*(cDs*(mc2*(-2*trs*gDs+trsgc+ms2*trg)-2*trs
     &   gc*gDs)+2*(mc2*trg-trsgc)*cDs**2+mc2*ms2*(trsgc-mc2*trg))+2*mc2
     &   *trsgc*cDs*gDs-2*cDg**2*cDs*(trs*(mc2-2*gDs)+ms2*trg))*(mc2*gDs
     &   **2-2*cDg*cDs*gDs+ms2*cDg**2)/(cDg**2*(cDs**2-mc2*ms2)*gDs)

      Bmmp = LsB1*(gDs**3*(cDs**3*(8*trc*gDs-4*(trsgc+mc2*trg))+3*mc2*m
     &   s2*cDs*(-2*trc*gDs+trsgc+mc2*trg)-8*trg*cDs**4+8*mc2*ms2*trg*cD
     &   s**2+mc2**2*ms2*(trsgc-ms2*trg))+2*cDg*gDs**2*(cDs**3*(4*trc*gD
     &   s-2*trsgc+4*ms2*trg)-3*mc2*ms2*cDs*(trc*gDs+ms2*trg)+2*ms2*cDs*
     &   *2*(-2*trc*gDs+trsgc+mc2*trg)-mc2*ms2**2*(-2*trc*gDs+trsgc+mc2*
     &   trg))-ms2*cDg**2*gDs*(2*cDs**2*(4*trc*gDs-3*trsgc+ms2*trg)+ms2*
     &   cDs*(-2*trc*gDs+trsgc+mc2*trg)+mc2*ms2*(-4*trc*gDs+trsgc-ms2*tr
     &   g))+2*ms2**2*cDg**3*cDs*(trc*gDs-trsgc))/((cDs**2-mc2*ms2)*gDs*
     &   *3)+Bmmp

      Bmmp = 3*tr3c002fs*(cDg*(trc*cDs*(2*gDs+ms2)**2-4*mc2*(trs+trg)*g
     &   Ds**2)-mc2*trc*gDs**2*(2*gDs+ms2))/gDs**3+Bmmp

      Bmmp = (cDg**4*(144*trs*gDs**6+24*ms2*(11*trs+6*trc)*gDs**5+cDs*(
     &   2*gDs+ms2)*(144*trs*gDs**4-40*ms2*(3*trg+2*trc)*gDs**3-2*ms2*(-
     &   36*trsgc+19*ms2*trs+23*mc2*trs+52*ms2*trg+4*ms2*trc+2*mc2*trc)*
     &   gDs**2+2*ms2*(mc2*(-9*trsgc-8*ms2*trs-17*ms2*trg+7*ms2*trc)+ms2
     &   **2*(8*trc-9*trg))*gDs+mc2*ms2**2*(-9*trsgc-9*ms2*trg+8*ms2*trc
     &   ))+8*ms2*(15*trsgc+3*ms2*trs+mc2*(-25*trs+6*trg+10*trc)+21*ms2*
     &   trg+44*ms2*trc)*gDs**4+2*ms2*(6*ms2*(-10*trsgc-5*ms2*trs+6*ms2*
     &   trg+15*ms2*trc)-2*mc2**2*(10*trs+19*trg)+3*mc2*ms2*(-7*trs+32*t
     &   rg+48*trc))*gDs**3-ms2**2*(6*ms2*(21*trsgc+ms2*trs+2*ms2*trc)-m
     &   c2*(-72*trsgc+39*ms2*trs+76*ms2*trg+180*ms2*trc)+4*mc2**2*(5*tr
     &   s+23*trg))*gDs**2+72*trs*cDs**2*gDs**3*(2*gDs+ms2)-ms2**3*(mc2*
     &   (54*trsgc+ms2*(-8*trs+trg-12*trc))+2*ms2*(9*trsgc+8*ms2*trc)+45
     &   *mc2**2*trg)*gDs-mc2*ms2**4*(9*(trsgc+mc2*trg)+8*ms2*trc))+cDg*
     &   *3*gDs*(-144*trs*gDs**6-72*(trsgc-7*ms2*trs-mc2*trs+5*ms2*trg)*
     &   gDs**5-cDs*(2*gDs+ms2)*(144*trs*gDs**4+72*(trsgc-mc2*trs+ms2*(t
     &   rg+trc))*gDs**3+2*ms2*(9*(-6*trsgc+ms2*trs+2*ms2*trc)+mc2*(19*t
     &   rs+30*trg+20*trc))*gDs**2+ms2*(-18*ms2*trsgc+mc2*(-54*trsgc+19*
     &   ms2*trs+16*ms2*trg+4*ms2*trc)+mc2**2*(23*trs-18*trg))*gDs+mc2*m
     &   s2**2*(9*ms2*trg-9*mc2*trg-8*ms2*trc))+4*ms2*(3*trsgc+2*mc2*(64
     &   *trs+13*trg+9*trc)+57*ms2*trs-9*ms2*trg)*gDs**4+2*ms2*(mc2*(48*
     &   trsgc+ms2*(123*trs+182*trg+44*trc))+3*ms2*(4*trsgc+ms2*(-7*trs+
     &   16*trg+4*trc))-5*mc2**2*(-3*trs-4*trg+2*trc))*gDs**3+ms2*(mc2**
     &   2*(36*trsgc+ms2*(37*trs+142*trg-10*trc))+2*mc2*ms2*(15*trsgc+ms
     &   2*(-4*trs+93*trg+23*trc))+6*ms2**3*(2*(trg+trc)-trs))*gDs**2-36
     &   *cDs**2*gDs**2*(2*gDs+ms2)*(2*trs*gDs+trsgc-mc2*trs+ms2*trg)+mc
     &   2*ms2**2*(36*mc2*trsgc+ms2*(-27*trsgc+11*mc2*trs+61*mc2*trg)-3*
     &   ms2**2*(trs-6*trg+2*trc))*gDs+mc2*ms2**3*(-9*ms2*trsgc+9*mc2*tr
     &   sgc-8*ms2**2*trc))+cDg**2*gDs**2*(2*gDs+ms2)*(72*(trg-2*trs)*gD
     &   s**5+cDs*(144*(trg-2*trs)*gDs**4-72*(mc2*(trs-trg)+ms2*trg)*gDs
     &   **3-4*(mc2*(9*trsgc-5*ms2*trs+15*ms2*trg+6*ms2*trc)+3*ms2**2*(2
     &   *trs+trg+trc))*gDs**2-mc2*ms2*(-54*trsgc+9*ms2*(trs-2*trg+2*trc
     &   )+19*mc2*trs+18*mc2*trg)*gDs+9*mc2*ms2*(mc2*(trsgc+2*ms2*trg)+m
     &   s2*trsgc))-36*(mc2*(trs-trg)+5*ms2*trg)*gDs**4+2*(mc2*(ms2*(174
     &   *trs-113*trg+22*trc)-9*trsgc)-6*ms2*(6*trsgc+ms2*(2*trs+trg+trc
     &   )))*gDs**3+18*cDs**2*gDs*(4*(trg-2*trs)*gDs**2-2*(mc2*(trs-trg)
     &   +ms2*trg)*gDs-mc2*(trsgc+ms2*trg))+ms2*(-mc2*(60*trsgc+ms2*(3*t
     &   rs-22*trg+32*trc))+mc2**2*(145*trs-12*trg+16*trc)+6*ms2**2*(-tr
     &   s+trg+trc))*gDs**2+mc2*ms2*(18*ms2*trsgc+mc2*(ms2*(11*trs+57*tr
     &   g-22*trc)-9*trsgc)+18*mc2**2*trs+ms2**2*(6*(trg+trc)-3*trs))*gD
     &   s+9*mc2**2*ms2**2*(2*trsgc+(ms2+mc2)*trg))+cDg*gDs**3*(2*gDs+ms
     &   2)*(72*trg*gDs**5+cDs*(144*trg*gDs**4+72*(trsgc+2*mc2*(trg-trs)
     &   )*gDs**3+12*(mc-ms)*(ms+mc)*(ms2*(trs+trg+trc)+3*mc2*trg)*gDs**
     &   2+2*mc2*ms2*(-9*trsgc+mc2*(5*(trs-3*trg)+3*trc)-3*ms2*(2*trs+tr
     &   g+trc))*gDs+9*mc2**2*ms2*(ms2-mc2)*trg)+36*(trsgc+2*mc2*(trg-tr
     &   s))*gDs**4+6*(mc2*ms2*(2*trs-25*trg+2*trc)-2*ms2**2*(trs+trg+tr
     &   c)+3*mc2**2*trg)*gDs**3+18*cDs**2*gDs*(4*trg*gDs**2+2*(trsgc+2*
     &   mc2*(trg-trs))*gDs+mc2*(mc-ms)*(ms+mc)*trg)-2*ms2*(-3*mc2*(ms2*
     &   (trg+trc)-15*trsgc)-2*mc2**2*(29*trs-28*trg+8*trc)+3*ms2**2*(tr
     &   s+trg+trc))*gDs**2+mc2*ms2*(-mc2*(54*trsgc+ms2*(-6*trs+25*trg+1
     &   6*trc))+mc2**2*(25*trs-22*trg+5*trc)+3*ms2**2*(-trs+trg+trc))*g
     &   Ds+9*mc2**2*ms2*(ms2-mc2)*trsgc)-ms2*cDg**5*(16*(16*trs+trg-18*
     &   trc)*gDs**4+4*(36*trsgc+ms2*(25*trs+4*trg-160*trc)+mc2*(trs+40*
     &   trg-18*trc))*gDs**3+2*cDs*(2*gDs+ms2)*(4*(9*trg+trc)*gDs**2+2*(
     &   9*trsgc+ms2*(8*trs+26*trg-7*trc))*gDs+ms2*(9*trsgc+9*ms2*trg-8*
     &   ms2*trc))+2*(3*mc2*(6*trsgc+ms2*(-5*trs+38*trg-24*trc))+144*ms2
     &   *trsgc+ms2**2*(-17*trs+10*trg-180*trc))*gDs**2+36*ms2*(4*ms2+mc
     &   2)*trsgc*gDs-2*ms2**2*(mc2*(16*trs-56*trg+45*trc)+ms2*(8*trs-tr
     &   g+12*trc))*gDs+9*ms2**2*(2*ms2+mc2)*trsgc+ms2**3*(mc2*(-8*trs+1
     &   9*trg-18*trc)+16*ms2*trc))+mc2*gDs**4*(2*gDs+ms2)*(36*trg*gDs**
     &   4+cDs*(72*trg*gDs**3+36*(trsgc+mc2*trg)*gDs**2+6*(mc-ms)*(ms+mc
     &   )*ms2*(trs+trg+trc)*gDs-9*mc2*ms2*(trsgc+mc2*trg))+18*(trsgc+mc
     &   2*trg)*gDs**3-6*ms2*(ms2*(trs+trg+trc)-mc2*(trs-5*trg+trc))*gDs
     &   **2+18*cDs**2*gDs*(2*trg*gDs+trsgc+mc2*trg)+ms2*(3*mc2*(2*ms2*(
     &   trs+trg+trc)-9*trsgc)+mc2**2*(7*trs-22*trg+5*trc)-3*ms2**2*(trs
     &   +trg+trc))*gDs-9*mc2**2*ms2*(trsgc+ms2*trg))+2*ms2*cDg**6*(2*gD
     &   s+ms2)*(4*(8*trs-trg+9*trc)*gDs**2+2*(ms2*(8*trs-10*trg+27*trc)
     &   -9*trsgc)*gDs+ms2*(ms2*(8*trs-trg+18*trc)-9*trsgc)))/(ms2*cDg**
     &   2*(2*cDg+mc2)*gDs**3*(gDs+cDg)*(2*gDs+ms2))/1.2e+1_dp+Bmmp

      Bmmp = tr3s00ft*(cDg*(8*trs*gDs**4+2*(8*trs*cDs+trsgc+5*ms2*trg)*
     &   gDs**3+2*(4*trs*cDs**2+2*(trsgc+ms2*trg)*cDs+ms2*(2*trsgc-5*mc2
     &   *trs+3*mc2*trg))*gDs**2+ms2*(2*trg*cDs**2+mc2*(-4*trs*cDs+2*trg
     &   *cDs-2*mc2*trs+5*ms2*trg+6*ms2*trc))*gDs+2*trsgc*cDs**2*gDs+mc2
     &   *ms2*(-trsgc*(cDs+ms2)-ms2*trg*(cDs+mc2)))+gDs*(-4*trg*gDs**4+c
     &   Ds*(-8*trg*gDs**3-4*(trsgc+mc2*trg)*gDs**2+mc2*ms2*(trsgc+mc2*t
     &   rg))-2*(trsgc+mc2*trg)*gDs**3+4*mc2*ms2*trg*gDs**2-2*cDs**2*gDs
     &   *(2*trg*gDs+trsgc+mc2*trg)+mc2*ms2*(3*trsgc+mc2*trs+3*mc2*trg)*
     &   gDs+mc2**2*ms2*(trsgc+ms2*trg))-cDg**2*(4*trs*gDs**3+8*trs*(cDs
     &   +2*ms2)*gDs**2+4*(trs*cDs**2+ms2*trsgc+ms2**2*trg)*gDs+mc2*ms2*
     &   (ms2*trg-2*trs*cDs))+8*ms2*trs*cDg**3*gDs)/(ms2*cDg**2*gDs)+Bm
     &   mp

      Bmmp = B0cgsf*(cDg**4*gDs*(4*trs*gDs**4+cDs*(8*trs*gDs**3+(-2*trs
     &   gc+8*ms2*trs-26*ms2*trg+4*ms2*trc)*gDs**2+ms2*(-6*trsgc-5*ms2*t
     &   rs-2*mc2*trs-8*ms2*trg+2*mc2*trg+2*ms2*trc)*gDs+ms2**2*(-trsgc-
     &   2*ms2*trg+mc2*trg+2*ms2*trc))+(ms2*(28*trs-2*trg+24*trc)-2*trsg
     &   c)*gDs**3+4*trs*cDs**2*gDs**2+ms2*(-14*trsgc-7*ms2*trs-12*ms2*t
     &   rg+4*mc2*(trc-trg)+20*ms2*trc)*gDs**2+ms2*(mc2*trsgc+2*ms2*(mc2
     &   *(trs+trg)-3*trsgc)+4*ms2**2*trc)*gDs-ms2**2*(2*ms2*trsgc-mc2*t
     &   rsgc+mc2*ms2*trg+2*ms2**2*trc))+cDg**5*(4*trs*gDs**4+4*(trs*cDs
     &   +4*ms2*(trs+trc))*gDs**3+ms2*(ms2*(-trs-8*trg+18*trc)-2*(2*trg*
     &   cDs+4*trsgc+mc2*(trs+trg)))*gDs**2+ms2*(2*(ms2*(-trs-3*trg+trc)
     &   -trsgc)*cDs-ms2*(4*trsgc+ms2*(-2*trs+trg-5*trc)+3*mc2*trg))*gDs
     &   -ms2**2*((trsgc+ms2*(trg-trc))*cDs+ms2*(trsgc+mc2*trg+ms2*trc))
     &   )+cDg*gDs**4*(4*(3*trg-2*trs)*gDs**4+cDs*(4*(5*trg-4*trs)*gDs**
     &   3+2*(3*trsgc-2*ms2*trg+5*mc2*trg)*gDs**2+2*ms2*(mc2*(3*trs+trg+
     &   2*trc)-ms2*(3*trs+2*trg+2*trc))*gDs-mc2*ms2*(trsgc-ms2*trg+2*mc
     &   2*trg))+(4*trsgc-6*ms**2*trg+6*mc**2*trg)*gDs**3-2*ms2*(trsgc+3
     &   *ms2*trs+mc2*(-8*trs+3*trg-3*trc)+2*ms2*trg+2*ms2*trc)*gDs**2+2
     &   *cDs**2*gDs*(4*(trg-trs)*gDs+trsgc-ms2*trg+2*mc2*trg)+ms2*(mc2*
     &   (-4*trsgc+4*ms2*trs+ms2*trg)+mc2**2*(2*trs-7*trg)-2*ms2**2*trs)
     &   *gDs+mc2*ms2*(ms2*trsgc-mc2*(2*trsgc+ms2*trg)))+cDg**3*gDs**2*(
     &   4*(trg-3*trs)*gDs**4+cDs*(4*(trg-3*trs)*gDs**3+2*(-3*trsgc+2*ms
     &   2*(2*trs-9*trg+2*trc)+mc2*trg)*gDs**2-2*ms2*(5*trsgc+ms2*(3*trs
     &   +7*trg+trc)+mc2*trs+5*mc2*trg)*gDs+ms2*(ms2*(trsgc+3*mc2*trg)+m
     &   c2*trsgc+ms2**2*(trc-trg)))+2*(-2*trsgc+ms2*(12*trs-5*trg+8*trc
     &   )+mc2*trg)*gDs**3+ms2*(-14*trsgc+2*mc2*(7*trs+5*trc)+ms2*(-9*tr
     &   s-4*trg+10*trc))*gDs**2-2*(trsgc+9*ms2*trg)*cDs**2*gDs+ms2*(2*m
     &   c2*(trsgc+2*ms2*(trs+trg-trc))-6*ms2*trsgc+mc2**2*(2*trs-3*trg)
     &   +ms2**2*(-2*trs+trg+3*trc))*gDs+ms2**2*(-ms2*trsgc+3*mc2*trsgc+
     &   mc2*ms2*trg+mc2**2*trg-ms2**2*trc))+cDg**2*gDs**3*(4*(3*trg-5*t
     &   rs)*gDs**4+cDs*(16*(trg-2*trs)*gDs**3+2*(-trsgc-9*ms2*trg+4*mc2
     &   *trg+2*ms2*trc)*gDs**2-ms2*(6*trsgc-2*mc2*(2*trs-2*trg+trc)+7*m
     &   s2*trs+6*ms2*trg+4*ms2*trc)*gDs+ms2*(mc2*(trsgc+3*ms2*trg)+ms2*
     &   trsgc-mc2**2*trg))+2*(ms2*(4*trs-7*trg+2*trc)+3*mc2*trg)*gDs**3
     &   +ms2*(-8*trsgc+2*mc2*(13*trs-trg+5*trc)-9*ms2*trs)*gDs**2+2*cDs
     &   **2*gDs*(trg*(2*gDs-6*ms2+mc2)-6*trs*gDs-trsgc)+ms2*(-3*ms2*trs
     &   gc+mc2**2*(4*trs-5*trg)+2*mc2*ms2*(2*trs+trg-3*trc)+ms2**2*(3*(
     &   trg+trc)-2*trs))*gDs+mc2*ms2*(3*ms2*trsgc-mc2*trsgc+ms2**2*trg+
     &   mc2*ms2*trg))+gDs**5*(4*trg*gDs**4+cDs*(8*trg*gDs**3+4*(trsgc+m
     &   c2*trg)*gDs**2+2*(mc-ms)*(ms+mc)*ms2*(trs+trg+trc)*gDs-mc2*ms2*
     &   (trsgc+mc2*trg))+2*(trsgc+mc2*trg)*gDs**3-2*ms2*(ms2*(trs+trg+t
     &   rc)-mc2*(trs-trg+trc))*gDs**2+2*cDs**2*gDs*(2*trg*gDs+trsgc+mc2
     &   *trg)-ms2*(mc2*(3*trsgc-2*ms2*(trs+trg+trc))+ms2**2*(trs+trg+tr
     &   c)+3*mc2**2*trg)*gDs-mc2**2*ms2*(trsgc+ms2*trg))+ms2*cDg**6*(4*
     &   (trs+trc)*gDs**2-2*trsgc*gDs+2*ms2*(trs-trg+3*trc)*gDs-ms2*trsg
     &   c+ms2**2*(trs+2*trc)))/(ms2*cDg**2*gDs**3*(gDs+cDg)**2)/4.0_dp+
     &   Bmmp

      Bmmp = lc*(cDg*(8*trg*gDs**4+4*(trsgc+mc2*(trg-2*trs))*gDs**3+2*c
     &   Ds*gDs*(4*trg*gDs**2+2*(trsgc+mc2*(trg-2*trs))*gDs-mc2*(trsgc+m
     &   s2*trg))-2*mc2*(trsgc+ms2*trg)*gDs**2+2*mc2*ms2*(trsgc-mc2*(trs
     &   +2*trg))*gDs+mc2**2*ms2*(trsgc+ms2*trg))-2*cDg**2*(8*trs*gDs**3
     &   +2*(4*trs*cDs+trsgc-mc2*trs+ms2*trg)*gDs**2+2*((trsgc+ms2*trg)*
     &   cDs-mc2*(trs*cDs-ms2*trs+ms2*trc))*gDs+mc2*ms2*(-trsgc+mc2*trs-
     &   ms2*trg))+mc2*gDs*(trg*(2*gDs*(gDs+cDs)*(2*gDs+mc2)-mc2**2*ms2)
     &   +trsgc*(2*gDs*(gDs+cDs)+mc2*ms2))+2*trs*cDg**3*(4*gDs*(gDs+cDs)
     &   -3*mc2*ms2))/(ms2*(2*cDg+mc2)**2*gDs)/2.0_dp+Bmmp

      Bmmp = BfunX*(cDg**3*(4*(3*trs+trg)*gDs**3+cDs*(4*(3*trs+2*trg+tr
     &   c)*gDs**2+2*ms2*(trs+trg+2*trc)*gDs+ms2**2*trc)+ms2*(7*trs+4*tr
     &   g+2*trc)*gDs**2+3*ms2**2*trc*gDs+ms2**3*trc)+cDg**4*((6*trs+4*t
     &   rg)*gDs**2-2*ms2*(trs+trg)*gDs-ms2**2*(trs+trg))+cDg**2*gDs**2*
     &   (trs*(ms2*(5*gDs+3*cDs)+2*(gDs+cDs)*(5*gDs+cDs)+ms2**2)-4*ms2*(
     &   trg+trc)*gDs)+(trs+trg+trc)*gDs**4*(2*(gDs+cDs)+ms2)**2+2*cDg*g
     &   Ds**3*((trg+trc)*(gDs-ms2)+2*trs*gDs+trs*cDs)*(2*(gDs+cDs)+ms2)
     &   )/(cDg**2*gDs**3)/4.0_dp+Bmmp

      Bmmp = B0csf*(cDg**2*(cDs*((4*mc2*trs-2*ms2*trs)*gDs**2+(2*ms2*(t
     &   rsgc+ms2*trc)+mc2*(ms2*(trs-3*trg+4*trc)-trsgc)+mc2**2*trs)*gDs
     &   +2*mc2*ms2**2*trg)+mc2*ms2*(2*trs*gDs**2+(trsgc-mc2*trg+ms2*trc
     &   +mc2*trc)*gDs+ms2*(ms2+mc2)*trg)-2*cDs**3*(2*(trs-trg+trc)*gDs+
     &   ms2*trg)-cDs**2*(2*(mc2*trc-ms2*(-trs+trg+trc))*gDs+ms2*(ms2+mc
     &   2)*trg))+cDg*gDs*(cDs*(2*mc2*trs*gDs**2+(ms2*(trsgc+ms2*trc)+mc
     &   2*(-2*trsgc+4*ms2*trs-5*ms2*trg+5*ms2*trc)+2*mc2**2*trs)*gDs-2*
     &   mc2*ms2*(trsgc+2*ms2*trg+3*mc2*trg))-mc2*ms2*(-2*trs*gDs**2+(tr
     &   sgc-mc2*(trs+2*trc)+ms2*(-trs+trg-2*trc))*gDs+(ms2+mc2)*(trsgc+
     &   mc2*trg))+cDs**2*(-2*(ms2*(trs-trg+trc)+mc2*(-trs+trg+2*trc))*g
     &   Ds+mc2*(trsgc-7*ms2*trg)+ms2*trsgc+mc2**2*trg)+2*cDs**3*(-2*(tr
     &   s-trg+2*trc)*gDs+trsgc+2*ms2*trg+3*mc2*trg)+8*trg*cDs**4)+cDg**
     &   3*(cDs*(mc2*(2*trs*gDs+ms2*(-2*trs+trg+trc))+ms2*(-4*trs*gDs+tr
     &   sgc+ms2*trc))+mc2*ms2*(-trs*(2*gDs+mc2)+trsgc+ms2*(trg-trs))+2*
     &   ms2*trc*cDs**2)+gDs**2*(mc2*cDs*((-trsgc+ms2*(trs-trg+2*trc)+mc
     &   2*trs)*gDs-2*ms2*(trsgc+(ms2+mc2)*trg))+cDs**2*(-2*(mc2*(-trs+t
     &   rg+trc)+ms2*trc)*gDs+mc2*(trsgc-4*ms2*trg)+ms2*trsgc)+mc2*ms2*(
     &   mc2*trg*gDs+(ms2+mc2)*trc*gDs+trsgc*(-gDs-ms2-mc2))+2*cDs**3*(-
     &   2*trc*gDs+trsgc+(ms2+mc2)*trg)+4*trg*cDs**4)-2*ms2*trs*cDg**4*(
     &   cDs+mc2))/((cDs-mc*ms)*(cDs+mc*ms)*gDs*(gDs+cDg)**2)/2.0_dp+Bm
     &   mp

      Bmmp = 3*mc2*tr3s002ft*((trg+trc)*((2*cDg+mc2)*gDs**2+ms2*cDg**2)
     &   -trs*cDg**2*(gDs+cDs))/(cDg**2*gDs)+Bmmp

      Bmmp = ls*(gDs*(2*gDs+ms2)*((2*gDs+ms2)*(2*trc*gDs-trg*(4*cDs+mc2
     &   ))-ms2*trsgc)+cDg*(-trsgc*(2*gDs+ms2)**2+8*trc*gDs**3-4*ms2*trs
     &   *gDs**2+8*ms2*trc*gDs**2+4*ms2**2*trg*gDs+2*ms2**2*trc*gDs+ms2*
     &   *3*trg))/(gDs*(2*gDs+ms2)**2)/2.0_dp+Bmmp

      Bmmp = 3*mc2*tr3s001ft*(mc2*trs*gDs-2*trs*cDg*cDs+2*ms2*(trg+trc)
     &   *cDg)/cDg**2+Bmmp

      Bmmp=Bmmp/zmmp

      return
      end
