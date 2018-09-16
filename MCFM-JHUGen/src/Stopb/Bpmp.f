      subroutine Bamp_pmp(q,mc,ms,Bpmp)
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
      complex(dp):: trc,trg,trs,trsgc,zpmp,Bpmp

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
      zpmp=za(2,4)**2*zb(4,3)

      Bpmp = -ms2*(cDg**3*gDs*(cDs*(8*trc*gDs**2-2*(2*trsgc+mc2*trg+3*m
     &   s2*trc)*gDs+ms2*(5*trsgc+mc2*trg))+mc2*(-4*trc*gDs**2+(2*trsgc-
     &   2*ms2*trc)*gDs+ms2*(trsgc-ms2*trg))+2*cDs**2*(2*trc*gDs-2*trsgc
     &   +ms2*trg))+cDg**2*gDs**2*(cDs*(4*trc*gDs**2-2*(trsgc+2*mc2*trg+
     &   3*ms2*trc)*gDs+4*ms2*trsgc+mc2*trsgc+5*mc2*ms2*trg)+mc2*(-8*trc
     &   *gDs**2+(4*trsgc+2*mc2*trg-6*ms2*trc)*gDs+ms2*(3*trsgc-2*ms2*tr
     &   g+mc2*trg))+2*cDs**2*(6*trc*gDs-5*trsgc+2*ms2*trg-mc2*trg)-4*tr
     &   g*cDs**3)+cDg*gDs**3*(mc2*(-4*trc*gDs**2+2*(trsgc+2*mc2*trg-ms2
     &   *trc)*gDs+ms2*(trsgc+ms2*trg+2*mc2*trg))+cDs*(-2*(mc2*trg+ms2*t
     &   rc)*gDs+ms2*trsgc+2*mc2*trsgc+7*mc2*ms2*trg)+cDs**2*(8*trc*gDs-
     &   6*trsgc-4*mc2*trg)-8*trg*cDs**3)+2*cDg**4*cDs*(2*trc*gDs**2-(tr
     &   sgc+ms2*trc)*gDs+ms2*trsgc)+mc2*gDs**4*(2*(mc2*trg+ms2*trc)*gDs
     &   +trsgc*(cDs-ms2)-ms2*trg*(cDs+mc2)))*tr5Xs/(cDg*(cDs-mc*ms)*(cD
     &   s+mc*ms)*gDs**2*(gDs+cDg))/2.0_dp

      Bpmp = mc2*(-2*mc2*trsgc*cDs*gDs**4-cDg**2*gDs**2*(cDs*(mc2*(-6*t
     &   rs*gDs+4*trsgc+ms2*trg)+ms2*trsgc)+mc2*ms2*(-2*trs*gDs+4*trc*gD
     &   s+trsgc+ms2*trg-4*mc2*trg)-2*cDs**2*(2*(trc-trs)*gDs+4*trsgc+ms
     &   2*trg-3*mc2*trg))-cDg**3*gDs*(cDs*(mc2*(-6*trs*gDs+trsgc+7*ms2*
     &   trg)+2*ms2*trsgc)+2*cDs**2*(6*trs*gDs-2*trc*gDs-3*trsgc-3*ms2*t
     &   rg+mc2*trg)+mc2*ms2*(-6*trs*gDs+4*trc*gDs+trsgc+4*ms2*trg-mc2*t
     &   rg)-8*trg*cDs**3)-cDg**4*(cDs*(-2*mc2*trs*gDs+ms2*trsgc+3*mc2*m
     &   s2*trg)-2*cDs**2*(-6*trs*gDs+trsgc+ms2*trg)+mc2*ms2*(-6*trs*gDs
     &   +trsgc+ms2*trg)-4*trg*cDs**3)-cDg*gDs**3*(mc2*cDs*(-2*trs*gDs+5
     &   *trsgc+ms2*trg)+2*(mc2*trg-2*trsgc)*cDs**2+mc2*ms2*(trsgc-mc2*t
     &   rg))-2*trs*cDg**5*(2*cDs**2-mc2*ms2))*tr5Xc/(cDg**3*(cDs**2-mc2
     &   *ms2)*(gDs+cDg))/2.0_dp+Bpmp

      Bpmp = ms2*(trsgc*(ms2*(cDg*cDs*(gDs+2*cDg)-mc2*gDs*(gDs+cDg))+gD
     &   s*(mc2*(cDs+2*cDg)*gDs-2*cDg*cDs*(cDs+cDg)))+mc2*trg*gDs*(gDs*(
     &   2*mc2*gDs-ms2*(cDs+mc2))+cDg*(cDs*(ms2-2*gDs)+ms2**2))-2*trc*gD
     &   s*(mc2*gDs-cDg*cDs)*(2*cDg*gDs-ms2*(gDs+cDg)))*tr4Xs/(cDg*(cDs-
     &   mc*ms)*(cDs+mc*ms)*gDs)/2.0_dp+Bpmp

      Bpmp = mc2*(trsgc*(ms2*cDg*(mc2*gDs+cDg*(cDs+mc2))-cDs*(2*gDs+cDg
     &   )*(2*cDg*cDs-mc2*gDs))+trg*cDg*(mc2*(2*cDs**2+ms2*cDs-mc2*ms2)*
     &   gDs+cDg*(-4*cDs**3-2*ms2*cDs**2+3*mc2*ms2*cDs+mc2*ms2**2))+2*tr
     &   s*cDg*(gDs+cDg)*(-mc2*cDs*gDs+2*cDg*cDs**2-mc2*ms2*cDg))*tr4Xc/
     &   (cDg**2*(cDs-mc*ms)*(cDs+mc*ms))/2.0_dp+Bpmp

      Bpmp = (-cDg**2*(mc2*ms2*(-4*(trg-5*trs)*gDs**2+2*(trsgc+mc2*trs+
     &   ms2*trg)*gDs+ms2*(mc2*trg-trsgc))+mc2*ms2*cDs*(-2*trs*gDs+trsgc
     &   +ms2*trg)+2*cDs**2*gDs*(-10*trs*gDs+2*trg*gDs-trsgc-ms2*trg))+c
     &   Dg*gDs*(mc2*ms2*(-4*(2*trg+3*trc)*gDs**2+2*(trsgc-4*mc2*trs+mc2
     &   *trg)*gDs+mc2*(trsgc+ms2*trg))+2*cDs**2*gDs*(4*trg*gDs+6*trc*gD
     &   s-trsgc+3*mc2*trs-mc2*trg)+mc2*ms2*(mc2*trg-3*trsgc)*cDs)+2*mc2
     &   **2*ms2*trsgc*gDs**2+2*cDg**3*(2*mc2*ms2*(trs*gDs-ms2*trg)+2*cD
     &   s**2*(ms2*trg-trs*gDs)+mc2*ms2*trs*cDs))*tr3Xs/(cDg**2*(cDs**2-
     &   mc2*ms2))/2.0_dp+Bpmp

      Bpmp = mc2*(trg*gDs*(cDs*(2*mc2*gDs**2-mc2*ms2*gDs+ms2**2*cDg)-mc
     &   2*ms2*(2*cDg+ms2)*gDs+ms2*cDg*cDs**2)+trsgc*(ms2*(mc2-cDs)*gDs*
     &   *2+cDg*gDs*(2*cDs*gDs-3*ms2*cDs+ms2**2)+2*ms2*cDg**2*(ms2-gDs))
     &   -2*trc*gDs*(cDs*gDs-ms2*cDg)*(2*cDg*gDs-ms2*(gDs+cDg)))*tr3Xc/(
     &   cDg*(cDs-mc*ms)*(cDs+mc*ms)*gDs)/2.0_dp+Bpmp

      Bpmp = (cDg**3*(cDs**2*(4*(trs+2*trc)*gDs**2+(ms2*(2*trs-4*trg-6*
     &   sck*trc+8*trc)-4*trsgc)*gDs-ms2*(-3*sck*trsgc+2*trsgc+2*ms2*trs
     &   +3*ms2*trg))+mc2*ms2*(-4*(trs+2*trc)*gDs**2+4*trsgc*gDs+2*ms2*(
     &   -2*trs+trg+3*(sck-1)*trc)*gDs+ms2*(-3*sck*trsgc+2*trsgc+2*ms2*t
     &   rs+3*ms2*trg))+ms2*cDs*(4*trc*gDs**2+2*(-trsgc-mc2*trs+ms2*trc)
     &   *gDs+mc2*ms2*trg)-ms2*trg*cDs**3)+cDg**2*gDs*(mc2*ms2*(-4*(-5*t
     &   rs+trg+trc)*gDs**2+2*(2*trsgc+mc2*(trs+2*trg)+ms2*trg)*gDs+ms2*
     &   (trsgc+mc2*(4-3*sck)*trg))+cDs**2*(4*(trg-5*trs)*gDs**2-2*(trsg
     &   c+ms2*(trg+trc)+2*mc2*trg)*gDs+ms2*(mc2*(3*sck-4)*trg-trsgc))+2
     &   *mc2*ms2*cDs*((trs+trg-trc)*gDs+2*trsgc)-4*trsgc*cDs**3)+ms2*(c
     &   Ds-mc*ms)*(cDs+mc*ms)*gDs**2*((trg+trc)*(cDs+ms2)*(2*gDs+2*cDs+
     &   ms2-mc2)-trs*(4*gDs**2+cDs*(6*gDs+3*ms2-mc2)+4*ms2*gDs-2*mc2*gD
     &   s+2*cDs**2+ms2**2-mc2*ms2+mc2**2))-cDg*(cDs-mc*ms)*(cDs+mc*ms)*
     &   gDs**2*(4*(2*trg+3*trc)*gDs**2-2*(trsgc+ms2*(-trs+trg+trc)+mc2*
     &   (trg-3*trs))*gDs+ms2*(4*trs*cDs+mc2*(trs+trg+trc)-ms2*(-trs+trg
     &   +trc))-2*mc2*trsgc)+ms2*trg*cDg**4*(mc*ms-cDs)*(cDs+mc*ms))*lVs
     &   /(cDg**3*(cDs-mc*ms)*(cDs+mc*ms)*gDs)/4.0_dp+Bpmp

      Bpmp = (cDg**3*(cDs**2*(4*trs*gDs**3+4*ms2*(trs-trg+trc)*gDs**2-2
     &   *ms2*(2*trsgc+mc2*trs+ms2*trg)*gDs+ms2**2*(2*trsgc+2*mc2*trs+mc
     &   2*trg))-mc2*ms2*(4*trs*gDs**3+4*ms2*(trs-trg)*gDs**2-2*ms2*(trs
     &   gc+ms2*(trg+trc))*gDs+ms2**2*(2*trsgc+mc2*(2*trs+trg)))+4*cDs**
     &   3*gDs*(trs*gDs-2*ms2*trg)+2*mc2*ms2*cDs*gDs*(ms2*(-trs+3*trg+tr
     &   c)-2*trs*gDs))+cDg**2*gDs*(cDs**2*(4*(trg-5*trs)*gDs**3-2*(trsg
     &   c+ms2*trg)*gDs**2-2*ms2*(trsgc+mc2*(trc-trs))*gDs-mc2*ms2*(trsg
     &   c+ms2*trg))+mc2*ms2*(-4*(trg-5*trs)*gDs**3+2*(trsgc+ms2*trg)*gD
     &   s**2+2*ms2*(trsgc+mc2*trg)*gDs+mc2*ms2*(trsgc+ms2*trg))+cDs**3*
     &   (4*(trg-5*trs)*gDs**2-2*(trsgc+ms2*trg)*gDs-4*ms2*trsgc)+2*mc2*
     &   ms2*cDs*(gDs*(-2*(-5*trs+trg+trc)*gDs+mc2*trs+ms2*(trg-trc))+2*
     &   trsgc*(gDs+ms2)))+cDg*(cDs-mc*ms)*(cDs+mc*ms)*gDs**2*(-4*(2*trg
     &   +3*trc)*gDs**3+2*cDs*(-2*(2*trg+3*trc)*gDs**2+(trsgc-3*mc2*trs+
     &   mc2*trg)*gDs+2*mc2*ms2*trs)+2*(trsgc+mc2*(trg-3*trs))*gDs**2-2*
     &   mc2*ms2*(-trs+trg+trc)*gDs+mc2*ms2*(trsgc+mc2*trs-ms2*(-trs+trg
     &   +trc)-mc2*trg))-mc2*ms2*(cDs-mc*ms)*(cDs+mc*ms)*gDs**2*((trg+tr
     &   c)*(cDs+ms2)*(2*gDs+2*cDs+ms2-mc2)-trs*(4*gDs**2+cDs*(6*gDs+3*m
     &   s2-mc2)+4*ms2*gDs-2*mc2*gDs+2*cDs**2+ms2**2-mc2*ms2+mc2**2))+4*
     &   ms2*trs*cDg**4*(cDs-mc*ms)*(cDs+mc*ms)*gDs)*lVc/(ms2*cDg**3*(cD
     &   s-mc*ms)*(cDs+mc*ms)*gDs)/4.0_dp+Bpmp

      Bpmp = xs*cDs*(-cDg**2*(4*trc*gDs**2-2*trsgc*gDs+ms2*trsgc)-mc2*t
     &   rsgc*gDs**2+2*cDg*gDs*(mc2*trg*gDs+trsgc*cDs))*lRcs/(mc*ms*(xs*
     &   *2-1)*cDg**2*gDs)+Bpmp

      Bpmp = (cDg**3*gDs*(-72*trs*gDs**5+12*ms2*(-63*trs+18*trg+4*trc)*
     &   gDs**4+cDs*(2*gDs+ms2)*(-72*trs*gDs**3+4*ms2*(10*trs+23*trg+8*t
     &   rc)*gDs**2+18*ms2*(-3*trsgc+mc2*trs+2*ms2*trc)*gDs+9*ms2**2*(ms
     &   2-mc2)*trg)-4*ms2*(15*trsgc+ms2*(82*trs+8*trg+3*(5*sck-4)*trc)+
     &   17*mc2*trs+6*mc2*trg)*gDs**3+9*cDs**2*(2*gDs+ms2)*(ms2**2*trg-4
     &   *trs*gDs**2)-2*ms2**2*(3*(7-5*sck)*trsgc+mc2*(7*trs-15*(sck-6)*
     &   trg)+ms2*(-11*trs+32*trg+(15*sck+28)*trc))*gDs**2+3*ms2**2*(ms2
     &   *(5*sck+4)-6*mc2)*trsgc*gDs+ms2**3*(ms2*(6*(trs+trg)-34*trc)+mc
     &   2*(10*trs+3*(5*sck-34)*trg))*gDs-9*ms2**3*(mc2*(trsgc+ms2*trg)-
     &   ms2*trsgc))+cDg*gDs**3*(2*gDs+ms2)*(36*(5*trs+trg+3*trc)*gDs**4
     &   +cDs*(72*(5*trs+trg+3*trc)*gDs**3+36*(3*mc2*trs+ms2*trg-mc2*trg
     &   )*gDs**2+2*ms2*(9*trsgc-3*ms2*(-5*trs+trg+trc)+mc2*(-22*trs+37*
     &   trg+19*trc))*gDs+ms2*(-9*ms2**2*(-trs+trg+trc)+3*mc2*ms2*(-4*tr
     &   s+trg+4*trc)+mc2**2*(-5*trs+14*trg+5*trc)))-18*(-3*mc2*trs+7*ms
     &   2*trg+mc2*trg+12*ms2*trc)*gDs**3+6*cDs**2*(6*(5*trs+trg+3*trc)*
     &   gDs**2+3*(3*mc2*trs+ms2*trg-mc2*trg)*gDs+(mc-ms)*(ms+mc)*ms2*(-
     &   trs+trg+trc))-6*ms2*(-9*trsgc+mc2*(60*trs-7*trg+17*trc)+ms2*(-3
     &   *trs+trg+trc))*gDs**2+ms2*(mc2*(54*trsgc+ms2*(-18*trs+5*trg-4*t
     &   rc))-3*ms2**2*(3*(trg+trc)-5*trs)+mc2**2*(27*trg-79*trs))*gDs+m
     &   s2*(mc2*(6*ms2**2*(-trs+trg+trc)-9*ms2*trsgc)+mc2**2*(9*trsgc+m
     &   s2*(6*trs-11*(trg+trc)))+5*mc2**3*trs-3*ms2**3*(-trs+trg+trc)))
     &   +ms2*cDg**4*(24*(3*trs+2*trc)*gDs**4-4*(6*(trsgc-3*trs*cDs)+ms2
     &   *(-9*trs+31*trg+(15*sck-11)*trc))*gDs**3+2*ms2*(18*(trs+2*trg+t
     &   rc)*cDs+15*(sck+1)*trsgc+mc2*(10*trs+trg)-ms2*(5*trs+27*trg+3*(
     &   5*sck+13)*trc))*gDs**2+ms2*(18*trg*cDs**2+18*(trsgc+ms2*(3*trg+
     &   trc))*cDs+ms2*(15*sck*trsgc+39*trsgc+10*mc2*trs+ms2*(-2*trs+7*t
     &   rg-44*trc)+mc2*trg))*gDs+9*ms2**2*(cDs+ms2)*(trg*cDs+trsgc))-cD
     &   g**2*gDs**2*(2*gDs+ms2)*(36*(trg-4*trs)*gDs**4+cDs*(72*(trg-4*t
     &   rs)*gDs**3-4*(9*trsgc+ms2*(trs+26*trg+8*trc))*gDs**2-2*ms2*(-18
     &   *trsgc+6*ms2*trs+2*mc2*(-2*trs+17*trg+8*trc)-9*ms2*trg+9*ms2*tr
     &   c)*gDs+9*ms2*(mc2*(trsgc+2*ms2*trg)+ms2*trsgc))+18*(ms2*(22*trs
     &   +trg+12*trc)-trsgc)*gDs**3+2*ms2*(-18*trsgc+mc2*(191*trs-45*trg
     &   -3*trc)+3*ms2*(-trs+4*trg+trc))*gDs**2+18*cDs**2*gDs*(-8*trs*gD
     &   s+2*trg*gDs-trsgc-ms2*trg)+ms2*(3*ms2*(6*trsgc+ms2*(-trs+trg+tr
     &   c))+mc2*(ms2*(-15*sck*trg+86*trg+10*trc)-27*trsgc)+13*mc2**2*tr
     &   s)*gDs+9*mc2*ms2**2*(2*trsgc+(ms2+mc2)*trg))+gDs**4*(2*gDs+ms2)
     &   *(36*(2*trg+3*trc)*gDs**4+cDs*(72*(2*trg+3*trc)*gDs**3-36*(trsg
     &   c+mc2*(trg-3*trs))*gDs**2+6*(mc-ms)*(ms+mc)*ms2*(-3*trs+trg+trc
     &   )*gDs+ms2*(3*mc2*(3*trsgc+4*ms2*(-trs+trg+trc))-9*ms2**2*(-trs+
     &   trg+trc)+mc2**2*(-5*trs+14*trg+5*trc)))-18*(trsgc+mc2*(trg-3*tr
     &   s))*gDs**3+6*cDs**2*(6*(2*trg+3*trc)*gDs**2-3*(trsgc+mc2*(trg-3
     &   *trs))*gDs+(mc-ms)*(ms+mc)*ms2*(-trs+trg+trc))-12*ms2*(mc2*(trs
     &   +6*trg+9*trc)-ms2*trs)*gDs**2+3*ms2*(mc2*(9*trsgc+2*ms2*(-3*trs
     &   +trg+trc))-2*ms2**2*(-2*trs+trg+trc)+mc2**2*(9*trg-22*trs))*gDs
     &   +ms2*(mc2**2*(9*trsgc+ms2*(6*trs-2*trg-11*trc))+5*mc2**3*trs-3*
     &   ms2**3*(-trs+trg+trc)+6*mc2*ms2**2*(-trs+trg+trc)))-ms2**2*cDg*
     &   *5*(2*gDs+ms2)*((16*trs-11*trg+2*trc)*gDs-9*(2*trg*cDs+trsgc)+m
     &   s2*(8*trs-trg+10*trc))+9*ms2**2*trg*cDg**6*(2*gDs+ms2))/(ms2*cD
     &   g**3*gDs**2*(gDs+cDg)*(2*gDs+ms2))/1.2e+1_dp+Bpmp

      Bpmp = tr3s00ft*(-4*(2*trg+3*trc)*gDs**5-cDs*(8*(2*trg+3*trc)*gDs
     &   **4-4*(2*(trg-5*trs)*cDg+trsgc+mc2*(trg-3*trs))*gDs**3+4*cDg*(-
     &   2*trs*cDg+trsgc+ms2*trg)*gDs**2+ms2*(2*cDg+mc2)*(2*(trg-2*trs)*
     &   cDg+trsgc+mc2*(trg-3*trs))*gDs+ms2*cDg*((2*cDg+mc2)*(2*trs*cDg-
     &   ms2*trg)-trsgc*(4*cDg+mc2)))+2*(2*(trg-5*trs)*cDg+trsgc+mc2*(tr
     &   g-3*trs))*gDs**4+2*(2*trs*cDg**2-trsgc*cDg+7*ms2*trg*cDg+12*ms2
     &   *trc*cDg+4*mc2*ms2*trg+6*mc2*ms2*trc)*gDs**3-2*cDs**2*gDs*((4*t
     &   rg+6*trc)*gDs**2-(2*(trg-5*trs)*cDg+trsgc-3*mc2*trs+mc2*trg)*gD
     &   s+cDg*(-2*trs*cDg+trsgc+ms2*trg))-ms2*(2*cDg+mc2)*(-22*trs*cDg+
     &   6*trg*cDg+3*trsgc-8*mc2*trs+3*mc2*trg)*gDs**2+ms2*(trsgc*(2*cDg
     &   **2-2*mc2*cDg-mc2**2)+(2*cDg+mc2)*(mc2*(5*trs*cDg+2*ms2*trg+3*m
     &   s2*trc)+cDg*(ms2*trg-2*trs*cDg)))*gDs+ms2**2*cDg*(-2*(3*trs-2*t
     &   rg+trc)*cDg**2+(2*trsgc-3*mc2*trs+5*mc2*trg)*cDg+mc2*(trsgc+mc2
     &   *trg)))/(ms2*cDg**3)+Bpmp

      Bpmp = B0cgsf*(cDg*gDs**4*(4*(5*(trs+trg)+9*trc)*gDs**4+cDs*(4*(1
     &   0*trs+8*trg+15*trc)*gDs**3+2*(-3*trsgc-5*mc2*(trg-3*trs)+2*ms2*
     &   trg)*gDs**2+2*ms2*(trsgc-2*ms2*(-4*trs+trg+trc)+mc2*(-9*trs+6*t
     &   rg+4*trc))*gDs+ms2*(mc2*(trsgc+ms2*(-8*trs+7*trg+8*trc))-6*ms2*
     &   *2*(-trs+trg+trc)+2*mc2**2*trg))-2*(2*trsgc+3*mc2*(trg-3*trs)+3
     &   *ms2*(trg+2*trc))*gDs**3+2*cDs**2*(2*(5*trs+3*trg+6*trc)*gDs**2
     &   +(-trsgc+6*mc2*trs+ms2*trg-2*mc2*trg)*gDs+2*(mc-ms)*(ms+mc)*ms2
     &   *(-trs+trg+trc))-2*ms2*(-2*trsgc+mc2*(21*trs+2*trg+11*trc)+ms2*
     &   (-5*trs+trg+trc))*gDs**2+ms2*(mc2*(6*trsgc+ms2*(-12*trs+5*trg+4
     &   *trc))+ms2**2*(9*trs-5*(trg+trc))+7*mc2**2*(trg-2*trs))*gDs+mc2
     &   *(2*mc2-ms2)*ms2*trsgc-ms2**2*(2*ms2**2*(-trs+trg+trc)-4*mc2*ms
     &   2*(-trs+trg+trc)+mc2**2*(-4*trs+3*trg+4*trc)))+cDg**2*gDs**3*(4
     &   *(14*trs+3*trg+9*trc)*gDs**4+cDs*(4*(23*trs+3*trg+12*trc)*gDs**
     &   3+2*(trsgc+12*mc2*trs+9*ms2*trg-4*mc2*trg+2*ms2*trc)*gDs**2+2*m
     &   s2*(3*trsgc+7*ms2*trs+mc2*(-8*trs+7*trg+5*trc)-2*ms2*trg)*gDs+m
     &   s2*(-ms2*(trsgc+3*ms2*(-trs+trg+trc))+mc2*(ms2*(-4*trs+trg+4*tr
     &   c)-trsgc)+mc2**2*trg))-2*(ms2*(12*trs+trg+12*trc)+3*mc2*(trg-3*
     &   trs))*gDs**3+2*cDs**2*(6*(3*trs+trc)*gDs**2+(trsgc+3*mc2*trs+2*
     &   ms2*trg-mc2*trg)*gDs+(mc2-ms2)*ms2*(-trs+trg+trc))-2*ms2*(-5*tr
     &   sgc+mc2*(37*trs-9*trg+3*trc)+ms2*(-4*trs+2*trg+trc))*gDs**2-ms2
     &   *(ms2*(trsgc+ms2*(-6*trs+4*trg+4*trc))+mc2*(ms2*(6*trs+4*trg-2*
     &   trc)-4*trsgc)+mc2**2*(10*trs-3*trg))*gDs+ms2*(mc2**2*(trsgc+ms2
     &   *(2*trs-3*trg-2*trc))+mc2*ms2*(ms2*(-2*trs+trg+2*trc)-3*trsgc)-
     &   ms2**3*(-trs+trg+trc)))+cDg**3*gDs**2*(4*(12*trs-trg+3*trc)*gDs
     &   **4+cDs*(4*(15*trs-2*trg+3*trc)*gDs**3+2*(3*trsgc-2*ms2*trs+3*m
     &   c2*trs+20*ms2*trg-mc2*trg)*gDs**2+ms2*(10*trsgc+ms2*(4*trs+11*t
     &   rg+6*trc)+2*mc2*(-trs+8*trg+2*trc))*gDs-ms2*(ms2+mc2)*trsgc+ms2
     &   **2*(ms2-3*mc2)*trg)-2*(-2*trsgc+ms2*(26*trs-9*trg+10*trc)+mc2*
     &   (trg-3*trs))*gDs**3+cDs**2*(-4*(trg-3*trs)*gDs**2+2*(trsgc+9*ms
     &   2*trg)*gDs+ms2**2*trg)-2*ms2*(-7*trsgc+ms2*(trs+trg+trc)-mc2*(1
     &   0*(trg-2*trs)+trc))*gDs**2+ms2*(mc2*(4*trsgc+mc2*(trg-2*trs))+m
     &   s2**2*(3*trs+4*trg-5*trc)-9*mc2*ms2*trg)*gDs-ms2**2*(-ms2*trsgc
     &   +3*mc2*trsgc+2*mc2*ms2*trg+mc2**2*trg))+cDg**5*(-4*trs*gDs**4-4
     &   *(trs*cDs+ms2*(3*trs-trg+2*trc))*gDs**3-ms2*(4*trs*cDs-8*trg*cD
     &   s-6*trsgc+8*ms2*trs+2*mc2*trs+ms2*trg)*gDs**2+ms2**2*(7*trg*cDs
     &   +2*trc*cDs+2*trsgc+3*ms2*trg-6*ms2*trc)*gDs+ms2**2*(cDs+ms2)*(t
     &   rg*cDs+trsgc))+cDg**4*gDs*(-4*(trg-2*trs)*gDs**4+cDs*(-4*(trg-t
     &   rs)*gDs**3+2*(trsgc+ms2*(-4*trs+17*trg-2*trc))*gDs**2+2*ms2*(3*
     &   trsgc+mc2*trs+5*ms2*trg+2*mc2*trg+3*ms2*trc)*gDs+ms2**2*(trsgc+
     &   2*ms2*trg-mc2*trg))+2*(trsgc+ms2*(9*(trg-2*trs)-8*trc))*gDs**3+
     &   cDs**2*(-4*trs*gDs**2+8*ms2*trg*gDs+2*ms2**2*trg)-ms2*(-14*trsg
     &   c+ms2*(10*trs+3*trg+2*trc)+2*mc2*(3*trs-3*trg+trc))*gDs**2+ms2*
     &   (ms2*(2*trsgc-6*mc2*trg)+3*mc2*trsgc+ms2**2*(3*trs+6*trg-9*trc)
     &   )*gDs-ms2**2*(mc2*(trsgc+ms2*trg)-2*ms2*trsgc))+gDs**5*(4*(2*tr
     &   g+3*trc)*gDs**4+cDs*(8*(2*trg+3*trc)*gDs**3-4*(trsgc+mc2*(trg-3
     &   *trs))*gDs**2+2*(mc**2-ms**2)*ms2*(-3*trs+trg+trc)*gDs+ms2*(mc2
     &   *(trsgc+4*ms2*(-trs+trg+trc))-3*ms2**2*(-trs+trg+trc)+mc2**2*tr
     &   g))-2*(trsgc+mc2*(trg-3*trs))*gDs**3+2*cDs**2*((4*trg+6*trc)*gD
     &   s**2-(trsgc-3*mc2*trs+mc2*trg)*gDs+(mc2-ms2)*ms2*(-trs+trg+trc)
     &   )-4*ms2*(mc2*(trs+2*trg+3*trc)-ms2*trs)*gDs**2+ms2*(mc2*(3*trsg
     &   c+2*ms2*(-3*trs+trg+trc))-2*ms2**2*(-2*trs+trg+trc)+3*mc2**2*(t
     &   rg-2*trs))*gDs+ms2*(mc2**2*(trsgc-ms2*(-2*trs+trg+2*trc))-ms2**
     &   3*(-trs+trg+trc)+2*mc2*ms2**2*(-trs+trg+trc)))-ms2*cDg**6*(4*tr
     &   s*gDs**2-ms2*(-2*trs*gDs+3*trg*gDs+2*trg*cDs+trsgc)+ms2**2*(trs
     &   +trc))+ms2**2*trg*cDg**7)/(ms2*cDg**3*gDs**2*(gDs+cDg)**2)/4.0d
     &   +0+Bpmp

      Bpmp = lc*(4*(2*trg+3*trc)*gDs**4+2*cDs*((4*trg+6*trc)*gDs**3-(2*
     &   (trg-5*trs)*cDg+trsgc-3*mc2*trs+mc2*trg)*gDs**2+cDg*(-2*trs*cDg
     &   +trsgc+ms2*trg)*gDs+2*ms2*trg*cDg*(2*cDg+mc2))-2*(2*(trg-5*trs)
     &   *cDg+trsgc+mc2*(trg-3*trs))*gDs**3+2*cDg*(-2*trs*cDg+trsgc+ms2*
     &   trg)*gDs**2+ms2*(2*cDg+mc2)*(2*(-trs+trg+trc)*cDg+trsgc-mc2*trg
     &   )*gDs+ms2*cDg*(mc2*trsgc-(2*cDg+mc2)*(2*trs*cDg-ms2*trg)))/(ms2
     &   *cDg*(2*cDg+mc2))/2.0_dp+Bpmp

      Bpmp = LsB2*(2*mc2**2*trsgc*cDs*gDs**3+cDg**3*(mc2*ms2*(-4*mc2*tr
     &   s*gDs-ms2*trsgc+mc2*ms2*trg)+mc2*cDs**2*(8*trs*gDs-8*ms2*trg)+4
     &   *cDs**3*(-2*trs*gDs+trsgc+ms2*trg)-3*mc2*ms2*cDs*(-2*trs*gDs+tr
     &   sgc+ms2*trg)+8*trg*cDs**4)+2*cDg**2*gDs*(mc2**2*cDs*(3*ms2*trg-
     &   trs*gDs)-2*mc2*cDs**2*(-2*trs*gDs+trsgc+ms2*trg)+mc2**2*ms2*(-2
     &   *trs*gDs+trsgc+ms2*trg)+(2*trsgc-4*mc2*trg)*cDs**3)+mc2*cDg*gDs
     &   **2*(mc2*cDs*(-2*trs*gDs+trsgc+ms2*trg)+2*(mc2*trg-3*trsgc)*cDs
     &   **2+mc2*ms2*(trsgc-mc2*trg))+cDg**4*(6*mc2*ms2*trs*cDs-8*trs*cD
     &   s**3))/(cDg**3*(cDs**2-mc2*ms2))+Bpmp

      Bpmp = B0csf*(-cDg*gDs*(cDs*(2*(ms2-2*mc2)*trc*gDs**2+(mc2*(2*trs
     &   gc+ms2*(4*trs-3*trg+trc))+ms2*(ms2*trc-trsgc)+2*mc2**2*trs)*gDs
     &   -2*mc2*ms2*(trsgc+3*ms2*trg+2*mc2*trg))-mc2*ms2*(2*trc*gDs**2+(
     &   -trsgc-mc2*trs+ms2*(trg-trs))*gDs+(ms2+mc2)*(trsgc+ms2*trg))+cD
     &   s**2*(2*(mc2*(trs+trg-trc)-ms2*trs)*gDs+ms2*(trsgc+ms2*trg)+mc2
     &   *(trsgc-7*ms2*trg))+2*cDs**3*(-2*(trs-trg+trc)*gDs+trsgc+3*ms2*
     &   trg+2*mc2*trg)+8*trg*cDs**4)-cDg**2*(cDs*(-2*(mc2-2*ms2)*trc*gD
     &   s**2+(mc2*(trsgc+ms2*(5*trs-5*trg+4*trc))+2*ms2*(ms2*trc-trsgc)
     &   +mc2**2*trs)*gDs-2*mc2*ms2*(trsgc+(ms2+mc2)*trg))+mc2*ms2*(gDs*
     &   (2*trc*gDs+mc2*(2*trs-trg+trc)+ms2*(2*trs+trc))-trsgc*(gDs+ms2+
     &   mc2))+2*cDs**3*(-2*(2*trs-trg+trc)*gDs+trsgc+(ms2+mc2)*trg)+cDs
     &   **2*(-2*(2*ms2+mc2)*trs*gDs-2*(mc-ms)*(ms+mc)*(trc-trg)*gDs+mc2
     &   *(trsgc-4*ms2*trg)+ms2*trsgc)+4*trg*cDs**4)+mc2*gDs**2*(cDs*(2*
     &   trc*gDs**2+(-trsgc-mc2*trs+ms2*(-trs-trg+2*trc))*gDs-2*mc2*ms2*
     &   trg)+ms2*(2*trc*gDs**2+(-trsgc+mc2*(trc-trg)+ms2*trc)*gDs-mc2*(
     &   ms2+mc2)*trg)+cDs**2*((ms2+mc2)*trg-2*trs*gDs)+2*trg*cDs**3)+cD
     &   g**3*(-mc2*ms2*(2*trc*gDs-trsgc+ms2*(trs+trg)+mc2*trs)+ms2*cDs*
     &   (-2*trc*gDs+trsgc+mc2*(-2*trs+trg-trc)-ms2*trc)+4*trs*cDs**3+2*
     &   (ms2*(trs+trg-trc)+mc2*trs)*cDs**2))/(cDg*(cDs-mc*ms)*(cDs+mc*m
     &   s)*(gDs+cDg)**2)/2.0_dp+Bpmp

      Bpmp = epinv*(cDg**2*(4*trc*(mc*ms*(xs**2-1)-4*lp*xs*cDs)*gDs**2+
     &   (8*lp*trsgc*xs*cDs-2*mc*ms*(trsgc+2*ms2*(sck-1)*trc)*(xs**2-1))
     &   *gDs+ms2*trsgc*(mc*ms*(2*sck-1)*(xs**2-1)-4*lp*xs*cDs))+mc2*trs
     &   gc*(mc*ms*(xs**2-1)-4*lp*xs*cDs)*gDs**2-2*cDg*gDs*(mc2*trg*(mc*
     &   ms*(xs**2-1)-4*lp*xs*cDs)*gDs-4*lp*trsgc*xs*cDs**2-mc*ms*(xs**2
     &   -1)*(mc2*ms2*(sck-1)*trg-trsgc*cDs)))/(mc*ms*(xs**2-1)*cDg**2*g
     &   Ds)/4.0_dp+Bpmp

      Bpmp = ls*ms2*(cDg*(4*(trc-trs)*gDs**2+2*(trg*cDs-trsgc+ms2*(2*tr
     &   g+trc))*gDs+ms2*(trg*(cDs+ms2)-trsgc))+gDs*(2*gDs+ms2)*(2*trc*g
     &   Ds-trsgc-mc2*trg)+trg*cDg**2*(2*gDs+ms2))/(cDg*(2*gDs+ms2)**2)/
     &   2.0_dp+Bpmp

      Bpmp = tr2fu*cDs*(-cDg**2*(4*trc*gDs**2-2*trsgc*gDs+ms2*trsgc)-mc
     &   2*trsgc*gDs**2+2*cDg*gDs*(mc2*trg*gDs+trsgc*cDs))/(cDg**2*gDs)+
     &   Bpmp

      Bpmp = 3*ms2*tr3c001fs*(mc2*(trs+trg)*gDs**2+2*trc*cDg*gDs*(2*gDs
     &   +ms2)+(trs+trg)*cDg**2*(2*gDs+ms2))/(cDg*gDs**2)+Bpmp

      Bpmp = 3*ms2*tr3c002fs*(mc2*trc*gDs**2+trc*cDg**2*(2*gDs+ms2)+2*m
     &   c2*(trs+trg)*cDg*gDs)/(cDg*gDs**2)+Bpmp

      Bpmp = LsB1*(gDs*(2*trc*cDs*(2*cDg*gDs-ms2*(gDs+cDg))+trg*(-2*mc2
     &   *cDs*gDs+ms2*cDs*(2*cDs+mc2)-mc2*ms2**2))+trsgc*(ms2*((cDs+mc2)
     &   *gDs+2*cDg*cDs)-2*cDs*(cDs+cDg)*gDs))*(mc2*gDs**2-2*cDg*cDs*gDs
     &   +ms2*cDg**2)/(cDg*(cDs-mc*ms)*(cDs+mc*ms)*gDs**2)+Bpmp

      Bpmp = BfunX*(cDg*gDs**3*(2*(gDs+cDs)+ms2)*((trg+trc)*(2*(gDs+cDs
     &   )+3*ms2)-3*trs*(2*(gDs+cDs)+ms2))+gDs**3*(2*(gDs+cDs)+ms2)**2*(
     &   (trg+trc)*(cDs+ms2)-trs*(2*gDs+cDs+ms2))-2*ms2*cDg**3*gDs*(2*(t
     &   rs+trg)*(gDs+cDs)+ms2*(trs+trg-trc))+2*(-trs+trg+trc)*cDg**2*gD
     &   s**3*(2*(gDs+cDs)+ms2)+ms2*cDg**4*(2*(-trs-trg+trc)*gDs+ms2*(tr
     &   s+trg+trc)))/(cDg**3*gDs**2)/4.0_dp+Bpmp

      Bpmp = ms2*tr3c00fs*(gDs*(mc2*(3*trs*gDs+4*trg*gDs+trsgc+ms2*trg)
     &   +cDs*(-2*trc*gDs+trsgc+mc2*trg))-cDg*(cDs*(2*(trg+trc)*gDs+trsg
     &   c+ms2*trg)-gDs*(10*trc*gDs+mc2*trg+8*ms2*trc)+trsgc*(2*gDs+ms2)
     &   +trg*cDs**2)+cDg**2*(2*(3*trs+2*trg+trc)*gDs-2*trg*cDs-trsgc+ms
     &   2*(3*trs+2*(trg+trc)))-trg*cDg**3)/(cDg*gDs**2)+Bpmp

      Bpmp = 3*tr3s002ft*(2*cDg+mc2)*((trg+trc)*((2*cDg+mc2)*cDs*gDs-ms
     &   2*cDg**2)+mc2*trs*(2*cDg+mc2)*gDs)/cDg**3+Bpmp

      Bpmp = 3*tr3s001ft*(2*cDg+mc2)*(trs*(2*cDg+mc2)*cDs*gDs+mc2*ms2*(
     &   trg+trc)*gDs-ms2*trs*cDg**2)/cDg**3+Bpmp

      Bpmp=Bpmp/zpmp

      return
      end
