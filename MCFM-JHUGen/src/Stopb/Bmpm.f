      subroutine Bamp_mpm(q,mc,ms,Bmpm)
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
      complex(dp):: trc,trg,trs,trsgc,zmpm,Bmpm

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

      Bmpm = -mc*ms*ms2*(gDs+cDg)*(cDs*(-4*trc*gDs**3+(-2*trs*cDg+trsgc
     &   +3*mc2*trg)*gDs**2+cDg*(2*trs*cDg-trsgc-ms2*trg)*gDs+2*ms2*trg*
     &   cDg**2)+gDs*(cDg*(-2*mc2*trs*gDs+4*ms2*trc*gDs-ms2*trsgc+mc2*ms
     &   2*trg)+mc2*gDs*(2*trs*gDs+trsgc-3*ms2*trg))+2*trg*cDs**2*gDs*(2
     &   *gDs-3*cDg))*tr5Xs/(cDg*(cDs**2-mc2*ms2)*gDs**2)/2.0_dp

      Bmpm = mc*ms*(gDs+cDg)*(-2*mc2**2*trg*cDs*gDs**2+cDg**2*(mc2*(2*(
     &   trsgc-2*mc2*trs+ms2*trc)*gDs+ms2*(3*mc2*trg-trsgc))+mc2*cDs*(2*
     &   (trg+trc)*gDs-trsgc-3*ms2*trg)+4*trg*cDs**2*(gDs-mc2))-2*cDg**3
     &   *(mc2*(2*trs*gDs-2*trs*cDs-ms2*trg+ms2*trc)+cDs*(2*trg*cDs+trsg
     &   c+ms2*trg))+mc2*cDg*gDs*(cDs*(-2*(trg+trc)*gDs+trsgc+mc2*trg)+6
     &   *trg*cDs**2+mc2*(trsgc-ms2*trg))+4*trs*cDg**4*cDs)*tr5Xc/(cDg**
     &   3*(cDs**2-mc2*ms2))/2.0_dp+Bmpm

      Bmpm = mc*ms*ms2*(cDs*(-4*trc*gDs**3+(-2*trs*cDg+trsgc+3*mc2*trg)
     &   *gDs**2+cDg*(2*trs*cDg-trsgc-ms2*trg)*gDs+2*ms2*trg*cDg**2)+gDs
     &   *(cDg*(-2*mc2*trs*gDs+4*ms2*trc*gDs-ms2*trsgc+mc2*ms2*trg)+mc2*
     &   gDs*(2*trs*gDs+trsgc-3*ms2*trg))+2*trg*cDs**2*gDs*(2*gDs-3*cDg)
     &   )*tr4Xs/(cDg*(cDs**2-mc2*ms2)*gDs)/2.0_dp+Bmpm

      Bmpm = mc*ms*(2*mc2**2*trg*cDs*gDs**2+cDg**2*(mc2*(ms2*(trsgc-3*m
     &   c2*trg)-2*(trsgc-2*mc2*trs+ms2*trc)*gDs)+mc2*cDs*(-2*(trg+trc)*
     &   gDs+trsgc+3*ms2*trg)+4*trg*cDs**2*(mc2-gDs))+2*cDg**3*(mc2*(2*t
     &   rs*gDs-2*trs*cDs-ms2*trg+ms2*trc)+cDs*(2*trg*cDs+trsgc+ms2*trg)
     &   )-mc2*cDg*gDs*(cDs*(-2*(trg+trc)*gDs+trsgc+mc2*trg)+6*trg*cDs**
     &   2+mc2*(trsgc-ms2*trg))-4*trs*cDg**4*cDs)*tr4Xc/(cDg**2*(cDs**2-
     &   mc2*ms2))/2.0_dp+Bmpm

      Bmpm = mc*ms*(2*mc2**2*ms2*trg*gDs**2+cDg**2*(cDs*(ms2*(trsgc+mc2
     &   *trg)-2*(trsgc-2*mc2*trs+2*ms2*trg+ms2*trc)*gDs)+mc2*ms2*(-2*tr
     &   c*gDs+trsgc+3*ms2*trg)-2*trg*cDs**2*gDs)+cDg*gDs*(2*(trg*cDs**2
     &   +mc2*ms2*trc)*gDs-mc2*(ms2*trg*(5*cDs+mc2)+trsgc*(cDs+ms2)))+2*
     &   cDg**3*(2*trs*cDs*gDs+ms2*(trg+trc)*cDs+ms2*(trsgc-2*mc2*trs+ms
     &   2*trg))-4*ms2*trs*cDg**4)*tr3Xs/(cDg**2*(cDs-mc*ms)*(cDs+mc*ms)
     &   )/2.0_dp+Bmpm

      Bmpm = mc*ms*(cDg*gDs*(cDs*(-2*mc2*trs*gDs+4*ms2*trc*gDs-ms2*trsg
     &   c-5*mc2*ms2*trg)-mc2*ms2*(2*trs*gDs+trsgc+ms2*trg))+mc2*gDs**2*
     &   (2*trs*cDs*gDs-4*ms2*trc*gDs+trsgc*cDs+ms2*trg*cDs+ms2*trsgc+3*
     &   mc2*ms2*trg)+2*mc2*ms2*cDg**2*(trs*gDs+ms2*trg))*tr3Xc/(cDg*(cD
     &   s**2-mc**2*ms**2)*gDs)/2.0_dp+Bmpm

      Bmpm = ms*(-cDg**3*gDs*(cDs**2*(4*trc*gDs**2+(-2*trsgc+6*mc2*sck*
     &   trs-4*mc2*trs+4*ms2*trc)*gDs-ms2*(trsgc+mc2*(3*sck-2)*trg))+mc2
     &   *ms2*(-4*trc*gDs**2+(4*trsgc+mc2*(2-6*sck)*trs-2*ms2*trc)*gDs+m
     &   s2*(trsgc+mc2*(3*sck-2)*trg))+2*mc2*cDs*gDs*(2*trs*gDs+ms2*(-tr
     &   s+trg+trc)))+ms2*cDg**4*(4*mc2*trs*gDs**2+trc*(ms2-2*cDs)*(mc2*
     &   ms2-cDs**2))+mc2*cDg**2*gDs**2*(cDs**2*(-2*trs*gDs+6*trg*gDs+3*
     &   sck*trsgc-2*trsgc)+mc2*ms2*(2*(trc-2*trg)*gDs+(2-3*sck)*trsgc)+
     &   2*(trsgc-mc2*trs+ms2*trc)*cDs*gDs)+mc2*(cDs-mc*ms)*(cDs+mc*ms)*
     &   gDs**3*(trs*(cDs+ms2)*(2*gDs+2*cDs+ms2-mc2)+ms2*(trg+trc)*gDs)+
     &   mc2*cDg*(cDs-mc*ms)*(cDs+mc*ms)*gDs**3*(-2*trg*gDs+ms2*(trs+2*(
     &   trg+trc))-mc2*trs)+2*ms2*trc*cDg**5*(cDs-mc*ms)*(cDs+mc*ms))*lV
     &   s/(mc*cDg**3*(cDs-mc*ms)*(cDs+mc*ms)*gDs**2)/4.0_dp+Bmpm

      Bmpm = mc*(cDg**2*gDs**2*(cDs**2*(2*trg*gDs**2-2*(-trsgc+3*mc2*tr
     &   s+ms2*trc)*gDs+ms2*(trsgc+2*mc2*trg))-mc2*ms2*(2*trg*gDs**2-4*(
     &   mc2*trs+ms2*trc)*gDs+ms2*(trsgc+2*mc2*trg))+2*trg*cDs**3*gDs+2*
     &   mc2*ms2*(trc-trs)*cDs*gDs)+ms2*cDg**4*(4*trs*cDs*gDs**2-trc*(cD
     &   s-mc*ms)*(cDs+mc*ms)*(2*cDs-ms2))-2*cDg**3*gDs**2*(mc2*ms2*(ms2
     &   *(-trs+trg+trc)-2*trs*gDs)+4*trs*cDs**2*gDs+ms2*(trsgc-mc2*trs+
     &   ms2*trc)*cDs)+mc2*(mc2*ms2-cDs**2)*gDs**3*(trs*(cDs+ms2)*(2*gDs
     &   +2*cDs+ms2-mc2)+ms2*(trg+trc)*gDs)-cDg*(cDs-mc*ms)*(cDs+mc*ms)*
     &   gDs**3*(2*trg*gDs*(gDs+cDs)+mc2*(ms2*(trs+trg+2*trc)-trsgc))+2*
     &   ms2*trc*cDg**5*(mc2*ms2-cDs**2))*lVc/(ms*cDg**3*(cDs-mc*ms)*(cD
     &   s+mc*ms)*gDs**2)/4.0_dp+Bmpm

      Bmpm = Bmpm-trg*xs*cDs*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)*lRc
     &   s/((xs**2-1)*cDg**2*gDs)

      Bmpm = (-ms2*cDg**3*(mc2*((60*sck*trs+72*trs+40*trg-52*trc)*gDs**
     &   3+2*(9*trsgc+ms2*(15*sck*trs+55*trs-15*sck*trg+71*trg-29*trc))*
     &   gDs**2+ms2*(27*trsgc+ms2*(34*trs-15*sck*trg+76*trg-16*trc))*gDs
     &   +9*ms2**2*(trsgc+ms2*trg))+9*cDs*(2*gDs+ms2)*(-4*trc*gDs**2+2*(
     &   trsgc+mc2*trs+2*mc2*trg-2*ms2*trc)*gDs+ms2*(trsgc+mc2*trg)))+mc
     &   2*cDg**2*gDs*((2*gDs+ms2)*(72*trs*gDs**3-8*ms2*(9*trc-2*trg)*gD
     &   s**2+ms2*(-36*mc2*trs+29*ms2*trg+23*mc2*trg-36*ms2*trc)*gDs+9*m
     &   c2*ms2**2*trg)+cDs*(2*gDs+ms2)*(68*trs*gDs**2+18*ms2*(trs+2*trg
     &   +trc)*gDs+9*ms2*(trsgc+ms2*trg))+3*ms2*trsgc*(2*(5*sck+6)*gDs**
     &   2+ms2*(5*sck+13)*gDs+3*ms2**2))-mc2*cDg*gDs**2*(2*gDs+ms2)*(18*
     &   trg*gDs**3+cDs*(36*trg*gDs**2+2*(9*trsgc-34*mc2*trs+9*ms2*(2*tr
     &   g+trc))*gDs-9*ms2*(trsgc+mc2*trg))-2*(-9*trsgc+27*mc2*trs-10*ms
     &   2*trg+17*ms2*trc)*gDs**2+18*trg*cDs**2*gDs+ms2*(9*trsgc+mc2*(10
     &   *trs+19*trg+28*trc)+3*ms2*(trs+2*(trg+trc)))*gDs-9*mc2*ms2*(trs
     &   gc+ms2*trg))+mc2*gDs**3*(2*gDs+ms2)*(18*trg*gDs**3+36*trg*cDs*g
     &   Ds**2+(18*trg*cDs**2+mc2*(6*trs*cDs-9*trsgc+ms2*(6*trs+2*trg+29
     &   *trc))-6*ms2*trs*cDs-3*ms2**2*(2*trs+trg+trc))*gDs+6*(mc-ms)*(m
     &   s+mc)*trs*cDs**2+(-3*mc2*(3*trsgc-4*ms2*trs+3*ms2*trg)-9*ms2**2
     &   *trs+5*mc2**2*trs)*cDs-ms2*(mc2*(9*trsgc-6*ms2*trs)+mc2**2*(11*
     &   trs+9*trg)+3*ms2**2*trs))+ms2*cDg**4*(2*gDs+ms2)*(40*trc*gDs**2
     &   +2*(-9*trsgc+17*mc2*trs+8*mc2*trg+11*ms2*trc)*gDs+6*(mc-ms)*(ms
     &   +mc)*trc*cDs+ms2*(-9*trsgc+mc2*(26*trs+17*trg-3*trc)-5*ms2*trc)
     &   )+6*(mc-ms)*(ms+mc)*ms2*trc*cDg**5*(2*gDs+ms2))/(mc*ms*cDg**3*g
     &   Ds**2*(2*gDs+ms2))/1.2e+1_dp+Bmpm

      Bmpm = B0cgsf*(-cDg**3*(-4*(2*mc2*trs+ms2*trc)*gDs**3+ms2*cDs*(-4
     &   *trc*gDs**2+2*(trsgc+mc2*trs+2*mc2*trg-2*ms2*trc)*gDs+ms2*(trsg
     &   c+mc2*trg))+2*ms2*(trsgc+mc2*(trs+3*trg-2*trc)-2*ms2*trc)*gDs**
     &   2+ms2*(mc2*(trsgc+4*ms2*trs+6*ms2*trg-2*ms2*trc)+ms2*trsgc)*gDs
     &   +mc2*ms2**2*(trsgc+ms2*trg))-mc2*cDg**2*gDs*(2*(trg-4*trs)*gDs*
     &   *3+cDs*(2*(trg-4*trs)*gDs**2-2*ms2*(trs+trc)*gDs-ms2*(trsgc+ms2
     &   *trg))+2*(trsgc-3*mc2*trs+ms2*(-trs+2*trg+trc))*gDs**2+ms2*(-2*
     &   trsgc+4*mc2*trs-ms2*trg+mc2*trg+4*ms2*trc)*gDs-ms2**2*(trsgc+mc
     &   2*trg))-mc2*gDs**3*(-2*trg*gDs**3-4*trg*cDs*gDs**2+(-2*trg*cDs*
     &   *2-mc2*(2*trs*cDs-trsgc+ms2*(2*trs+trg+4*trc))+2*ms2*trs*cDs+ms
     &   2**2*(2*trs+trg+trc))*gDs+2*(ms2-mc2)*trs*cDs**2+(mc2*(trsgc+ms
     &   2*(trg-4*trs))+3*ms2**2*trs)*cDs+ms2*(mc2*(trsgc-2*ms2*trs)+mc2
     &   **2*(2*trs+trg)+ms2**2*trs))-mc2*cDg*gDs**2*(cDs*(2*trg*gDs**2+
     &   2*(trsgc-4*mc2*trs+ms2*(2*trg+trc))*gDs-ms2*(trsgc+mc2*trg))+(2
     &   *trsgc-6*mc2*trs-4*ms2*trc)*gDs**2+2*trg*cDs**2*gDs+((ms2+mc2)*
     &   trsgc+ms2**2*trs+2*ms2*(ms2+mc2)*(trg+trc))*gDs-mc2*ms2*(trsgc+
     &   ms2*trg))-ms2*cDg**4*(-4*trc*gDs**2-2*(-trsgc+2*mc2*trs+mc2*trg
     &   +ms2*trc)*gDs+2*(ms**2-mc**2)*trc*cDs+ms2*(trsgc+mc2*(-3*trs-2*
     &   trg+trc)))-2*ms2*(ms2-mc2)*trc*cDg**5)/(mc*ms*cDg**3*gDs**2)/4.
     &   0_dp+Bmpm

      Bmpm = mc*tr3s00ft*(cDg*(2*trg*gDs**3+cDs*(4*trg*gDs**2+2*(trsgc-
     &   2*mc2*trs+2*ms2*trg+ms2*trc)*gDs-ms2*(trsgc+mc2*trg))-2*(-trsgc
     &   +3*mc2*trs+ms2*(trg+4*trc))*gDs**2+2*trg*cDs**2*gDs+ms2*(trsgc+
     &   2*mc2*(3*trs+trg+trc))*gDs-mc2*ms2*(trsgc+ms2*trg))+gDs*(-2*trg
     &   *gDs**3-4*trg*cDs*gDs**2+(mc2*(trsgc-2*ms2*trg-5*ms2*trc)-2*trg
     &   *cDs**2)*gDs+mc2*((trsgc+ms2*trg)*cDs+ms2*(trsgc+3*mc2*trs+mc2*
     &   trg)))-cDg**2*(8*trs*gDs**2+2*cDs*(2*trs*gDs+ms2*(trg+trc))+ms2
     &   *(6*trg+4*trc)*gDs+ms2*(trsgc-4*mc2*trs+5*ms2*trg+3*ms2*trc))+4
     &   *ms2*trs*cDg**3)/(ms*cDg**3)+Bmpm

      Bmpm = Bmpm-LsB1*mc*ms*(ms2*gDs*(4*mc2*trc*gDs**3-mc2*(3*trsgc+m
     &   c2*trg)*gDs**2+4*ms2*cDg*(mc2*trg-trc*cDg)*gDs+ms2*(trsgc-mc2*t
     &   rg)*cDg**2)-2*cDs**2*gDs*(4*trc*gDs**3-2*(-trs*cDg+trsgc+mc2*tr
     &   g)*gDs**2+cDg*(-2*trs*cDg+trsgc+3*ms2*trg)*gDs-4*ms2*trg*cDg**2
     &   )+cDs*(2*cDg*gDs**2*(-mc2*trs*gDs+4*ms2*trc*gDs-ms2*trsgc)+ms2*
     &   cDg**2*gDs*(2*trs*gDs+trsgc+ms2*trg)+mc2*gDs**3*(2*trs*gDs+trsg
     &   c-7*ms2*trg)-2*ms2*cDg**3*(trs*gDs+ms2*trg))+8*trg*cDs**3*gDs**
     &   2*(gDs-cDg))/(cDg*(cDs**2-mc2*ms2)*gDs**2)

      Bmpm = ls*ms*(cDg*(8*trc*gDs**3+4*(-trsgc-2*mc2*trs+3*ms2*trc)*gD
     &   s**2+2*ms2*(-2*trsgc-mc2*trs+2*mc2*trg+2*ms2*trc)*gDs+ms2**2*(m
     &   c2*trg-trsgc))+mc2*gDs*(trsgc*(4*gDs+ms2)-(2*gDs+ms2)*(-2*trs*g
     &   Ds+4*trg*gDs+3*ms2*trg)))/(mc*cDg*(2*gDs+ms2)**2)/2.0_dp+Bmpm

      Bmpm = Bmpm-LsB2*mc*ms*(-2*mc2**2*trg*cDs*gDs**3+cDg**3*(mc2*(-4
     &   *trs*gDs**2-4*ms2*trc*gDs+ms2*(3*trsgc+ms2*trg))+cDs*(ms2*(trsg
     &   c+mc2*trg)-2*(2*trsgc-4*mc2*trs+3*ms2*trg+ms2*trc)*gDs)-2*cDs**
     &   2*(2*trg*gDs+trsgc-ms2*trg))+2*cDg**2*gDs*(mc2*((trsgc-2*mc2*tr
     &   s+ms2*trg+2*ms2*trc)*gDs+ms2*(mc2*trg-trsgc))+mc2*cDs*((trg+trc
     &   )*gDs-trsgc-2*ms2*trg)+2*trg*cDs**2*(gDs-mc2)-2*trg*cDs**3)+2*c
     &   Dg**4*(cDs*(4*trs*gDs+ms2*(trg+trc))+ms2*(trsgc-2*mc2*trs+ms2*t
     &   rg))+mc2*cDg*gDs**2*(cDs*(-2*(trg+trc)*gDs+trsgc+mc2*trg)+8*trg
     &   *cDs**2+mc2*(trsgc-ms2*trg))-4*ms2*trs*cDg**5)/(cDg**3*(cDs**2-
     &   mc2*ms2))

      Bmpm = ms*tr3c00fs*(cDg*(mc2*(10*(trs+trg)*gDs**2+(trsgc+8*ms2*tr
     &   s+10*ms2*trg)*gDs+ms2*(trsgc+ms2*trg))+cDs*(-4*trc*gDs**2+2*(tr
     &   sgc+mc2*trs+2*mc2*trg-2*ms2*trc)*gDs+ms2*(trsgc+mc2*trg)))+mc2*
     &   gDs*(12*trc*gDs**2+(-2*(trs+2*trg)*cDs-3*trsgc-2*mc2*trg+7*ms2*
     &   trc)*gDs-ms2*trg*(cDs+mc2)+trsgc*(-cDs-ms2))+cDg**2*(4*trc*gDs*
     &   *2+2*(trsgc-mc2*trs+3*ms2*trc)*gDs+ms2*(trsgc-mc2*(2*trs+trg)+3
     &   *ms2*trc)))/(mc*cDg*gDs**2)+Bmpm

      Bmpm = 3*mc*ms*tr3s001ft*(mc2*trs*(2*cDg+mc2)*gDs-(trg+trc)*((2*c
     &   Dg+mc2)*gDs**2+ms2*cDg**2))/cDg**3+Bmpm

      Bmpm = lc*mc*(-cDg*(2*trg*gDs**2+2*(trg*cDs+trsgc-3*mc2*trs+ms2*t
     &   rc)*gDs+ms2*(3*mc2*trg-trsgc))+gDs*(2*trg*gDs**2+2*trg*cDs*gDs-
     &   mc2*trsgc+mc2*ms2*trg)+2*cDg**2*(4*trs*gDs+ms2*(trc-2*trg)))/(m
     &   s*cDg*(2*cDg+mc2))/2.0_dp+Bmpm

      Bmpm = epinv*(mc2*trg*(mc*ms*(xs**2-1)-4*lp*xs*cDs)*gDs**2+ms*cDg
     &   **2*(mc*(xs**2-1)*(ms2*(2*sck-1)*trg-4*(sck-1)*trs*gDs)-4*lp*ms
     &   *trg*xs*cDs)+2*cDg*(4*lp*trg*xs*cDs**2-mc*ms*trg*(xs**2-1)*cDs+
     &   mc*ms*(sck-1)*trsgc*(xs**2-1))*gDs)/((xs**2-1)*cDg**2*gDs)/4.0d
     &   +0+Bmpm

      Bmpm = 3*mc*ms*tr3c002fs*(mc2*(trs+trg)*gDs**2+2*trc*cDg*gDs*(2*g
     &   Ds+ms2)+(trs+trg)*cDg**2*(2*gDs+ms2))/(cDg*gDs**2)+Bmpm

      Bmpm = 3*ms*tr3c001fs*(2*gDs+ms2)*(mc2*trc*gDs**2+trc*cDg**2*(2*g
     &   Ds+ms2)+2*mc2*(trs+trg)*cDg*gDs)/(mc*cDg*gDs**2)+Bmpm

      Bmpm = Bmpm-mc*ms*tr2fu*trg*cDs*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cD
     &   g**2)/(cDg**2*gDs)

      Bmpm = BfunX*(mc2*gDs**3*(2*(gDs+cDs)+ms2)*(trs*(cDs+ms2)*(2*(gDs
     &   +cDs)+ms2)+ms2*(trg+trc)*gDs)+ms2*cDg**4*(mc2*(2*(trs+trg)*gDs+
     &   ms2*(trs+trg-trc))+2*trc*cDs*(2*gDs+mc2)+4*trc*cDs**2)+mc2*cDg*
     &   (2*trs*cDs+ms2*(3*trs+2*(trg+trc)))*gDs**3*(2*(gDs+cDs)+ms2)+2*
     &   mc2*ms2*cDg**2*gDs**2*((trs+trc)*gDs-trg*(cDs+ms2))+2*mc2*ms2*c
     &   Dg**3*gDs*((trs+trc)*gDs+ms2*(trs+trg+trc))+2*ms2*trc*cDg**5*(2
     &   *gDs+4*cDs+mc2)+4*ms2*trc*cDg**6)/(mc*ms*cDg**3*gDs**2)/4.0_dp+
     &   Bmpm

      Bmpm = Bmpm-B0csf*mc*ms*(-cDg*(mc2*(2*trs*gDs+ms2*(trc-trg))+ms2
     &   *(trsgc+ms2*trc))+cDs*((-2*trs*cDg+trsgc+mc2*(-2*trs+trg+trc)+m
     &   s2*trc)*gDs+2*trs*cDg**2-(trsgc-mc2*trs+ms2*(-trs+trg+2*trc))*c
     &   Dg+2*mc2*ms2*trg)+mc2*((trsgc-mc2*trs+ms2*(trg-trs))*gDs+ms2*(m
     &   s2+mc2)*trg)+cDs**2*(2*trc*gDs-trg*(2*cDg+ms2+mc2)+2*trs*cDg)-2
     &   *trg*cDs**3+2*ms2*trs*cDg**2)/(cDg*(cDs-mc*ms)*(cDs+mc*ms))/2.0
     &   _dp

      Bmpm = 3*mc*tr3s002ft*(2*cDg+mc2)*(trs*((2*cDg+mc2)*cDs*gDs-ms2*c
     &   Dg**2)-2*ms2*(trg+trc)*cDg*gDs)/(ms*cDg**3)+Bmpm

      Bmpm=Bmpm/zmpm

      return
      end
