      subroutine Bamp_mpp(q,mc,ms,Bmpp)
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
      complex(dp):: trc,trg,trs,trsgc,zmpp,Bmpp

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
      zmpp=za(4,3)**2*zb(2,4)*zb(2,3)

      Bmpp = mc*ms2*(gDs+cDg)*(cDs*(-mc2*(trsgc+3*ms2*trg)*gDs**3+2*cDg
     &   *gDs**2*(mc2*trs*gDs+2*ms2*trc*gDs-ms2*trsgc)+2*ms2*cDg**3*(trs
     &   *gDs-ms2*trg)+ms2*(ms2*trg-trsgc)*cDg**2*gDs)+ms2*gDs*(2*mc2*tr
     &   c*gDs**3-mc2*(trsgc+mc2*trg)*gDs**2+2*ms2*cDg*(mc2*trg-trc*cDg)
     &   *gDs+ms2*(trsgc-mc2*trg)*cDg**2)+2*cDs**2*gDs*(-2*trc*gDs**3+(t
     &   rsgc+mc2*trg)*gDs**2+cDg*(-2*trs*cDg+trsgc-2*ms2*trg)*gDs+3*ms2
     &   *trg*cDg**2)+4*trg*cDs**3*gDs**2*(gDs-cDg))*tr5Xs/(cDg*(cDs**2-
     &   mc2*ms2)*gDs**3)/2.0_dp

      Bmpp = mc*(gDs+cDg)*(-2*mc2**2*trsgc*cDs*gDs**3+cDg**3*(cDs*(4*mc
     &   2*trs*gDs**2+2*ms2*(-trsgc+3*mc2*trs-2*mc2*trg+2*mc2*trc)*gDs-m
     &   c2*ms2*(trsgc+3*ms2*trg))-4*cDs**2*(2*trs*gDs**2-(trsgc-2*mc2*t
     &   rs+ms2*trg)*gDs+mc2*ms2*trg)+8*trg*cDs**3*gDs+mc2*ms2**2*(3*mc2
     &   *trg-trsgc))+2*cDg**2*gDs*(-mc2*cDs*(-2*trs*gDs**2+(trsgc-2*mc2
     &   *trs+ms2*trg)*gDs+3*mc2*ms2*trg)+cDs**2*(mc2*(2*trsgc+3*ms2*trg
     &   )-2*(mc2*(3*trs+trg+trc)-trsgc)*gDs)+mc2**2*ms2*((trg+trc)*gDs-
     &   trsgc)+4*mc2*trg*cDs**3)-2*cDg**4*(cDs*(2*cDs*(2*trs*gDs+ms2*tr
     &   g)+ms2*(-2*trs*gDs+trsgc+ms2*trg))+mc2*ms2*(ms2*(trc-trg)-2*trs
     &   *cDs))-mc2*cDg*gDs**2*(cDs*(mc2*(-6*trs*gDs+trsgc+3*ms2*trg)+2*
     &   trsgc*gDs)+2*(mc2*trg-2*trsgc)*cDs**2+mc2*ms2*(trsgc-mc2*trg))+
     &   4*ms2*trs*cDg**5*cDs)*tr5Xc/(cDg**3*(cDs**2-mc**2*ms**2)*gDs)/2
     &   .0_dp+Bmpp

      Bmpp = Bmpp-mc*ms2*(cDs*(-mc2*(trsgc+3*ms2*trg)*gDs**3+2*cDg*gDs
     &   **2*(mc2*trs*gDs+2*ms2*trc*gDs-ms2*trsgc)+2*ms2*cDg**3*(trs*gDs
     &   -ms2*trg)+ms2*(ms2*trg-trsgc)*cDg**2*gDs)+ms2*gDs*(2*mc2*trc*gD
     &   s**3-mc2*(trsgc+mc2*trg)*gDs**2+2*ms2*cDg*(mc2*trg-trc*cDg)*gDs
     &   +ms2*(trsgc-mc2*trg)*cDg**2)+2*cDs**2*gDs*(-2*trc*gDs**3+(trsgc
     &   +mc2*trg)*gDs**2+cDg*(-2*trs*cDg+trsgc-2*ms2*trg)*gDs+3*ms2*trg
     &   *cDg**2)+4*trg*cDs**3*gDs**2*(gDs-cDg))*tr4Xs/(cDg*(cDs**2-mc2*
     &   ms2)*gDs**2)/2.0_dp

      Bmpp = mc*(2*mc2**2*trsgc*cDs*gDs**3+cDg**3*(cDs*(-4*mc2*trs*gDs*
     &   *2+2*ms2*(trsgc-3*mc2*trs+2*mc2*trg-2*mc2*trc)*gDs+mc2*ms2*(trs
     &   gc+3*ms2*trg))+4*cDs**2*(2*trs*gDs**2-(trsgc-2*mc2*trs+ms2*trg)
     &   *gDs+mc2*ms2*trg)-8*trg*cDs**3*gDs+mc2*ms2**2*(trsgc-3*mc2*trg)
     &   )+2*cDg**2*gDs*(mc2*(cDs**2*(2*(3*trs+trg+trc)*gDs-3*ms2*trg)+m
     &   s2*trg*cDs*(gDs+3*mc2)-2*trs*cDs*gDs*(gDs+mc2)-mc2*ms2*(trg+trc
     &   )*gDs-4*trg*cDs**3)+trsgc*(-2*cDs**2*(gDs+mc2)+mc2*cDs*gDs+mc2*
     &   *2*ms2))+2*cDg**4*(cDs*(2*cDs*(2*trs*gDs+ms2*trg)+ms2*(-2*trs*g
     &   Ds+trsgc+ms2*trg))+mc2*ms2*(ms2*(trc-trg)-2*trs*cDs))+mc2*cDg*g
     &   Ds**2*(cDs*(mc2*(-6*trs*gDs+trsgc+3*ms2*trg)+2*trsgc*gDs)+2*(mc
     &   2*trg-2*trsgc)*cDs**2+mc2*ms2*(trsgc-mc2*trg))-4*ms2*trs*cDg**5
     &   *cDs)*tr4Xc/(cDg**2*(cDs-mc*ms)*(cDs+mc*ms)*gDs)/2.0_dp+Bmpp

      Bmpp = Bmpp-mc*ms2*(-2*mc2**2*trsgc*gDs**3+cDg**3*(4*trs*(mc2-2*
     &   cDs)*gDs**2+2*(trg*cDs**2+mc2*(ms2*(3*trs+trg+2*trc)-4*trs*cDs)
     &   +2*(trsgc+ms2*trg)*cDs-ms2*trsgc)*gDs-ms2*(mc2*trg*(cDs+3*ms2)+
     &   trsgc*(cDs+mc2)))+2*cDg**2*gDs*(mc2*(-cDs*((6*trs+trg+trc)*gDs-
     &   3*ms2*trg)+2*trs*gDs*(gDs+mc2)-ms2*trg*(gDs-mc2))+trsgc*(cDs*(2
     &   *gDs+mc2)-mc2*gDs))-2*cDg**4*(ms2*(-2*trs*(gDs+mc2)+trsgc+ms2*t
     &   rg)+4*trs*cDs*gDs+ms2*(trg+trc)*cDs)+mc2*cDg*gDs**2*(-mc2*(-6*t
     &   rs*gDs+trg*cDs+3*ms2*trg)-trsgc*(2*gDs-3*cDs+mc2))+4*ms2*trs*cD
     &   g**5)*tr3Xs/(cDg**2*(cDs-mc*ms)*(cDs+mc*ms)*gDs)/2.0_dp

      Bmpp = mc*(4*ms2*cDs*gDs*(2*mc2*trc*gDs**3-mc2*(trsgc+mc2*trg)*gD
     &   s**2+2*cDg*(2*mc2*trs*cDg+ms2*trc*cDg-mc2*trsgc+mc2*ms2*trg)*gD
     &   s-ms2*(trsgc+5*mc2*trg)*cDg**2)+mc2*ms2*(-4*cDg*gDs**2*((-4*trs
     &   gc+4*mc2*trs+mc2*trg+4*ms2*trc)*gDs-2*ms2*(trsgc+mc2*trg))+8*ms
     &   2*cDg**3*(ms2*trg-trs*gDs)-mc2*gDs**3*(2*(trg-trs)*gDs-8*trsgc+
     &   ms2*(trs+trg))-4*ms2*(3*ms2*(trs+trg)-trsgc)*cDg**2*gDs)+cDs**2
     &   *gDs*(4*cDg*gDs*((-4*trsgc+2*mc2*trs+mc2*trg)*gDs+2*mc2*ms2*trg
     &   )+mc2*gDs**2*(2*(trg-trs)*gDs-4*trsgc+ms2*trs-3*ms2*trg)+4*ms2*
     &   *2*(3*trs+2*trg)*cDg**2))*tr3Xc/(cDg*(cDs**2-mc**2*ms**2)*gDs**
     &   2)/8.0_dp+Bmpp

      Bmpp = (mc2*cDg**3*gDs*(mc2*ms2*(6*(3*trs+trg)*gDs**3+(ms2*(-3*tr
     &   s+trg+16*trc)-4*trsgc)*gDs**2-ms2*(8*trsgc+5*ms2*(trs+trg))*gDs
     &   +5*ms2**3*(trs+trg))+cDs**2*(-2*(trs+3*trg)*gDs**3+(4*trsgc-13*
     &   ms2*trs+15*ms2*trg)*gDs**2+ms2*(8*trsgc+5*ms2*(trs+trg))*gDs-5*
     &   ms2**3*(trs+trg))+4*ms2*cDs*gDs*((4*trsgc-10*mc2*trs-3*mc2*trg+
     &   4*ms2*trc)*gDs+2*mc2*ms2*trg)+4*cDs**3*gDs*(3*(2*trs+trg)*gDs-2
     &   *ms2*trg))-4*cDg**4*gDs*(cDs**2*((-4*trsgc+2*mc2*trs+mc2*trg+6*
     &   ms2*trc)*gDs**2+2*ms2*(-2*trsgc+2*mc2*trs+mc2*trg+ms2*trc)*gDs-
     &   ms2**2*(trsgc+mc2*trg+6*mc2*trc))+mc2*ms2*(-(-4*trsgc+2*mc2*trs
     &   +mc2*trg+6*ms2*trc)*gDs**2-2*ms2*(-3*trsgc+3*mc2*trs+mc2*trg)*g
     &   Ds+ms2**2*(trsgc+mc2*trg))+mc2*ms2*cDs*(8*trs*gDs**2+2*ms2*(-tr
     &   s+trg+trc)*gDs+ms2**2*(-3*trs-2*trg+3*trc))+6*ms2*trc*cDs**4+ms
     &   2**2*(3*trs+2*trg-3*trc)*cDs**3)-mc2*cDg**2*gDs**3*(4*mc2*ms2*(
     &   2*trc*gDs**2+(2*trsgc-2*mc2*(3*trs+trg)+ms2*(trs+trg+trc))*gDs+
     &   ms2*(-trsgc+3*ms2*(trs+trg+trc)-mc2*trg))-4*cDs**2*(2*trc*gDs**
     &   2+(ms2*(trs+trg-trc)-2*mc2*(2*trs+trg))*gDs+ms2*(-trsgc+3*ms2*(
     &   trs+trg+trc)-mc2*trg))+mc2*ms2*cDs*((8*trc-2*(trs+3*trg))*gDs-1
     &   2*trsgc-ms2*trs+3*ms2*trg)+cDs**3*(-6*trs*gDs+14*trg*gDs+12*trs
     &   gc+ms2*trs-3*ms2*trg))+4*mc2*ms2*(cDs-mc*ms)*(cDs+mc*ms)*gDs**4
     &   *((trg+trc)*(cDs+ms2)*(2*gDs+2*cDs+ms2-mc2)-trs*(4*gDs**2+cDs*(
     &   6*gDs+3*ms2-mc2)+4*ms2*gDs-2*mc2*gDs+2*cDs**2+ms2**2-mc2*ms2+mc
     &   2**2))-4*ms2*cDg**5*((cDs-mc*ms)*(cDs+mc*ms)*(trc*(4*cDs*(3*gDs
     &   +cDs)-ms2*(gDs+2*cDs)+ms2**2)+2*ms2*trg*gDs)-ms2*trs*gDs*(4*mc2
     &   *gDs-3*cDs**2+3*mc2*ms2))-4*mc2*cDg*(cDs-mc*ms)*(cDs+mc*ms)*gDs
     &   **4*(ms2*(2*trs*(2*gDs+ms2)+(trg+trc)*(-2*gDs-2*ms2+mc2)+6*trs*
     &   cDs)-2*mc2*trsgc)-8*ms2*trc*cDg**6*(mc2*ms2-cDs**2)*(-3*gDs-4*c
     &   Ds+ms2)-16*ms2*trc*cDg**7*(cDs-mc*ms)*(cDs+mc*ms))*lVs/(mc*cDg*
     &   *3*(cDs-mc*ms)*(cDs+mc*ms)*gDs**3)/1.6e+1_dp+Bmpp

      Bmpp = Bmpp-mc*(cDg**3*gDs*(cDs**2*(8*trg*gDs**3+4*(-6*trsgc+6*m
     &   c2*trs+mc2*trg)*gDs**2+4*ms2*(mc2*trg-3*trsgc)*gDs-5*mc2*ms2**2
     &   *(trs+trg))+mc2*ms2*(-8*trg*gDs**3-4*(-2*trsgc+2*mc2*trs+mc2*tr
     &   g+4*ms2*trc)*gDs**2-4*ms2*(mc2*trg-3*trsgc)*gDs+5*mc2*ms2**2*(t
     &   rs+trg))+8*(trg+2*trc)*cDs**3*gDs**2-8*mc2*cDs*gDs**2*(2*trs*gD
     &   s+ms2*(-2*trs+3*trg+4*trc)))+4*cDg**4*gDs*(cDs**2*(2*(3*trs+trc
     &   )*gDs**2+ms2*trc*gDs+ms2*(3*ms2*trs+2*ms2*trg+6*mc2*trc))+mc2*m
     &   s2*(-2*(trc-trs)*gDs**2+ms2*(-2*trs+2*trg+trc)*gDs-ms2**2*(3*tr
     &   s+2*trg))+ms2*cDs*(2*(trsgc-mc2*trs+ms2*trc)*gDs-3*mc2*ms2*trc)
     &   -6*trc*cDs**4+3*ms2*trc*cDs**3)+4*mc2*(cDs-mc*ms)*(cDs+mc*ms)*g
     &   Ds**4*((trg+trc)*(cDs+ms2)*(2*gDs+2*cDs+ms2-mc2)-trs*(4*gDs**2+
     &   cDs*(6*gDs+3*ms2-mc2)+4*ms2*gDs-2*mc2*gDs+2*cDs**2+ms2**2-mc2*m
     &   s2+mc2**2))+4*cDg**5*(-4*ms2*trs*cDs*gDs**2-trc*(cDs-mc*ms)*(cD
     &   s+mc*ms)*(4*cDs*(3*gDs+cDs)-ms2*(gDs+2*cDs)+ms2**2))+cDg**2*gDs
     &   **3*(4*trsgc*(cDs*(2*cDs*gDs+2*mc2*gDs+4*cDs**2-mc2*cDs)+mc2*ms
     &   2*(-2*gDs-4*cDs+mc2))+mc2*(-mc2*ms2*(2*(5*trs+11*trg+4*trc)*gDs
     &   +ms2*(13*trs+9*trg+12*trc))+cDs**2*(2*(trs+15*trg+8*trc)*gDs+ms
     &   2*(13*trs+9*trg+12*trc))+8*(ms2*trc-mc2*trs)*cDs*gDs))-4*mc2*cD
     &   g*(cDs-mc*ms)*(cDs+mc*ms)*gDs**4*(-2*(-2*trs+trg+trc)*gDs+6*trs
     &   *cDs+trsgc-2*ms2*(-trs+trg+trc)-mc2*trg)-8*trc*cDg**6*(cDs-mc*m
     &   s)*(cDs+mc*ms)*(3*gDs+4*cDs-ms2)-16*trc*cDg**7*(cDs-mc*ms)*(cDs
     &   +mc*ms))*lVc/(cDg**3*(cDs-mc*ms)*(cDs+mc*ms)*gDs**3)/1.6e+1_dp

      Bmpp = Bmpp-xs*cDs*(-2*trs*cDg*gDs+trsgc*gDs+ms2*trg*cDg)*(mc2*g
     &   Ds**2-2*cDg*cDs*gDs+ms2*cDg**2)*lRcs/(ms*(xs**2-1)*cDg**2*gDs**
     &   2)

      Bmpp = (-mc2*cDg**4*(128*(2*trs+3*trc)*gDs**5+4*(-288*trsgc+ms2*(
     &   98*trs+66*trg+516*trc)+mc2*(-19*trs-25*trg+4*trc))*gDs**4+4*cDs
     &   *(2*gDs+ms2)*((-246*trs-34*trg+48*trc)*gDs**3+(mc2*(68*trs+50*t
     &   rg)-ms2*(41*trs+165*trg+84*trc))*gDs**2+ms2*(mc2*(68*(trs+trg)-
     &   9*trc)+ms2*(20*trs-16*trg-33*trc))*gDs+ms2**2*(9*trsgc+8*mc2*tr
     &   s+17*mc2*trg))+4*(mc2*(18*trsgc+ms2*(310*trs+25*trg+56*trc))-54
     &   0*ms2*trsgc+ms2**2*(30*trs-74*trg+340*trc))*gDs**3-36*cDs**2*(2
     &   *gDs+ms2)*(2*(2*trs+5*trg)*gDs**2+(2*ms2*(2*trg+trc)-2*mc2*trc)
     &   *gDs+ms2**2*(3*trs+2*trg))+ms2*(mc2*(180*trsgc+ms2*(1051*trs+48
     &   7*trg+108*trc))+2*ms2*(ms2*(-15*trs-119*trg+112*trc)-540*trsgc)
     &   )*gDs**2+2*ms2**2*(72*(mc-ms)*(ms+mc)*trsgc+(179*mc2-15*ms2)*ms
     &   2*(trs+trg))*gDs+4*mc2*ms2**3*(9*trsgc+19*ms2*(trs+trg)))+mc2*c
     &   Dg**3*gDs*(288*trg*gDs**5+2*cDs*(32*(10*trs+8*trg+23*trc)*gDs**
     &   4-4*(180*trsgc+ms2*(-37*trs+43*trg-128*trc)+mc2*(-267*trs+19*tr
     &   g+20*trc))*gDs**3+4*(mc2*(18*trsgc+ms2*(154*trs+73*trg-16*trc))
     &   +18*ms2*(ms2*(trc-2*trg)-5*trsgc))*gDs**2+mc2*ms2*(108*trsgc+ms
     &   2*(trs+197*trg-12*trc))*gDs+4*mc2*ms2**2*(9*trsgc-5*ms2*trs+4*m
     &   s2*trg))+16*(18*trsgc-2*mc2*(55*trs+6*trc)+ms2*(-3*trs+22*trg+1
     &   3*trc))*gDs**4-4*(72*ms2*trsgc+mc2*(ms2*(193*trs+309*trg+456*tr
     &   c)-144*trsgc)+2*ms2**2*(21*trs+5*(trg+trc)))*gDs**3+2*cDs**2*(2
     &   *gDs+ms2)*((-50*trs-230*trg+32*trc)*gDs**2+9*(-4*trsgc+ms2*trs+
     &   4*mc2*trs-3*ms2*trg+10*mc2*trg)*gDs+36*mc2*ms2*trg)-8*ms2*(9*ms
     &   2*(3*trsgc+ms2*(trs+trg+trc))+mc2*(ms2*(-3*trs+37*trg+136*trc)-
     &   153*trsgc)+36*mc2**2*trs)*gDs**2+mc2*ms2**2*(612*trsgc-94*mc2*t
     &   rs+ms2*(-3*trs+173*trg-112*trc)+50*mc2*trg)*gDs+mc2*ms2**3*(72*
     &   trsgc+5*(3*ms2+5*mc2)*(trs+trg)))-2*cDg**5*(4*(mc2*(-7*trs-25*t
     &   rg+4*trc)+12*ms2*trc)*gDs**4+4*cDs*(2*gDs+ms2)*(8*trc*gDs**3+(3
     &   6*trsgc+68*mc2*trs+50*mc2*trg-96*ms2*trc)*gDs**2+(36*ms2*trsgc+
     &   mc2*ms2*(68*(trs+trg)-27*trc)-33*ms2**2*trc+18*mc2**2*trc)*gDs+
     &   9*ms2**2*trsgc-mc2*ms2*(ms2*(19*trs+trg-3*trc)+3*mc2*trc))+4*mc
     &   2*(54*trsgc+ms2*(172*trs+25*trg-16*trc)-9*mc2*(2*trs+trg))*gDs*
     &   *3+mc2*ms2*(396*trsgc-6*mc2*(15*(2*trs+trg)+2*trc)+ms2*(863*trs
     &   +587*trg-72*trc))*gDs**2-12*cDs**2*(2*gDs+ms2)*(6*(ms**2-mc**2)
     &   *trc*gDs+3*ms2**2*(3*trs+2*trg)+2*mc2*ms2*trc-2*mc2**2*trc)+2*m
     &   c2*ms2**2*(126*trsgc+3*mc2*(trc-12*(2*trs+trg))+ms2*(204*(trs+t
     &   rg)+trc))*gDs+2*mc2*ms2**3*(27*trsgc+3*mc2*(trc-3*(2*trs+trg))+
     &   ms2*(38*(trs+trg)+5*trc)))+mc2*cDg**2*gDs**2*(-3*(2*gDs+ms2)*(-
     &   8*(6*trsgc-2*ms2*(-2*trs+trg+trc)+mc2*(-16*trs+23*trg+8*trc))*g
     &   Ds**3+4*(-mc2*(30*trsgc+ms2*(-5*trs+17*trg+5*trc))+58*mc2**2*tr
     &   s+4*ms2**2*(-trs+trg+trc))*gDs**2+4*mc2*ms2*(21*trsgc+3*ms2*(tr
     &   s+trg+trc)+mc2*(-5*trs+35*trg+17*trc))*gDs+3*mc2**2*ms2*(ms2*(t
     &   rs-3*trg)-4*trsgc))+12*cDs*gDs*(8*(3*(trsgc+ms2*trs)+mc2*(trs+1
     &   7*trg+16*trc))*gDs**2+(12*ms2*(trsgc+ms2*trs)+mc2*(-84*trsgc+3*
     &   ms2*trs+19*ms2*trg+76*ms2*trc)-12*mc2**2*(trg-4*trs))*gDs-6*mc2
     &   *ms2*(7*trsgc-4*mc2*trs+4*ms2*trg+mc2*trg-ms2*trc))+mc2*cDs**2*
     &   (2*gDs+ms2)*((-50*trs-230*trg+32*trc)*gDs-36*trsgc+9*ms2*(trs-3
     &   *trg)))+4*mc2*cDg*gDs**3*(2*gDs+ms2)*(-24*(mc-ms)*(ms+mc)*trs*g
     &   Ds**3+6*cDs*(2*(mc-ms)*(ms+mc)*(-3*trs+trg+trc)*gDs**2+(mc2*(6*
     &   trsgc+ms2*(4*(trg+trc)-trs))-3*ms2**2*(-trs+trg+trc)+mc2**2*(-4
     &   *trs+19*trg+10*trc))*gDs-3*mc2**2*(trsgc+ms2*trg))+6*(mc2*(6*tr
     &   sgc+ms2*(-4*trs+trg+trc))-2*ms2**2*(-2*trs+trg+trc)+mc2**2*(-12
     &   *trs+13*trg+4*trc))*gDs**2+12*(mc-ms)*(ms+mc)*(-trs+trg+trc)*cD
     &   s**2*gDs+(3*mc2**2*(18*trsgc+ms2*(2*trs+5*trg-7*trc))-26*mc2**3
     &   *trs-6*ms2**3*(-trs+trg+trc)+6*mc2*ms2**2*(-trs+trg+trc))*gDs-1
     &   8*mc2**2*ms2*(trsgc+mc2*trg))+4*cDg**6*(2*gDs+ms2)*(18*(-4*trsg
     &   c+mc2*(2*trs+trg)+8*ms2*trc)*gDs**2-12*cDs*(6*(mc**2-ms**2)*trc
     &   *gDs+ms2**2*(-9*trs-6*trg+trc)-3*mc2*ms2*trc+2*mc2**2*trc)-6*(-
     &   3*ms2*(ms2*trc-4*trsgc)-2*mc2*ms2*(6*trs+3*trg+2*trc)+3*mc2**2*
     &   trc)*gDs-24*(mc**2-ms**2)*trc*cDs**2+ms2*(-2*ms2*(9*trsgc+5*ms2
     &   *trc)+3*mc2*ms2*(21*trs+12*trg-4*trc)+6*mc2**2*trc))+4*mc2**2*g
     &   Ds**4*(2*gDs+ms2)*(-12*(mc-ms)*(ms+mc)*trs*gDs**2+9*mc2*trsgc*(
     &   gDs+cDs+mc2)+6*(mc-ms)*(ms+mc)*(-3*trs+trg+trc)*cDs*gDs-6*ms2**
     &   2*(-2*trs+trg+trc)*gDs+6*mc2*ms2*(-3*trs+trg+trc)*gDs+3*mc2**2*
     &   (3*trg-4*trs)*gDs+6*(mc-ms)*(ms+mc)*(-trs+trg+trc)*cDs**2-9*ms2
     &   **2*(-trs+trg+trc)*cDs+12*mc2*ms2*(-trs+trg+trc)*cDs+mc2**2*(-5
     &   *trs+14*trg+5*trc)*cDs+3*ms2**3*trs-6*mc2*ms2**2*trs+6*mc2**2*m
     &   s2*trs+5*mc2**3*trs-3*ms2**3*trg+6*mc2*ms2**2*trg-2*mc2**2*ms2*
     &   trg-3*ms2**3*trc+6*mc2*ms2**2*trc-11*mc2**2*ms2*trc)+24*cDg**7*
     &   (2*gDs+ms2)*(6*(ms**2-mc**2)*trc*gDs+8*(ms**2-mc**2)*trc*cDs+9*
     &   ms2**2*trs+6*ms2**2*trg-2*ms2**2*trc+4*mc2*ms2*trc-2*mc2**2*trc
     &   )-96*(mc-ms)*(ms+mc)*trc*cDg**8*(2*gDs+ms2))/(mc*cDg**3*(2*cDg+
     &   mc2)*gDs**3*(2*gDs+ms2))/4.8e+1_dp+Bmpp

      Bmpp = LsB1*mc*(ms2*cDs*(6*mc2*trc*gDs**5-3*mc2*(trsgc+mc2*trg)*g
     &   Ds**4+3*cDg*(mc2*(trsgc+3*ms2*trg)-2*(mc2*trs+ms2*trc)*cDg)*gDs
     &   **3+3*ms2*(trsgc-mc2*trg)*cDg**2*gDs**2+ms2*cDg**3*(-2*trs*cDg+
     &   trsgc-ms2*trg)*gDs+2*ms2**2*trg*cDg**4)+ms2*gDs*(mc2**2*(ms2*tr
     &   g-trsgc)*gDs**3+ms2*cDg**3*(2*mc2*trs*gDs+2*ms2*trc*gDs-ms2*trs
     &   gc+mc2*ms2*trg)+mc2*cDg*gDs**2*(2*mc2*trs*gDs-6*ms2*trc*gDs+3*m
     &   s2*trsgc+mc2*ms2*trg)-mc2*ms2*(trsgc+3*ms2*trg)*cDg**2*gDs)+2*m
     &   s2*cDs**2*gDs*(-4*mc2*trg*gDs**3+cDg**3*(2*trs*gDs-4*ms2*trg)+c
     &   Dg*gDs**2*(6*trc*gDs-3*trsgc+mc2*trg)+(3*ms2*trg-trsgc)*cDg**2*
     &   gDs)+4*cDs**3*gDs**2*(-2*trc*gDs**3+(trsgc+mc2*trg)*gDs**2-3*ms
     &   2*trg*cDg*gDs+3*ms2*trg*cDg**2)+8*trg*cDs**4*gDs**3*(gDs-cDg))/
     &   (cDg*(cDs**2-mc2*ms2)*gDs**3)+Bmpp

      Bmpp = ls*(-16*mc2*trc*gDs**5+4*mc2*(-3*trs*cDs+7*trg*cDs+(trs+3*
     &   trg)*cDg+4*mc2*trg-6*ms2*trc)*gDs**4+4*(mc2*(ms2*(11*trg*cDs+2*
     &   trsgc+4*mc2*trg)-2*trsgc*cDs-2*ms2**2*trc)+mc2*cDg*(-(3*trg-2*t
     &   rs)*(2*cDs+ms2)-2*trsgc)+2*(-4*trsgc+mc2*(2*trs+trg)+8*ms2*trc)
     &   *cDg**2)*gDs**3+ms2*(mc2*(-4*trsgc*cDs+ms2*(trs+13*trg)*cDs+4*m
     &   s2*(trsgc+mc2*trg))-mc2*cDg*(-8*trs*cDs+28*trg*cDs+4*trsgc+9*ms
     &   2*trs+37*ms2*trg)+4*(-12*trsgc+2*mc2*trs+5*mc2*trg+10*ms2*trc)*
     &   cDg**2)*gDs**2+8*ms2**2*cDg*(cDg*(2*trg*(cDs+mc2)+3*trs*cDs-3*t
     &   rsgc+ms2*trc)-mc2*trg*(cDs+ms2)+(3*trs+2*trg)*cDg**2)*gDs+4*ms2
     &   **3*cDg**2*(trg*(2*(cDs+cDg)+mc2)+3*trs*(cDs+cDg)-trsgc))/(mc*c
     &   Dg*gDs*(2*gDs+ms2)**2)/8.0_dp+Bmpp

      Bmpp = BfunX*(-cDg**4*(-16*trc*gDs**4+4*cDs*(-4*trc*gDs**3+4*mc2*
     &   (trs+trg)*gDs**2+mc2*ms2*(4*(trs+trg)-3*trc)*gDs+mc2*ms2**2*(tr
     &   s+trg))-4*mc2*(trs-trg+trc)*gDs**3+4*mc2*ms2*(2*trs+3*trg+4*trc
     &   )*gDs**2+24*trc*cDs**2*gDs*(2*gDs+mc2)+48*trc*cDs**3*gDs+13*mc2
     &   *ms2**2*(trs+trg)*gDs+4*mc2*ms2**3*(trs+trg))+mc2*cDg**3*gDs*(4
     &   *(3*trs+trg+trc)*gDs**3+4*(2*trs*cDs+ms2*(4*trs+trg+trc))*gDs**
     &   2-ms2*(8*trs*cDs+12*(trg+trc)*cDs+ms2*(5*trs+9*trg+4*trc))*gDs-
     &   5*ms2**2*(trs+trg)*(cDs+ms2))-4*trc*cDg**5*(4*(-gDs**3+cDs**2*(
     &   11*gDs+mc2)+3*cDs*gDs*(2*gDs+mc2)+2*cDs**3)-mc2*ms2*(gDs+2*cDs)
     &   +mc2*ms2**2)-8*trc*cDg**6*(6*gDs**2+26*cDs*gDs+3*mc2*gDs+12*cDs
     &   **2+4*mc2*cDs-mc2*ms2)+8*mc2*cDg*gDs**4*(2*(gDs+cDs)+ms2)*((trg
     &   +trc)*(gDs+cDs+2*ms2)-2*trs*(2*(gDs+cDs)+ms2))+4*mc2*cDg**2*gDs
     &   **3*(cDs*(-5*trs*gDs+7*(trg+trc)*gDs+6*ms2*(trs+trg+trc))+3*(2*
     &   gDs+ms2)*((-trs+trg+trc)*gDs+ms2*(trs+trg+trc))+(trs+trg+trc)*c
     &   Ds**2)+4*mc2*gDs**4*(2*(gDs+cDs)+ms2)**2*((trg+trc)*(cDs+ms2)-t
     &   rs*(2*gDs+cDs+ms2))-16*trc*cDg**7*(5*gDs+6*cDs+mc2)-32*trc*cDg*
     &   *8)/(mc*cDg**3*gDs**3)/1.6e+1_dp+Bmpp

      Bmpp = LsB2*mc*(2*mc2**2*trsgc*cDs*gDs**4+cDg**3*gDs*(cDs*(-4*mc2
     &   *trs*gDs**2+2*ms2*(3*trsgc+mc2*(-5*trs+trg-3*trc))*gDs+mc2*ms2*
     &   (3*trsgc+5*ms2*trg))+cDs**2*(8*trs*gDs**2-4*(trsgc-2*mc2*trs+ms
     &   2*trg)*gDs+8*mc2*ms2*trg)+mc2*ms2*(4*trs*gDs**2-2*(trsgc-2*mc2*
     &   trs+ms2*trg)*gDs+ms2*(trsgc-5*mc2*trg))+cDs**3*(4*ms2*trg-8*(tr
     &   s+trg)*gDs))-cDg**2*gDs**2*(mc2*cDs*(4*trs*gDs**2-2*(trsgc-2*mc
     &   2*trs+ms2*trg)*gDs+ms2*(trsgc-5*mc2*trg))+mc2*ms2*(2*(trsgc+mc2
     &   *(-trs+trg+trc))*gDs+mc2*(ms2*trg-trsgc))-4*cDs**2*((mc2*(4*trs
     &   +trg+trc)-trsgc)*gDs-mc2*(trsgc+2*ms2*trg))+(8*mc2*trg-4*trsgc)
     &   *cDs**3)+cDg**4*(4*trs*(2*cDs**2-3*ms2*cDs+mc2*ms2)*gDs**2+2*ms
     &   2*(2*(trs+2*trg)*cDs**2+mc2*(ms2*(trs-trg+3*trc)-6*trs*cDs)+3*(
     &   trsgc+ms2*trg)*cDs-ms2*trsgc)*gDs-ms2*(2*(ms2*trg-trsgc)*cDs**2
     &   +ms2*(trsgc+mc2*trg)*cDs+mc2*ms2*(3*trsgc+ms2*trg)))+mc2*cDg*gD
     &   s**3*(cDs*(mc2*(-6*trs*gDs+trsgc+3*ms2*trg)+2*trsgc*gDs)+2*(mc2
     &   *trg-3*trsgc)*cDs**2+mc2*ms2*(trsgc-mc2*trg))-2*ms2*cDg**5*(cDs
     &   *(6*trs*gDs+ms2*(trg+trc))+ms2*(-2*trs*gDs+trsgc-2*mc2*trs+ms2*
     &   trg))+4*ms2**2*trs*cDg**6)/(cDg**3*(cDs**2-mc**2*ms**2)*gDs)+B
     &   mpp

      Bmpp = tr3c00fs*(cDg**2*(4*cDs*(-16*trc*gDs**3+(4*trsgc+4*mc2*trs
     &   +2*mc2*trg-36*ms2*trc)*gDs**2+ms2*(4*(trsgc+mc2*(trs+trg))-11*m
     &   s2*trc)*gDs+ms2**2*(trsgc+mc2*trg))+mc2*(4*trsgc*(gDs+ms2)**2-g
     &   Ds*(2*(5*trs+7*trg)*gDs**2+ms2*(21*trs+17*trg)*gDs-ms2**2*(15*t
     &   rs+19*trg)))-4*ms2**2*(3*trs+2*trg)*cDs**2)+mc2*cDg*gDs*(12*trc
     &   *gDs**3+4*(-7*trs*cDs-9*trg*cDs-8*trsgc+11*ms2*trc)*gDs**2+2*(-
     &   ms2*(19*trs*cDs+31*trg*cDs+18*trsgc)-2*cDs*(2*trs*cDs+5*trg*cDs
     &   +2*trsgc)+8*ms2**2*trc)*gDs+ms2*(cDs*(-8*trg*cDs+15*ms2*trs+7*m
     &   s2*trg)-8*trsgc*(cDs+ms2)))+mc2*gDs**2*(4*cDs*(-7*trc*gDs**2+(4
     &   *trsgc+2*mc2*trg-2*ms2*trc)*gDs+ms2*(trsgc+mc2*trg))+cDs**2*(14
     &   *trs*gDs+26*trg*gDs+4*trsgc-ms2*trs+3*ms2*trg)+mc2*ms2**2*(trs+
     &   trg))-4*cDg**3*((mc2*(2*trs+trg)+8*ms2*trc)*gDs**2+ms2*(2*mc2*(
     &   2*trs+trg)-5*ms2*trc)*gDs+ms2**2*(6*trs*cDs+4*trg*cDs+2*mc2*trs
     &   +mc2*trg-3*ms2*trc))+4*trsgc*cDg**3*(2*gDs+ms2)**2-4*ms2**2*(3*
     &   trs+2*trg)*cDg**4)/(mc*cDg*gDs**3)/4.0_dp+Bmpp

      Bmpp = 3.0_dp*tr3c001fs*(-cDg**2*(mc2*(trs+trg)*(4*gDs**3+12*ms2*
     &   gDs**2+7*ms2**2*gDs+4*ms2**3)+4*trc*cDs*gDs*(8*gDs**2+12*ms2*gD
     &   s+3*ms2**2))+4*mc2*(trs+trg)*cDs**2*gDs**3-mc2*(trs+trg)*cDg*cD
     &   s*gDs*(16*gDs**2+12*ms2*gDs-5*ms2**2)+4*ms2**2*trc*cDg**3*(3*gD
     &   s+ms2))/(4.0_dp*mc*cDg*gDs**3)+Bmpp

      Bmpp = B0cgsf*(-cDg**4*(-2*(-8*trsgc+mc2*(trs+5*trg+2*trc)+12*ms2
     &   *trc)*gDs**3+4*cDs*((4*trsgc+8*mc2*trs+6*mc2*trg-12*ms2*trc)*gD
     &   s**2+ms2*(4*trsgc+mc2*(8*(trs+trg)-3*trc)-ms2*(3*trs+2*(trg+trc
     &   )))*gDs+ms2**2*(trsgc+mc2*(trs+2*trg)))+(mc2*(4*trsgc+ms2*(7*tr
     &   s+15*trg+8*trc))-8*ms2*(ms2*trc-2*trsgc))*gDs**2+4*cDs**2*(6*(m
     &   c-ms)*(ms+mc)*trc*gDs-ms2**2*(3*trs+2*trg))+4*ms2*(ms2+2*mc2)*t
     &   rsgc*gDs+mc2*ms2**2*(17*trs+25*trg)*gDs+4*mc2*ms2**2*(trsgc+2*m
     &   s2*(trs+trg)))+mc2*cDg**3*gDs*(-2*(3*trs-7*trg+6*trc)*gDs**3+(4
     &   *(-trs*cDs+8*trg*cDs+5*trsgc+2*mc2*trs+mc2*trg)-ms2*(7*(trs+trg
     &   )+20*trc))*gDs**2+(4*(2*trs+5*trg)*cDs**2+8*trsgc*cDs-2*ms2*(2*
     &   trc-5*(trs+3*trg))*cDs-ms2*(-24*trsgc+ms2*(5*trs+trg+12*trc)+4*
     &   mc2*(4*trs+trg)))*gDs+ms2*(8*trg*cDs**2+(8*trsgc-5*ms2*trs+3*ms
     &   2*trg)*cDs+ms2*(8*trsgc+5*ms2*(trs+trg))))+mc2*cDg**2*gDs**2*(-
     &   8*(trc-trg)*gDs**3-2*(-8*trsgc+2*ms2*(trs-trg+trc)+7*mc2*trs-11
     &   *mc2*trg)*gDs**2+cDs*gDs*(2*(5*(trs+3*trg)+6*trc)*gDs+4*trsgc+8
     &   *mc2*(4*trs+trg)+ms2*(trs-3*trg+8*trc))+cDs**2*((-6*trs+6*trg+4
     &   *trc)*gDs-4*trsgc+ms2*(trs-3*trg))-(mc2*(ms2*(-17*trs+23*trg+16
     &   *trc)-4*trsgc)+12*ms2**2*(trs+trg+trc))*gDs+mc2*ms2*(4*trsgc-ms
     &   2*trs+3*ms2*trg))+4*mc2*cDg*gDs**3*(2*(trsgc+2*mc2*(-2*trs+2*tr
     &   g+trc)-ms2*(-2*trs+trg+trc))*gDs**2+2*cDs*((trsgc+3*ms2*trs+mc2
     &   *(-3*trs+5*trg+3*trc))*gDs-mc2*(trsgc+ms2*trg))+(mc2*(trsgc+2*m
     &   s2*(-trs+2*trg+trc))-2*ms2**2*(-trs+trg+trc)+mc2**2*(trg-4*trs)
     &   )*gDs-2*mc2*ms2*(trsgc+mc2*trg))+4*cDg**5*((-4*trsgc+2*mc2*trs+
     &   mc2*trg+8*ms2*trc)*gDs**2+2*cDs*(6*(ms**2-mc**2)*trc*gDs+ms2*(3
     &   *ms2*trs+2*ms2*trg-ms2*trc+mc2*trc))+ms2*(-4*trsgc+mc2*(4*trs+2
     &   *trg+trc)+3*ms2*trs+2*ms2*trg)*gDs+4*(ms2-mc2)*trc*cDs**2-ms2**
     &   2*(trsgc+mc2*(-2*trs-trg+trc)))+4*mc2*gDs**4*(-4*(mc-ms)*(ms+mc
     &   )*trs*gDs**2+cDs*(2*(mc-ms)*(ms+mc)*(-3*trs+trg+trc)*gDs+mc2*(t
     &   rsgc+4*ms2*(-trs+trg+trc))-3*ms2**2*(-trs+trg+trc)+mc2**2*trg)+
     &   (mc2*(trsgc+2*ms2*(-3*trs+trg+trc))-2*ms2**2*(-2*trs+trg+trc)+m
     &   c2**2*trg)*gDs+2*(mc-ms)*(ms+mc)*(-trs+trg+trc)*cDs**2+mc2**2*t
     &   rsgc-ms2*(ms2**2*(-trs+trg+trc)-2*mc2*ms2*(-trs+trg+trc)+mc2**2
     &   *(-2*trs+trg+2*trc)))+4*cDg**6*(6*(ms**2-mc**2)*trc*gDs+8*(ms**
     &   2-mc**2)*trc*cDs+ms2*(3*ms2*trs+2*ms2*trg-2*ms2*trc+2*mc2*trc))
     &   -16*(mc-ms)*(ms+mc)*trc*cDg**7)/(mc*cDg**3*gDs**3)/1.6e+1_dp+Bmpp

      Bmpp = mc*tr3s00ft*(-cDg**2*(2*trg*gDs**3+2*(-2*trs*cDs+trsgc-5*m
     &   c2*trs+3*ms2*(trg+trc))*gDs**2-(2*trg*cDs**2+mc2*(ms2*(-7*trs+7
     &   *trg+4*trc)-8*trs*cDs)+4*(trsgc+ms2*trg)*cDs+3*ms2*trsgc)*gDs+m
     &   s2*(mc2*trg*(cDs+ms2)+trsgc*(cDs+mc2)))+cDg*gDs*(-2*(trsgc+mc2*
     &   (-2*trs+3*trg+trc))*gDs**2+2*cDs*(mc2*(trsgc+ms2*trg)-(trsgc+mc
     &   2*(-4*trs+3*trg+trc))*gDs)+mc2*(-3*trsgc+7*mc2*trs+ms2*(trg+3*t
     &   rc))*gDs+2*mc2*ms2*(trsgc+mc2*trg))+mc2*gDs**2*(mc2*(trs*(2*gDs
     &   +3*cDs)-trg*(gDs+cDs)+ms2*(2*trg+3*trc))-trsgc*(gDs+cDs+mc2))-c
     &   Dg**3*(ms2*(5*trs*gDs+4*(trg+trc)*gDs+trsgc-4*mc2*trs+ms2*(5*tr
     &   g+3*trc))+trs*cDs*(8*gDs-3*ms2)+2*ms2*(trg+trc)*cDs)+7*ms2*trs*
     &   cDg**4)/(cDg**3*gDs)+Bmpp

      Bmpp = (-3.0_dp)*mc*tr3c002fs*(4*cDg**2*(trc*gDs**3+(4*(trs+trg)*
     &   cDs+3*ms2*trc)*gDs**2+4*ms2*(trs+trg)*cDs*gDs+ms2**2*(trs+trg)*
     &   cDs)-4*trc*cDs**2*gDs**3+cDg*gDs*(4*trc*cDs*gDs*(4*gDs+3*ms2)-5
     &   *mc2*ms2**2*(trs+trg)))/(4.0_dp*cDg*gDs**3)+Bmpp

      Bmpp = Bmpp-B0csf*mc*(cDs*(2*cDg*(mc2*trs*gDs**2+ms2*(trsgc+mc2*
     &   (trc-trg)+ms2*trc)*gDs+mc2*ms2**2*trg)+mc2*gDs*((-trsgc+ms2*(tr
     &   s-trg+2*trc)+mc2*trs)*gDs-2*ms2*(trsgc+(ms2+mc2)*trg))-ms2*cDg*
     &   *2*(4*trs*gDs+trsgc-mc2*trs+ms2*(-trs+trg+2*trc))+2*ms2*trs*cDg
     &   **3)+ms2*(mc2*(2*trs*cDg-trsgc+mc2*(trg+trc)+ms2*trc)*gDs**2-mc
     &   2*(ms2+mc2)*trsgc*gDs+ms2*cDg*(2*trs*cDg**2-(trsgc-mc2*trg+ms2*
     &   trc+mc2*trc)*cDg+mc2*(ms2+mc2)*trg))+cDs**2*(-2*(mc2*(-trs+trg+
     &   trc)+ms2*trc)*gDs**2+2*ms2*((trg+2*trc)*cDg-2*mc2*trg)*gDs-2*tr
     &   s*cDg*(2*cDg+ms2+mc2)*gDs+trsgc*(2*cDg+ms2+mc2)*gDs-ms2*cDg*(tr
     &   g*(2*cDg+ms2+mc2)-2*trs*cDg))+2*cDs**3*(-2*trc*gDs**2+(trg*(2*c
     &   Dg+ms2+mc2)-2*trs*cDg+trsgc)*gDs-ms2*trg*cDg)+4*trg*cDs**4*gDs)
     &   /(cDg*(cDs-mc*ms)*(cDs+mc*ms)*gDs)/2.0_dp

      Bmpp = 3*mc*tr3s001ft*(trs*((2*cDg+mc2)**2*cDs*gDs**2+ms2*cDg**2*
     &   (cDg*(cDs-3*gDs)-3*mc2*gDs))-ms2*(trg+trc)*((cDg-mc2)*(2*cDg+mc
     &   2)*gDs**2+ms2*cDg**3))/(cDg**3*gDs)+Bmpp

      Bmpp = lc*mc*(cDg**2*((4*trsgc+6*mc2*trg)*gDs**2+26*mc2*trg*cDs*g
     &   Ds+6*mc2*trsgc*gDs+mc2*ms2*(trsgc-3*mc2*trg))+mc2**2*(trsgc-mc2
     &   *trg)*gDs**2+2*cDg**3*(mc2*(ms2*(trc-5*trg)-2*trs*gDs)+trsgc*(2
     &   *gDs+ms2)+2*trg*gDs*(gDs+5*cDs))+2*mc2*cDg*gDs*(trsgc*(2*gDs+mc
     &   2)+4*mc2*trg*cDs)+4*cDg**4*(ms2*(trc-2*trg)-trs*gDs))/(cDg*(2*c
     &   Dg+mc2)**2*gDs)/2.0_dp+Bmpp

      Bmpp = epinv*(mc*ms*(xs**2-1)-4*lp*xs*cDs)*(-2*trs*cDg*gDs+trsgc*
     &   gDs+ms2*trg*cDg)*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)/(ms*(xs*
     &   *2-1)*cDg**2*gDs**2)/4.0_dp+Bmpp

      Bmpp = Bmpp-mc*tr2fu*cDs*(-2*trs*cDg*gDs+trsgc*gDs+ms2*trg*cDg)*
     &   (mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)/(cDg**2*gDs**2)

      Bmpp = 3*mc*tr3s002ft*((trg+trc)*(2*cDg+mc2)*((2*cDg+mc2)*cDs*gDs
     &   -3*ms2*cDg**2)+trs*(3*mc2**2*cDg*gDs+mc2**3*gDs+2*cDg**3*(cDs-g
     &   Ds)))/cDg**3+Bmpp

      Bmpp=Bmpp/zmpp

      return
      end
