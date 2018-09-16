      subroutine Bamp_pmm(q,mc,ms,Bpmm)
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
      complex(dp):: trc,trg,trs,trsgc,zpmm,Bpmm

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

      Bpmm = ms*(gDs+cDg)*(mc2*gDs**3*(cDs*(4*trc*gDs**2-2*(trsgc+mc2*t
     &   rg-2*ms2*trc)*gDs-ms2*(trsgc+3*mc2*trg))+mc2*ms2*(2*(trg-trs)*g
     &   Ds-trsgc+3*ms2*trg)-4*trg*cDs**2*(gDs+ms2))-cDg**2*gDs*(2*cDs**
     &   2*(4*trc*gDs**2-2*trsgc*gDs+2*ms2*(trs+trg+3*trc)*gDs-2*ms2*trs
     &   gc+ms2**2*trg)+ms2*cDs*(-4*trc*gDs**2+2*(trsgc+mc2*trg-2*ms2*tr
     &   c)*gDs+ms2*(trsgc+3*mc2*trg))+mc2*ms2**2*(-2*(trs+trg)*gDs+trsg
     &   c-ms2*trg))+2*cDg*gDs**2*(cDs*(mc2*((2*gDs+3*ms2)*(trc*gDs-ms2*
     &   trg)+2*ms2*trs*gDs)+mc2*trg*cDs*(2*gDs+3*ms2)-4*trc*cDs*gDs*(gD
     &   s+ms2)+4*trg*cDs**2*(gDs+ms2))+trsgc*(2*cDs**2*(gDs+ms2)-mc2*cD
     &   s*gDs-mc2*ms2**2))+2*ms2*cDg**3*cDs*(trc*gDs*(2*gDs+3*ms2)-trsg
     &   c*(gDs+ms2)))*tr5Xs/(cDg*(cDs-mc*ms)*(cDs+mc*ms)*gDs**3)/2.0_dp

      Bpmm = Bpmm-mc2*ms*(gDs+cDg)*(2*mc2**2*trg*cDs*gDs**3+cDg**3*(cD
     &   s*(-4*mc2*trs*gDs-2*ms2*trc*gDs+ms2*trsgc+3*mc2*ms2*trg)-4*trg*
     &   cDs**3-2*(trsgc+ms2*trg)*cDs**2+mc2*ms2*(trsgc+ms2*trg))+2*cDg*
     &   *2*gDs*(mc2**2*(trs*gDs-ms2*trg)+cDs**2*(2*trc*gDs-trsgc+2*mc2*
     &   trg)+2*trg*cDs**3+mc2*trsgc*cDs)-mc2*cDg*gDs**2*(cDs*(2*trc*gDs
     &   -trsgc+mc2*trg)+6*trg*cDs**2+mc2*(trsgc-ms2*trg))+cDg**4*(4*trs
     &   *cDs**2-2*mc2*ms2*trs))*tr5Xc/(cDg**3*(cDs**2-mc2*ms2)*gDs)/2.0
     &   _dp

      Bpmm = ms*(cDg**2*gDs*(2*cDs**2*(4*trc*gDs**2-2*trsgc*gDs+2*ms2*(
     &   trs+trg+3*trc)*gDs-2*ms2*trsgc+ms2**2*trg)+ms2*cDs*(-4*trc*gDs*
     &   *2+2*(trsgc+mc2*trg-2*ms2*trc)*gDs+ms2*(trsgc+3*mc2*trg))+mc2*m
     &   s2**2*(-2*(trs+trg)*gDs+trsgc-ms2*trg))+mc2*gDs**3*(cDs*(-4*trc
     &   *gDs**2+2*(trsgc+mc2*trg-2*ms2*trc)*gDs+ms2*(trsgc+3*mc2*trg))+
     &   mc2*ms2*(2*trs*gDs-2*trg*gDs+trsgc-3*ms2*trg)+4*trg*cDs**2*(gDs
     &   +ms2))+2*cDg*gDs**2*(trsgc*(-2*cDs**2*(gDs+ms2)+mc2*cDs*gDs+mc2
     &   *ms2**2)-cDs*(mc2*((2*gDs+3*ms2)*(trc*gDs-ms2*trg)+2*ms2*trs*gD
     &   s)+mc2*trg*cDs*(2*gDs+3*ms2)-4*trc*cDs*gDs*(gDs+ms2)+4*trg*cDs*
     &   *2*(gDs+ms2)))+2*ms2*cDg**3*cDs*(trsgc*(gDs+ms2)-trc*gDs*(2*gDs
     &   +3*ms2)))*tr4Xs/(cDg*(cDs-mc*ms)*(cDs+mc*ms)*gDs**2)/2.0_dp+Bp
     &   mm

      Bpmm = mc2*ms*(2*mc2**2*trg*cDs*gDs**3+cDg**3*(cDs*(-4*mc2*trs*gD
     &   s-2*ms2*trc*gDs+ms2*trsgc+3*mc2*ms2*trg)-4*trg*cDs**3-2*(trsgc+
     &   ms2*trg)*cDs**2+mc2*ms2*(trsgc+ms2*trg))+2*cDg**2*gDs*(mc2**2*(
     &   trs*gDs-ms2*trg)+cDs**2*(2*trc*gDs-trsgc+2*mc2*trg)+2*trg*cDs**
     &   3+mc2*trsgc*cDs)-mc2*cDg*gDs**2*(cDs*(2*trc*gDs-trsgc+mc2*trg)+
     &   6*trg*cDs**2+mc2*(trsgc-ms2*trg))+cDg**4*(4*trs*cDs**2-2*mc2*ms
     &   2*trs))*tr4Xc/(cDg**2*(cDs**2-mc2*ms2)*gDs)/2.0_dp+Bpmm

      Bpmm = Bpmm-ms*(-2*mc2**3*ms2*trg*gDs**3-2*cDg**2*gDs*(mc2*ms2*(
     &   -13*(trg+trc)*gDs**2+2*(trsgc+mc2*trg)*gDs+mc2*(trsgc+2*ms2*trg
     &   ))+mc2*cDs*(mc2*trs*gDs+2*ms2*trc*gDs-ms2*trsgc+mc2*ms2*trg)+cD
     &   s**2*gDs*(13*(trg+trc)*gDs-2*(trsgc+mc2*trg)))+cDg**3*(6*(2*trg
     &   -3*trs)*(cDs**2-mc**2*ms**2)*gDs**2+2*((trsgc+ms2*trg)*cDs**2+m
     &   c2*ms2*(-trsgc-ms2*trg+ms2*trc)+2*mc2**2*ms2*trs)*gDs+mc2*ms2*(
     &   trsgc*cDs+ms2*trg*cDs-ms2*trsgc+mc2*ms2*trg))-2*cDg**4*(2*mc2*m
     &   s2*(trs*gDs-ms2*trg)+2*cDs**2*(ms2*trg-trs*gDs)+mc2*ms2*trs*cDs
     &   )+mc2*cDg*gDs**2*(2*(mc2*ms2*(trg+trc)-trg*cDs**2)*gDs+mc2*(trs
     &   gc*cDs+5*ms2*trg*cDs-ms2*trsgc+mc2*ms2*trg)))*tr3Xs/(cDg**2*(cD
     &   s**2-mc**2*ms**2)*gDs)/2.0_dp

      Bpmm = ms*(2*cDg**2*gDs*(mc2*cDs*(8*trc*gDs**2-4*trsgc*gDs+12*ms2
     &   *trc*gDs-3*ms2*trsgc+ms2**2*trg)+mc2*ms2*(3*trc*gDs**2-(trsgc+1
     &   3*ms2*trc)*gDs+ms2*(trsgc+3*mc2*trg))+cDs**2*gDs*(-7*trc*gDs+3*
     &   trsgc+2*mc2*trg+9*ms2*trc)+2*(trs+trg)*cDs**3*gDs)-2*mc2*cDg*gD
     &   s**2*(mc2*(4*trc*gDs**2+6*ms2*trc*gDs+trsgc*(-2*gDs-ms2)+ms2**2
     &   *trg)+cDs*(-8*trc*gDs**2+4*(trsgc+mc2*trg-2*ms2*trc)*gDs+2*ms2*
     &   (trsgc+3*mc2*trg))+cDs**2*(4*(trs+trg)*gDs+trsgc+ms2*trg))+mc2*
     &   gDs**3*(mc2*(-8*trc*gDs**2+(4*trsgc+4*mc2*trg-11*ms2*trc)*gDs+m
     &   s2*(2*trsgc+6*mc2*trg+3*ms2*trc))+2*mc2*cDs*(2*(trs+trg)*gDs+tr
     &   sgc+ms2*trg)+3*trc*cDs**2*(gDs-ms2))+4*cDg**3*(mc2*ms2*(gDs*(2*
     &   (trs+trg)*gDs+ms2*(trg-3*trc))+trsgc*(2*gDs+ms2))-cDs**2*gDs*(2
     &   *(trs+trg+trc)*gDs+trsgc+ms2*trg))-8*trg*cDg**4*(cDs-mc*ms)*(cD
     &   s+mc*ms)*gDs)*tr3Xc/(cDg*(cDs-mc*ms)*(cDs+mc*ms)*gDs**2)/4.0_dp
     &   +Bpmm

      Bpmm = ms*(-mc2*cDg**3*gDs*(mc2*ms2*((36*trs-24*trg+5*trc)*gDs**2
     &   +(ms2*(2*trs-2*trg+trc)+8*mc2*trs)*gDs+2*ms2**2*(trs+trg))-cDs*
     &   *2*(3*(12*trs-8*trg+7*trc)*gDs**2-8*trsgc*gDs+ms2*(2*trs-2*trg+
     &   9*trc)*gDs+2*ms2**2*(trs+trg))-2*mc2*ms2*cDs*(4*(-trs-trg+trc)*
     &   gDs-3*trsgc+ms2*(trs+3*trg))+2*(ms2*(trs+3*trg)-3*trsgc)*cDs**3
     &   )+2*cDg**4*(mc2**2*ms2*(8*trs*gDs**2-(trsgc+ms2*(trs+7*trg+2*tr
     &   c))*gDs+ms2*(trsgc+ms2*trs))+mc2*cDs**2*(ms2*(trs*gDs+7*trg*gDs
     &   -trsgc)+gDs*(trsgc-8*trs*gDs)-ms2**2*(trs+3*trc))+cDs**4*(2*(tr
     &   s+trg)*gDs+3*ms2*trc)+mc2*ms2*cDs*gDs*(3*trc*gDs-trsgc+2*mc2*tr
     &   s-2*mc2*trg-11*ms2*trc)+cDs**3*gDs*(-7*trc*gDs+3*trsgc+2*mc2*tr
     &   g+9*ms2*trc))+mc2*cDg**2*gDs**2*(cDs**2*(52*(trg+trc)*gDs**2+(-
     &   8*trsgc+ms2*(2*trs-trg+trc)+4*mc2*trs-8*mc2*trg)*gDs-2*(3*mc2*t
     &   rsgc+ms2**2*(trs+trg)-mc2*ms2*trs))+mc2*ms2*(-52*(trg+trc)*gDs*
     &   *2+(8*trsgc+ms2*(-2*trs+trg-trc)-4*mc2*(trc-3*trg))*gDs+2*(3*mc
     &   2*trsgc+ms2**2*(trs+trg)-mc2*ms2*trs))-mc2*cDs*(8*trc*gDs**2+(5
     &   *ms2*trc-4*(trsgc+mc2*trs))*gDs+ms2**2*(2*(trs+trg)+3*trc))+cDs
     &   **3*(trc*gDs+ms2*(2*(trs+trg)+3*trc)))+2*mc2*(cDs-mc*ms)*(cDs+m
     &   c*ms)*gDs**3*(trs*(cDs+ms2)*(4*gDs**2+ms2*(4*(gDs+cDs)-mc2)+8*c
     &   Ds*gDs-2*mc2*gDs+4*cDs**2-2*mc2*cDs+ms2**2+mc2**2)+ms2*(trg+trc
     &   )*gDs*(2*gDs+2*cDs+ms2-mc2))-2*cDg**5*(cDs-mc*ms)*(cDs+mc*ms)*(
     &   cDs*(2*(trs+trg+2*trc)*gDs+2*trsgc+2*ms2*trg-3*ms2*trc)+gDs*(7*
     &   trc*gDs-3*trsgc-2*mc2*trg-9*ms2*trc))+2*mc2*cDg*(cDs-mc*ms)*(cD
     &   s+mc*ms)*gDs**3*((ms2*(4*trs+trg+trc)+2*mc2*trg)*gDs+2*(2*ms2+m
     &   c2)*trs*cDs+3*ms2*(trg+trc)*cDs+(2*ms2**2+mc2**2)*trs)-4*cDg**6
     &   *(cDs-mc*ms)*(cDs+mc*ms)*(2*(trs+trg+trc)*gDs+2*trg*cDs+trsgc+m
     &   s2*trg)-8*trg*cDg**7*(cDs-mc*ms)*(cDs+mc*ms))*lVs/(mc2*cDg**3*(
     &   mc2*ms2-cDs**2)*gDs**2)/8.0_dp+Bpmm

      Bpmm = (-2*cDg**4*(cDs**2*(4*trs*gDs**3+ms2*(6*trc-5*trg)*gDs**2-
     &   ms2*(3*trsgc+2*mc2*(trs+trg)+ms2*(2*trg+13*trc))*gDs+ms2**2*(2*
     &   trsgc-mc2*(trs+trg+3*trc)))+mc2*ms2*(-4*trs*gDs**3+ms2*(5*trg-2
     &   *trc)*gDs**2+ms2*(trsgc+2*ms2*trg+2*mc2*trg+15*ms2*trc)*gDs+ms2
     &   **2*(mc2*(trs+trg)-2*trsgc))+cDs**3*gDs*(4*trs*gDs+ms2*(3*trc-2
     &   *(trs+5*trg)))-mc2*ms2*cDs*gDs*(4*trs*gDs+ms2*(trc-8*trg))+3*ms
     &   2*trc*cDs**4)-2*cDg**3*gDs*(cDs**2*(6*(2*trg-3*trs)*gDs**3+2*(t
     &   rsgc-4*mc2*trs+ms2*trg)*gDs**2+mc2*ms2*(5*trs+6*trg-4*trc)*gDs+
     &   mc2*ms2**2*(trs+trg))-mc2*ms2*(6*(2*trg-3*trs)*gDs**3+2*(trsgc-
     &   4*mc2*trs+ms2*trg)*gDs**2+mc2*ms2*(trs+2*trg)*gDs+mc2*ms2**2*(t
     &   rs+trg))+cDs**3*(6*(2*trg-3*trs)*gDs**2+(-2*trsgc+6*ms2*trg+8*m
     &   s2*trc)*gDs-ms2*(4*trsgc+mc2*trs+2*mc2*trg))+mc2*ms2*cDs*(-2*(-
     &   9*trs+6*trg+4*trc)*gDs**2+(6*trsgc+4*mc2*trs-6*ms2*(trg+2*trc))
     &   *gDs+ms2*(4*trsgc+mc2*(trs+2*trg))))+cDg**2*gDs**2*(cDs**2*(52*
     &   (trg+trc)*gDs**3-8*(trsgc+mc2*trg)*gDs**2+mc2*(-8*trsgc+8*mc2*t
     &   rs+ms2*(-2*trs+trg+10*trc))*gDs-mc2*ms2*(4*trsgc+ms2*(3*trc-2*(
     &   trs+trg))+2*mc2*(trs+trg)))+mc2*ms2*(-52*(trg+trc)*gDs**3+8*(tr
     &   sgc+mc2*(trg-trc))*gDs**2-mc2*(4*(mc2*trs-3*trsgc)+ms2*(-2*trs+
     &   trg+14*trc))*gDs+mc2*ms2*(4*trsgc+ms2*(3*trc-2*(trs+trg))+2*mc2
     &   *(trs+trg)))+cDs**3*(52*(trg+trc)*gDs**2-8*(trsgc+mc2*trg)*gDs-
     &   2*mc2*ms2*(trs+trg+3*trc))+2*mc2*ms2*cDs*(-26*(trg+trc)*gDs**2+
     &   (4*trsgc+2*mc2*(trs+3*trg-trc))*gDs+mc2*ms2*(trs+trg+3*trc)))-2
     &   *mc2*(cDs-mc*ms)*(cDs+mc*ms)*gDs**3*(trs*(cDs+ms2)*(4*gDs**2+ms
     &   2*(4*(gDs+cDs)-mc2)+8*cDs*gDs-2*mc2*gDs+4*cDs**2-2*mc2*cDs+ms2*
     &   *2+mc2**2)+ms2*(trg+trc)*gDs*(2*gDs+2*cDs+ms2-mc2))-2*mc2*cDg*(
     &   cDs-mc*ms)*(cDs+mc*ms)*gDs**3*(-2*trg*gDs**2+cDs*(-2*trg*gDs+ms
     &   2*(4*trs+3*(trg+trc))+2*mc2*trs)+ms2*(4*trs+trg+trc)*gDs+mc2*tr
     &   sgc+2*ms2**2*trs+mc2*ms2*trg)+2*ms2*cDg**5*(mc2*ms2-cDs**2)*(2*
     &   (3*trs*gDs+2*(trg+trc)*gDs+trsgc+ms2*trg)+3*trc*cDs)+8*ms2*trg*
     &   cDg**6*(mc2*ms2-cDs**2))*lVc/(cDg**3*(mc2*ms*ms2-ms*cDs**2)*gDs
     &   **2)/8.0_dp+Bpmm

      Bpmm = xs*cDs*(2*trc*cDg*gDs-mc2*trg*gDs-trsgc*cDg)*(mc2*gDs**2-2
     &   *cDg*cDs*gDs+ms2*cDg**2)*lRcs/(mc*(xs**2-1)*cDg**2*gDs**2)+Bpm
     &   m

      Bpmm = (mc2*cDg**4*(-288*trs*gDs**6-2*cDs*(320*trs*gDs**5-4*ms2*(
     &   -68*trs+88*trg+199*trc)*gDs**4-2*ms2*(-180*trsgc+mc2*(2*trs-132
     &   *trg+51*trc)-34*ms2*trs+24*ms2*trg+298*ms2*trc)*gDs**3+ms2**2*(
     &   252*trsgc+6*ms2*trs+mc2*(-66*trs+104*trg-123*trc)+64*ms2*trg-90
     &   *ms2*trc)*gDs**2-2*ms2**2*(mc2*(9*trsgc+2*ms2*(8*(trs+trg)+9*tr
     &   c))-18*ms2*trsgc)*gDs-9*mc2*ms2**3*(trsgc+ms2*trg))+16*ms2*(-17
     &   5*trs+132*trg+6*trc)*gDs**5+2*ms2*(144*trsgc-652*ms2*trs-4*mc2*
     &   trs+756*ms2*trg+2*mc2*trg+186*ms2*trc+153*mc2*trc)*gDs**4-ms2*(
     &   mc2*(216*trsgc+104*ms2*trs+854*ms2*trg+785*ms2*trc)-12*ms2**2*(
     &   3*trs+21*trg+14*trc)+144*mc2**2*trg)*gDs**3-6*cDs**2*gDs*(2*gDs
     &   +ms2)*(24*trs*gDs**2+ms2*(24*(trs+trg)+13*trc)*gDs+ms2*(-3*trsg
     &   c+6*ms2*trg+3*mc2*(trc-2*trg)-22*ms2*trc))-ms2**2*(mc2*(144*trs
     &   gc+ms2*(50*trs+500*trg+937*trc))-12*ms2*(ms2*(trs+trg)-6*trsgc)
     &   +4*mc2**2*(5*trs+23*trg))*gDs**2+36*ms2*(trs+trg)*cDs**3*gDs*(2
     &   *gDs+ms2)-2*mc2*ms2**3*(-9*trsgc+3*ms2*(trs+7*trg+39*trc)+5*mc2
     &   *trs-13*mc2*trg)*gDs+18*mc2*ms2**4*(trsgc+mc2*trg))-2*mc2*cDg**
     &   3*gDs*(216*(2*trg-3*trs)*gDs**6+cDs*(432*(2*trg-3*trs)*gDs**5+8
     &   *(ms2*(-81*trs+104*trg+50*trc)-28*mc2*trs)*gDs**4-2*ms2*(72*trs
     &   gc+mc2*(10*trs+274*trg+229*trc)+6*ms2*trs-94*ms2*trg-85*ms2*trc
     &   )*gDs**3-ms2*(6*ms2*(12*trsgc+ms2*(trs+trg+3*trc))+mc2*(ms2*(-5
     &   2*trs+210*trg+319*trc)-180*trsgc)+mc2**2*(26*trs-92*trg))*gDs**
     &   2+mc2*ms2**2*(126*trsgc+ms2*(3*trs+32*trg-45*trc)-13*mc2*trs+82
     &   *mc2*trg)*gDs+18*mc2*ms2**3*(trsgc+mc2*trg))+4*(18*(trsgc-3*mc2
     &   *trs)+ms2*(-81*trs+536*trg+464*trc))*gDs**5+2*ms2*(-126*trsgc+m
     &   c2*(656*trs-611*trg+217*trc)-6*ms2*trs+485*ms2*trg+461*ms2*trc)
     &   *gDs**4+cDs**2*(2*gDs+ms2)*(108*(2*trg-3*trs)*gDs**3-3*(12*(trs
     &   gc-mc2*trs)+ms2*(12*trg+trc))*gDs**2+ms2*(19*ms2*trc-12*mc2*(tr
     &   c-3*(trs+trg)))*gDs+9*mc2*ms2*(trsgc+ms2*trg))+ms2*(-3*ms2*(48*
     &   trsgc+ms2*(-2*trs-5*trg+trc))+mc2*(ms2*(680*trs-1001*trg+154*tr
     &   c)-396*trsgc)+6*mc2**2*(29*trs-16*trg+12*trc))*gDs**3+ms2*(-3*m
     &   c2*ms2*(54*trsgc+7*ms2*trs+67*ms2*trg+29*ms2*trc)+mc2**2*(18*tr
     &   sgc+109*ms2*trs+100*ms2*trg+72*ms2*trc)+6*ms2**3*(trs+trg))*gDs
     &   **2+8*ms2*trc*cDs**3*gDs*(2*gDs+ms2)+mc2*ms2**2*(mc2*(27*trsgc+
     &   ms2*(11*trs+92*trg+18*trc))-3*ms2*(ms2*(trs+trg+9*trc)-6*trsgc)
     &   )*gDs+9*mc2**2*ms2**3*(trsgc+ms2*trg))+mc2*cDg**2*gDs**2*((2*gD
     &   s+ms2)*(936*(trg+trc)*gDs**5-36*(4*trsgc-9*mc2*trs+10*mc2*trg)*
     &   gDs**4+4*(-mc2*(45*trsgc+12*ms2*trs+506*ms2*trg+425*ms2*trc)+3*
     &   ms2**2*(4*trs+trg+trc)+72*mc2**2*trs)*gDs**3+ms2*(3*mc2*(72*trs
     &   gc+ms2*(-6*trs-trg+trc))+24*ms2**2*trs+mc2**2*(-462*trs+419*trg
     &   -229*trc))*gDs**2+mc2*ms2*(3*mc2*(78*trsgc+4*ms2*trs+70*ms2*trg
     &   -3*ms2*trc)-6*ms2**2*(trs+trg)+mc2**2*(46*trg-62*trs))*gDs+27*m
     &   c2**2*ms2**3*trc)+2*cDs*gDs*(1872*(trg+trc)*gDs**4+72*(-4*trsgc
     &   +9*mc2*trs+13*ms2*(trg+trc)-10*mc2*trg)*gDs**3+4*(-mc2*(36*trsg
     &   c+ms2*(-75*trs+245*trg+119*trc))-36*ms2*trsgc+3*ms2**2*(4*trs+3
     &   *(trg+trc))+108*mc2**2*trs)*gDs**2+ms2*(mc2*(72*trsgc-ms2*(6*tr
     &   s+304*trg+223*trc))+6*mc2**2*(33*trs+43*trg+5*trc)+6*ms2**2*(4*
     &   trs+3*(trg+trc)))*gDs+3*mc2*ms2**2*(24*trsgc+ms2*(trs+trg+3*trc
     &   )+mc2*(-3*trs+43*trg+5*trc)))+cDs**2*(2*gDs+ms2)*(936*(trg+trc)
     &   *gDs**3-36*(4*trsgc-9*mc2*trs+10*mc2*trg)*gDs**2+3*mc2*(12*trsg
     &   c+ms2*(12*trg+trc))*gDs-19*mc2*ms2**2*trc)-8*mc2*ms2*trc*cDs**3
     &   *(2*gDs+ms2))+2*ms2*cDg**5*(mc2*(6*(32*trs+101*trc)*gDs**4+(-14
     &   4*trsgc+8*mc2*(9*trs+19*trc)+84*ms2*trs-264*ms2*trg-491*ms2*trc
     &   )*gDs**3+(2*mc2*(45*trsgc+62*ms2*trs+107*ms2*trg+46*ms2*trc)-ms
     &   2*(36*trsgc+6*ms2*trs+132*ms2*trg+865*ms2*trc))*gDs**2+ms2*(mc2
     &   *(99*trsgc+ms2*(28*trs+181*trg-12*trc))-6*ms2*(ms2*(trs+trg+39*
     &   trc)-9*trsgc))*gDs+ms2**2*(mc2*(27*trsgc-8*ms2*trs+37*ms2*trg-1
     &   0*ms2*trc)+18*ms2*trsgc))+cDs*(-4*(4*mc2*(6*trs+10*trg+3*trc)-9
     &   *ms2*trc)*gDs**3+2*mc2*(54*trsgc+40*ms2*trs+60*ms2*trg-9*mc2*(t
     &   rc-4*trg)+207*ms2*trc)*gDs**2+mc2*ms2*(90*trsgc+64*ms2*trs+136*
     &   ms2*trg-9*mc2*(trc-4*trg)+219*ms2*trc)*gDs+18*mc2*ms2**2*(trsgc
     &   +ms2*trg))-6*cDs**2*(2*gDs+ms2)*(21*trc*gDs**2+(-9*trsgc-6*mc2*
     &   trg-22*ms2*trc+9*mc2*trc)*gDs+3*mc2*(trsgc+ms2*trg))+36*(trs+tr
     &   g)*cDs**3*gDs*(2*gDs+ms2))+2*mc2*cDg*gDs**3*(2*gDs+ms2)*(18*mc2
     &   *(15*trg+13*trc)*gDs**4+cDs*(12*(2*ms2*trs+mc2*(-2*trs+45*trg+3
     &   9*trc))*gDs**3-12*(mc2*(6*trsgc+ms2*(7*trs+trg+trc))-ms2**2*(6*
     &   trs+trg+trc)+mc2**2*(6*trg-trs))*gDs**2+(-3*mc2**2*(18*trsgc-6*
     &   ms2*trs+41*ms2*trg+23*ms2*trc)+30*ms2**3*trs+70*mc2**3*trs+9*mc
     &   2*ms2**2*(-4*trs+trg+trc))*gDs+18*mc2**2*ms2*(trsgc+mc2*trg))-1
     &   2*(mc2*(3*trsgc+ms2*(2*trs+trg+trc))-ms2**2*(2*trs+trg+trc)+3*m
     &   c2**2*trg)*gDs**3+3*(-3*mc2**2*(6*trsgc+29*ms2*trg+15*ms2*trc)+
     &   2*ms2**3*(4*trs+trg+trc)+12*mc2**3*trs+mc2*ms2**2*(-8*trs-3*(tr
     &   g+trc)))*gDs**2+6*cDs**2*gDs*((8*ms2*trs+mc2*(-8*trs+45*trg+39*
     &   trc))*gDs-2*mc2*(3*trsgc+5*ms2*trs)+mc2**2*(2*trs-6*trg)+8*ms2*
     &   *2*trs)-24*(mc-ms)*(ms+mc)*trs*cDs**3*gDs+ms2*(6*ms2**3-6*mc2*m
     &   s2**2+6*mc2**2*ms2-43*mc2**3)*trs*gDs+18*mc2**3*ms2*(trsgc+ms2*
     &   trg))+2*mc2**2*gDs**4*(2*gDs+ms2)*(18*mc2*trg*gDs**3+12*(ms2*tr
     &   s-mc2*trs+3*mc2*trg)*cDs*gDs**2+6*ms2*(ms2-mc2)*(2*trs+trg+trc)
     &   *gDs**2+6*cDs**2*(3*mc2*trg*gDs+(mc-ms)*(ms+mc)*trs*(-4*gDs-4*m
     &   s2+mc2))-9*mc2**2*trsgc*(gDs+cDs+ms2)+6*(mc-ms)*(ms+mc)*(mc2*tr
     &   s-ms2*(6*trs+trg+trc))*cDs*gDs+mc2**2*ms2*(6*trs+2*trg+29*trc)*
     &   gDs+3*ms2**3*(4*trs+trg+trc)*gDs-6*mc2*ms2**2*(3*trs+trg+trc)*g
     &   Ds-12*(mc-ms)*(ms+mc)*trs*cDs**3+15*ms2**3*trs*cDs-24*mc2*ms2**
     &   2*trs*cDs+5*mc2**3*trs*cDs-3*mc2**2*ms2*(3*trg-4*trs)*cDs-mc2**
     &   3*ms2*(11*trs+9*trg)+3*ms2**4*trs-6*mc2*ms2**3*trs+6*mc2**2*ms2
     &   **2*trs)-2*ms2*cDg**6*(2*gDs+ms2)*(-mc2*((72*(trs+2*trg)+89*trc
     &   )*gDs**2+(117*trsgc+ms2*(108*trs+234*trg+97*trc)+90*mc2*trg)*gD
     &   s+2*ms2*(27*trsgc-8*ms2*trs+19*ms2*trg+18*mc2*trg-10*ms2*trc))+
     &   6*cDs*(42*trc*gDs**2+(-18*trsgc+3*mc2*(3*trs-trg+5*trc)-49*ms2*
     &   trc)*gDs+6*mc2*(trsgc+ms2*trg))+36*cDs**2*(2*trc*gDs+trsgc+(ms2
     &   +mc2)*trg))-36*ms2*cDg**7*(2*gDs+ms2)*(7*trc*gDs**2+cDs*((6*(tr
     &   s+trg)+8*trc)*gDs+4*(trsgc+(ms2+mc2)*trg))+(2*mc2*(trs-4*trg+tr
     &   c)-3*(trsgc+3*ms2*trc))*gDs+4*trg*cDs**2+mc2*(trsgc-3*ms2*trg))
     &   -72*ms2*cDg**8*(2*gDs+ms2)*(2*(trs+trg+trc)*gDs+trg*(4*cDs+ms2+
     &   mc2)+trsgc)-144*ms2*trg*cDg**9*(2*gDs+ms2))/(mc2*ms*cDg**3*(2*c
     &   Dg+mc2)*gDs**3*(2*gDs+ms2))/2.4e+1_dp+Bpmm

      Bpmm = B0cgsf*(cDg**4*(-2*mc2*(4*trs*gDs**4+ms2*(3*trs-5*trg-11*t
     &   rc)*gDs**3+ms2*(2*mc2*(trs+trg)+ms2*(trs+5*trc))*gDs**2+ms2**2*
     &   (2*trsgc+ms2*(trs+13*trc))*gDs-ms2**3*(trsgc+mc2*trg))+cDs*(-2*
     &   (4*mc2*trs+7*ms2*trc)*gDs**3+ms2*(6*(trsgc+3*ms2*trc)+mc2*(-16*
     &   trs+16*trg+7*trc))*gDs**2+4*mc2*ms2**2*(2*trs+3*trg+2*trc)*gDs+
     &   2*mc2*ms2**2*(trsgc+ms2*trg))+2*ms2*cDs**2*gDs*((2*(trs+trg)-7*
     &   trc)*gDs+3*trsgc+2*mc2*trg+9*ms2*trc-3*mc2*trc)+4*ms2*(trs+trg)
     &   *cDs**3*gDs)+mc2*cDg**2*gDs**2*(4*(9*trs+7*trg+13*trc)*gDs**4-(
     &   12*trsgc-16*mc2*trs+55*ms2*trg+8*mc2*trg+51*ms2*trc)*gDs**3+cDs
     &   *gDs*(4*(18*trs+trg+13*trc)*gDs**2-(8*(trsgc+mc2*(trg-3*trs))+m
     &   s2*(23*trg+22*trc))*gDs+ms2*(8*trsgc+ms2*(2*(trs+trg)+3*trc)-4*
     &   mc2*trs+12*mc2*trg))+(-mc2*(8*trsgc+ms2*(40*trs-24*trg+7*trc))+
     &   ms2*(8*trsgc+ms2*(2*trs-trg+trc))+8*mc2**2*trs)*gDs**2+cDs**2*(
     &   4*gDs*(9*trs*gDs-6*trg*gDs+trsgc)+4*ms2*trg*gDs-2*ms2**2*trc)-2
     &   *ms2*(2*mc2*(ms2*(-trs-3*trg+trc)-2*trsgc)+ms2**2*(trs+trg)+4*m
     &   c2**2*trs)*gDs-ms2*trc*cDs**3+3*mc2*ms2**3*trc)+mc2*cDg**3*gDs*
     &   (-4*(6*trg-7*trs)*gDs**4-2*cDs*(2*(6*trg-5*trs)*gDs**3-(2*trsgc
     &   -5*ms2*trs+11*ms2*trg+6*ms2*trc)*gDs**2-ms2*(trsgc-ms2*trs+2*mc
     &   2*trs-2*ms2*trg+5*ms2*trc)*gDs+2*ms2**2*(trsgc+mc2*trg))+(-4*tr
     &   sgc+16*mc2*trs+ms2*(-38*trs+31*trg+8*trc))*gDs**3-cDs**2*(8*trs
     &   *gDs**2+ms2*(8*trs-8*trg-3*trc)*gDs+2*ms2*(trsgc+ms2*trg))+ms2*
     &   (8*trsgc+2*ms2*trs-4*mc2*(trg+2*trc)+10*ms2*trg+9*ms2*trc)*gDs*
     &   *2+2*ms2*(mc2*(trsgc-ms2*(2*trs+9*trg+2*trc))+ms2**2*(trs+trg))
     &   *gDs-2*mc2*ms2**2*(trsgc+ms2*trg))+2*mc2*cDg*gDs**3*(26*(trg+tr
     &   c)*gDs**4+cDs*(52*(trg+trc)*gDs**3-2*(4*trsgc+3*mc2*trg)*gDs**2
     &   +(-2*mc2*(2*trsgc+ms2*(trs+7*trg+5*trc))+ms2**2*(4*trs+3*(trg+t
     &   rc))+6*mc2**2*trs)*gDs+2*mc2*ms2*(trsgc+mc2*trg))-2*(2*trsgc+mc
     &   2*trg)*gDs**3+(-2*mc2*(2*trsgc+ms2*(2*trs+14*trg+11*trc))+ms2**
     &   2*(4*trs+trg+trc)+4*mc2**2*trs)*gDs**2+2*cDs**2*gDs*(13*(trg+tr
     &   c)*gDs-2*(trsgc+mc2*trg))+(2*mc2*ms2*(trsgc-ms2*trs)+mc2**2*(ms
     &   2*(trg-2*trs)-trsgc)+2*ms2**3*trs)*gDs+2*mc2**2*ms2*(trsgc+ms2*
     &   trg))+2*mc2*gDs**4*(2*mc2*trg*gDs**3+4*(ms2*trs+mc2*(trg-trs))*
     &   cDs*gDs**2+2*ms2*(ms2-mc2)*(2*trs+trg+trc)*gDs**2+2*cDs**2*(mc2
     &   *trg*gDs+(mc-ms)*(ms+mc)*trs*(-4*gDs-4*ms2+mc2))-mc2**2*trsgc*(
     &   gDs+cDs+ms2)+2*(mc-ms)*(ms+mc)*(mc2*trs-ms2*(6*trs+trg+trc))*cD
     &   s*gDs+ms2**3*(4*trs+trg+trc)*gDs-2*mc2*ms2**2*(3*trs+trg+trc)*g
     &   Ds+mc2**2*ms2*(2*trs+trg+4*trc)*gDs+4*(ms2-mc2)*trs*cDs**3+5*ms
     &   2**3*trs*cDs-8*mc2*ms2**2*trs*cDs-mc2**2*ms2*(trg-4*trs)*cDs-mc
     &   2**3*ms2*(2*trs+trg)+ms2**4*trs-2*mc2*ms2**3*trs+2*mc2**2*ms2**
     &   2*trs)+ms2*cDg**5*(-14*trc*gDs**3+(6*(trsgc+3*ms2*trc)+mc2*(-4*
     &   trs+12*trg+9*trc))*gDs**2-2*cDs*gDs*(2*(trs+trg+9*trc)*gDs-4*tr
     &   sgc+2*ms2*trg-4*mc2*trg-18*ms2*trc+3*mc2*trc)-4*cDs**2*(2*trc*g
     &   Ds+trsgc+ms2*trg)+2*mc2*(3*trsgc+ms2*(6*trs+11*trg+trc))*gDs-2*
     &   mc2*ms2*(ms2*(trs-2*trg+trc)-3*trsgc))-2*ms2*cDg**6*((4*(trs+tr
     &   g)+11*trc)*gDs**2+2*cDs*((3*trs+5*trg+4*trc)*gDs+2*(trsgc+ms2*t
     &   rg))-(trsgc-2*ms2*trg+6*mc2*trg+9*ms2*trc)*gDs+4*trg*cDs**2-4*m
     &   c2*ms2*trg)-4*ms2*cDg**7*(2*(trs+2*trg+trc)*gDs+4*trg*cDs+trsgc
     &   +ms2*trg)-8*ms2*trg*cDg**8)/(mc2*ms*cDg**3*gDs**3)/8.0_dp+Bpmm

      Bpmm = tr3s00ft*(cDg**2*(12*(2*trg-3*trs)*gDs**4+cDs*(24*(2*trg-3
     &   *trs)*gDs**3+37*ms2*(trg+trc)*gDs**2-2*ms2*(4*trsgc-9*mc2*trs+8
     &   *mc2*trg)*gDs+2*mc2*ms2*(trsgc+ms2*trg))+(4*trsgc-16*mc2*trs+10
     &   5*ms2*trg+101*ms2*trc)*gDs**3+2*ms2*(mc2*(27*trs-19*trg+trc)-8*
     &   trsgc)*gDs**2+4*cDs**2*gDs*(-9*trs*gDs+6*trg*gDs-trsgc-ms2*trg)
     &   +2*mc2*ms2*(-7*trsgc+4*mc2*trs-7*ms2*trg-2*ms2*trc)*gDs+2*mc2*m
     &   s2**2*(trsgc+mc2*trg))-2*cDg*gDs*(26*(trg+trc)*gDs**4+cDs*(52*(
     &   trg+trc)*gDs**3-8*(trsgc+mc2*trg)*gDs**2+mc2*(-4*trsgc+2*mc2*tr
     &   s-17*ms2*trg-13*ms2*trc)*gDs+2*mc2*ms2*(trsgc+mc2*trg))-4*(trsg
     &   c+mc2*trg)*gDs**3+mc2*(-4*trsgc+4*mc2*trs-27*ms2*trg-19*ms2*trc
     &   )*gDs**2+2*cDs**2*gDs*(13*(trg+trc)*gDs-2*(trsgc+mc2*trg))+mc2*
     &   ms2*(2*trsgc-9*mc2*trs+mc2*trg)*gDs+2*mc2**2*ms2*(trsgc+ms2*trg
     &   ))+cDg**3*(8*trs*gDs**3+2*cDs*(16*trs*gDs**2+ms2*(15*trs-4*trg+
     &   4*trc)*gDs+ms2*(4*trsgc-2*mc2*trs+5*ms2*trg+3*ms2*trc))+ms2*(78
     &   *trs-59*trg-3*trc)*gDs**2+8*trs*cDs**2*gDs-2*ms2*(4*trsgc-7*mc2
     &   *trs+6*ms2*trg+4*ms2*trc)*gDs+2*ms2**2*(2*trsgc+5*mc2*trg))+2*m
     &   c2*gDs**2*(-2*trg*gDs**3-4*trg*cDs*gDs**2+(mc2*(trsgc-2*ms2*trg
     &   -5*ms2*trc)-2*trg*cDs**2)*gDs+mc2*((trsgc+ms2*trg)*cDs+ms2*(trs
     &   gc+3*mc2*trs+mc2*trg)))+2*ms2*cDg**4*(-5*trs*gDs-4*trs*cDs+7*ms
     &   2*trg+ms2*trc))/(ms*cDg**3*gDs)/2.0_dp+Bpmm

      Bpmm = lc*(cDg**2*(52*(trg+trc)*gDs**4+2*cDs*(26*(trg+trc)*gDs**3
     &   +(-4*trsgc+9*mc2*trs-10*mc2*trg)*gDs**2+mc2*(trsgc-7*ms2*trg)*g
     &   Ds+2*mc2**2*ms2*trg)-2*(4*trsgc-9*mc2*trs+10*mc2*trg)*gDs**3+2*
     &   mc2*(-5*trsgc+8*mc2*trs-ms2*trg)*gDs**2-6*mc2**2*ms2*trg*gDs+mc
     &   2**2*ms2*(trsgc+ms2*trg))-2*cDg**3*(6*(2*trg-3*trs)*gDs**3+2*(t
     &   rg*(6*cDs+ms2)-3*trs*(3*cDs+mc2)+trsgc)*gDs**2-(mc2*(ms2*(trg+t
     &   rc)-2*trs*cDs)+2*(trsgc-3*ms2*trg)*cDs)*gDs+mc2*ms2*(-2*trg*(4*
     &   cDs+ms2)-trsgc+mc2*trs))+2*mc2*cDg*gDs*((15*trg+13*trc)*gDs**3+
     &   cDs*((15*trg+13*trc)*gDs**2-2*(trsgc+mc2*trg)*gDs-2*mc2*ms2*trg
     &   )-2*(trsgc+mc2*trg)*gDs**2+mc2*(-3*trsgc+2*mc2*trs+ms2*trg)*gDs
     &   -mc2**2*ms2*trg)+mc2**2*gDs**2*(2*trg*gDs**2+2*trg*cDs*gDs-mc2*
     &   trsgc+mc2*ms2*trg)+2*cDg**4*(ms2*((5*trg+trc)*gDs+2*trg*(4*cDs+
     &   ms2))-trs*(4*gDs*(gDs+cDs)+3*mc2*ms2))-4*ms2*trs*cDg**5)/(ms*cD
     &   g*(2*cDg+mc2)**2*gDs)/2.0_dp+Bpmm

      Bpmm = LsB2*ms*(2*mc2**3*trg*cDs*gDs**4+mc2*cDg**2*gDs**2*(mc2*(2
     &   *(mc2*trs+ms2*trc)*gDs-ms2*(trsgc+3*mc2*trg))+cDs**2*(4*trc*gDs
     &   -2*trsgc+6*mc2*trg)+12*trg*cDs**3+3*mc2*(trsgc-ms2*trg)*cDs)+cD
     &   g**3*gDs*(3*mc2*cDs*(ms2*(trsgc+3*mc2*trg)-2*(mc2*trs+ms2*trc)*
     &   gDs)-8*trg*cDs**4-12*mc2*trg*cDs**3+2*mc2*(ms2*trg-3*trsgc)*cDs
     &   **2+mc2**2*ms2*(3*trsgc+ms2*trg))-mc2**2*cDg*gDs**3*(cDs*(2*trc
     &   *gDs-trsgc+mc2*trg)+8*trg*cDs**2+mc2*(trsgc-ms2*trg))+cDg**4*(2
     &   *mc2*(6*trs*cDs**2-3*mc2*ms2*trs+ms2**2*trc)*gDs+8*trg*cDs**4+4
     &   *(trsgc+ms2*trg)*cDs**3-8*mc2*ms2*trg*cDs**2-3*mc2*ms2*(trsgc+m
     &   s2*trg)*cDs+mc2*ms2**2*(mc2*trg-trsgc))+cDg**5*(6*mc2*ms2*trs*c
     &   Ds-8*trs*cDs**3))/(cDg**3*(cDs**2-mc2*ms2)*gDs)+Bpmm

      Bpmm = Bpmm-ls*ms*(-2*cDg**2*(mc2*(-8*trs*gDs**3+2*(trsgc-2*ms2*
     &   trs+ms2*trg)*gDs**2+ms2*(3*trsgc-2*ms2*trs-3*ms2*trg)*gDs+ms2**
     &   2*(trsgc-ms2*trg))+cDs*gDs*(-14*trc*gDs**2+(6*trsgc+4*mc2*trg+5
     &   *ms2*trc)*gDs+ms2*(3*trsgc+2*mc2*trg+9*ms2*trc))+2*(trs+trg)*cD
     &   s**2*gDs*(2*gDs+ms2))-mc2*cDg*gDs*(10*trc*gDs**3+(8*(-2*trs*cDs
     &   +4*trg*cDs+trsgc)+ms2*trc)*gDs**2+(4*ms2*(-2*trs*cDs+11*trg*cDs
     &   +3*trsgc)-4*trsgc*cDs-3*ms2**2*trc)*gDs+2*ms2*(-trsgc*cDs+7*ms2
     &   *trg*cDs+2*ms2*trsgc))+mc2*gDs**2*(trc*cDs*(-2*gDs**2+3*ms2*gDs
     &   +3*ms2**2)+2*mc2*(2*gDs+ms2)*(-2*trs*gDs+4*trg*gDs-trsgc+3*ms2*
     &   trg))+2*cDg**3*(2*gDs+ms2)*(2*cDs*((trs+trg+2*trc)*gDs+trsgc+ms
     &   2*trg)+gDs*(7*trc*gDs-3*trsgc-2*mc2*trg-9*ms2*trc))+4*cDg**4*(2
     &   *gDs+ms2)*(2*(trs+trg+trc)*gDs+2*trg*cDs+trsgc+ms2*trg)+8*trg*c
     &   Dg**5*(2*gDs+ms2))/(mc2*cDg*gDs*(2*gDs+ms2)**2)/4.0_dp

      Bpmm = ms*tr3c00fs*(-2*cDg**2*(mc2*(17*trc*gDs**3-2*(2*trsgc+2*mc
     &   2*trg+9*ms2*trc)*gDs**2-13*ms2**2*trc*gDs+ms2**2*(trsgc+mc2*trg
     &   ))+mc2*cDs*(4*(trg+2*trc)*gDs**2+2*ms2*(6*trs+7*trg+2*trc)*gDs+
     &   ms2*(trsgc+ms2*trg))+cDs**2*gDs*(-7*trc*gDs+3*trsgc+2*mc2*trg)+
     &   2*(trs+trg)*cDs**3*gDs)+2*cDg**3*(-mc2*(4*(trs+2*trg+2*trc)*gDs
     &   **2+(5*trsgc+6*ms2*trs+13*ms2*trg)*gDs+ms2*(3*trsgc-3*ms2*trs-2
     &   *ms2*trc))+cDs*gDs*(14*trc*gDs-6*trsgc-4*mc2*trg-9*ms2*trc)+2*c
     &   Ds**2*(2*trc*gDs+trsgc+ms2*trg))+2*mc2*cDg*gDs*(mc2*(4*(trg+trc
     &   )*gDs**2+(trsgc+3*ms2*trs+6*ms2*trg+2*ms2*trc)*gDs+ms2*(trsgc+m
     &   s2*trg))+cDs*(-17*trc*gDs**2+4*(trsgc+mc2*trg-2*ms2*trc)*gDs+2*
     &   ms2*(trsgc+mc2*trg))+cDs**2*(4*(trs+trg)*gDs+trsgc+ms2*trg))+mc
     &   2**2*gDs**2*(20*trc*gDs**2+(5*ms2*trc-2*(2*(trs+trg)*cDs+trsgc+
     &   2*mc2*trg))*gDs-2*ms2*(trg*(cDs+mc2)+trsgc)-2*trsgc*cDs-3*ms2**
     &   2*trc)+2*cDg**4*(7*trc*gDs**2+cDs*((6*(trs+trg)+8*trc)*gDs+4*(t
     &   rsgc+ms2*trg))-(3*trsgc+10*mc2*trg+9*ms2*trc)*gDs+4*trg*cDs**2-
     &   4*mc2*ms2*trg)+4*cDg**5*(2*(trs+trg+trc)*gDs+4*trg*cDs+trsgc+ms
     &   2*trg)+8*trg*cDg**6)/(mc2*cDg*gDs**3)/2.0_dp+Bpmm

      Bpmm = BfunX*(-ms2*cDg**4*gDs*(mc2*((2*(trs+trg)+3*trc)*gDs**2+ms
     &   2*(2*(trs+trg)+trc)*gDs-2*ms2**2*(trs+trg))+cDs*(12*trc*gDs**2+
     &   mc2*(4*trs+9*trc)*gDs+4*mc2*ms2*(trs+trg))+6*trc*cDs**2*(4*gDs+
     &   mc2)+12*trc*cDs**3)-mc2*ms2*cDg**2*gDs**2*(gDs*((20*trs+3*trg+7
     &   *trc)*gDs**2+ms2*(6*trs-5*trg+trc)*gDs-2*ms2**2*(trs+trg))+cDs*
     &   *2*((4*(trs+trg)+11*trc)*gDs+ms2*trc)+cDs*gDs*(24*trs*gDs+7*trg
     &   *gDs+17*trc*gDs-2*ms2*trs-2*ms2*trg+6*ms2*trc)+trc*cDs**3)-ms2*
     &   cDg**5*(-mc2*(trc*gDs**2+2*ms2*(2*(trs+trg)+trc)*gDs+2*ms2**2*(
     &   trs+trg+trc))+6*trc*cDs*gDs*(4*gDs+mc2)+24*trc*cDs**2*gDs)-mc2*
     &   ms2*cDg**3*gDs**2*((6*trs+trg+7*trc)*gDs**2+2*cDs*((trs-trg+5*t
     &   rc)*gDs+ms2*(trs-trg))+2*ms2*(trs+trg+2*trc)*gDs-(4*trs+8*trg+3
     &   *trc)*cDs**2+2*ms2**2*(trs+trg))-2*mc2*cDg*gDs**4*(2*(gDs+cDs)+
     &   ms2)*(2*trs*(cDs+2*ms2)*(2*(gDs+cDs)+ms2)+3*ms2*(trg+trc)*(gDs+
     &   cDs))-2*mc2*gDs**4*(2*(gDs+cDs)+ms2)**2*(trs*(cDs+ms2)*(2*(gDs+
     &   cDs)+ms2)+ms2*(trg+trc)*gDs)-12*ms2*trc*cDg**6*cDs*gDs)/(mc2*ms
     &   *cDg**3*gDs**3)/8.0_dp+Bpmm

      Bpmm = Bpmm-LsB1*ms*(cDg**2*gDs**2*(mc2*cDs*(12*trc*gDs**2+2*(-3
     &   *trsgc+3*ms2*trs-ms2*trg+5*ms2*trc)*gDs+ms2*(trsgc-5*ms2*trg))+
     &   4*cDs**2*(-2*trc*gDs**2+(trsgc+mc2*trg-2*ms2*trc)*gDs+ms2*(trsg
     &   c+2*mc2*trg))+mc2*ms2*(-4*trc*gDs**2+2*(trsgc+mc2*trg-2*ms2*trc
     &   )*gDs+ms2*(mc2*trg-trsgc))+cDs**3*(8*(trg+trc)*gDs-4*trsgc+8*ms
     &   2*trg))-cDg**3*gDs*(2*cDs**2*(4*trc*gDs**2-2*trsgc*gDs+2*ms2*(t
     &   rs+trg+4*trc)*gDs-3*ms2*trsgc+ms2**2*trg)+mc2*ms2*(4*trc*gDs**2
     &   +2*(ms2*(-trs-trg+trc)-trsgc)*gDs+ms2*(trsgc-ms2*trg))+ms2*cDs*
     &   (-4*trc*gDs**2+2*(trsgc+mc2*trg-2*ms2*trc)*gDs+ms2*(trsgc+3*mc2
     &   *trg)))-mc2*cDg*gDs**3*(mc2*(4*trc*gDs**2+2*(ms2*(3*trs-trg+trc
     &   )-trsgc)*gDs+ms2*(trsgc-5*ms2*trg))+cDs*(-12*trc*gDs**2+6*(trsg
     &   c+mc2*trg-2*ms2*trc)*gDs+ms2*(3*trsgc+5*mc2*trg))+4*cDs**2*(2*t
     &   rg*(gDs+ms2)+trc*gDs)+4*trg*cDs**3)+mc2*gDs**4*(-4*mc2*trc*gDs*
     &   *2+2*mc2*((trs+trg)*cDs+trsgc+mc2*trg-2*ms2*trc)*gDs+2*(mc2*trg
     &   -trsgc)*cDs**2+mc2*(trsgc+ms2*trg)*cDs+mc2*ms2*(3*trsgc+mc2*trg
     &   ))+2*ms2*cDg**4*cDs*(trc*gDs*(2*gDs+3*ms2)-trsgc*(gDs+ms2)))/(c
     &   Dg*(cDs-mc*ms)*(cDs+mc*ms)*gDs**3)

      Bpmm = Bpmm-B0csf*ms*(mc2*gDs*(mc2*(2*trc*gDs**2+(-trsgc-mc2*trs
     &   +ms2*(trg-trs))*gDs+ms2*(ms2+mc2)*trg)+cDs*(2*trc*gDs**2+(-trsg
     &   c+mc2*(-2*trs-trg+trc)+ms2*trc)*gDs+2*mc2*ms2*trg)+cDs**2*(2*tr
     &   c*gDs-trg*(2*gDs+ms2+mc2))-2*trg*cDs**3)+cDg*(2*mc2*cDs*(-2*trc
     &   *gDs**2+(trsgc+ms2*(trs-trg)+mc2*trs)*gDs-ms2*(trsgc+(ms2+mc2)*
     &   trg))+cDs**2*(-4*trc*gDs**2+2*(trsgc+mc2*(2*trs+trg-trc)-ms2*tr
     &   c)*gDs+mc2*(trsgc-4*ms2*trg)+ms2*trsgc)+2*cDs**3*(trg*(2*gDs+ms
     &   2+mc2)-2*trc*gDs+trsgc)+4*trg*cDs**4-mc2*ms2*(ms2+mc2)*trsgc)+c
     &   Dg**2*(ms2*cDs*(2*trc*gDs-trsgc+mc2*(2*trs-trg+trc)+ms2*trc)+mc
     &   2*ms2*(2*trc*gDs-trsgc+ms2*(trs+trg)+mc2*trs)-4*trs*cDs**3+2*(m
     &   s2*(-trs-trg+trc)-mc2*trs)*cDs**2))/(cDg*(cDs-mc*ms)*(cDs+mc*ms
     &   )*gDs)/2.0_dp

      Bpmm = tr3s002ft*(6*trs*(2*cDg+mc2)*cDs*((2*cDg+mc2)**2*gDs**2+ms
     &   2*cDg**3)+3*ms2*(trg+trc)*cDg**2*gDs*(mc2*(gDs+6*cDs)+10*cDg*cD
     &   s))/(ms*cDg**3*gDs)/2.0_dp+Bpmm

      Bpmm = 3.0_dp*ms*tr3c002fs*(trc*cDg**3*(gDs**2+2*ms2*gDs+2*ms2**2
     &   )+cDg*cDs*(3*trc*cDs-4*mc2*(trs+trg))*gDs**2-cDg**2*gDs*(3*trc*
     &   cDs*gDs+2*mc2*ms2*(trs+trg))+gDs**2*(2*mc2**2*(trs+trg)*gDs-trc
     &   *cDs**3))/(2.0_dp*cDg*gDs**3)+Bpmm

      Bpmm = 3.0_dp*ms*tr3c001fs*(-mc2*ms2*trc*cDs**2*gDs**2+2*mc2*cDg*
     &   gDs**2*(mc2*(trs+trg)*(2*gDs+ms2)-ms2*trc*cDs)-cDg**2*gDs*(8*mc
     &   2*(trs+trg)*cDs*(gDs+ms2)+mc2*ms2*trc*gDs-6*ms2*trc*cDs**2)+2*m
     &   s2*cDg**3*(3*trc*cDs*gDs+mc2*ms2*(trs+trg)))/(2.0_dp*mc2*cDg*gD
     &   s**3)+Bpmm

      Bpmm = epinv*(4*lp*xs*cDs+mc*(ms-ms*xs**2))*(2*trc*cDg*gDs-mc2*tr
     &   g*gDs-trsgc*cDg)*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2)/(mc*(xs*
     &   *2-1)*cDg**2*gDs**2)/4.0_dp+Bpmm

      Bpmm = ms*tr2fu*cDs*(2*trc*cDg*gDs-mc2*trg*gDs-trsgc*cDg)*(mc2*gD
     &   s**2-2*cDg*cDs*gDs+ms2*cDg**2)/(cDg**2*gDs**2)+Bpmm

      Bpmm = tr3s001ft*(3*ms*(trg+trc)*(2*ms2*cDg**3*cDs-gDs**2*(2*mc2*
     &   cDg*(gDs-3*cDs)+cDg**2*(gDs-7*cDs)+2*mc2**2*gDs))+6*ms*trs*gDs*
     &   (cDg**3*(gDs+5*cDs)+3*mc2*cDg**2*(gDs+cDs)+3*mc2**2*cDg*gDs+mc2
     &   **3*gDs))/(cDg**3*gDs)/2.0_dp+Bpmm

      Bpmm=Bpmm/zpmm

      return
      end
