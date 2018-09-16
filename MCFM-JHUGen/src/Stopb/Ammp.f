      subroutine Aamp_mmp(q,mc,ms,Ammp)
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
      complex(dp):: trc,trg,trs,trsgc,zmmp,Ammp

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

      Ammp = ms2*(cDg**2*(2*trs*gDs**2-ms2*trg*gDs+ms2*trsgc)+mc2*trsgc
     &   *gDs**2-trsgc*cDg*gDs*(gDs+2*cDs))*tr3Xs/gDs**3

      Ammp = mc2*(cDg**2*(2*(trs+trc)*gDs**2-(trsgc+ms2*trg)*gDs+ms2*tr
     &   sgc)+mc2*trsgc*gDs**2-cDg*gDs*((trsgc+mc2*trg)*gDs+2*trsgc*cDs)
     &   )*tr3Xc/(cDg**2*gDs)+Ammp

      Ammp = (cDg*(4*trs*gDs+ms2*trsgc/gDs-2*ms2*trg)-2*trsgc*(gDs+cDs)
     &   +mc2*trsgc*gDs/cDg)*tr1Xfs+Ammp

      Ammp = (cDg**2*(4*trs*gDs**2-2*ms2*trg*gDs+ms2*trsgc)+mc2*trsgc*g
     &   Ds**2-2*trsgc*cDg*gDs*(gDs+cDs))*tr1Xfc/gDs**2+Ammp

      Ammp = (-mc2*cDg**3*(4*(trs+trg+trc)*gDs**2+2*(ms2*(trc-trg)-trsg
     &   c)*gDs+ms2*(trsgc-2*ms2*trs-3*ms2*trg))+2*mc2*cDg*gDs**2*(ms2*t
     &   rs*(2*(gDs+cDs)-mc2)-mc2*trsgc-ms2**2*(-trs+trg+trc))+4*cDg**5*
     &   ((trs+trg)*gDs-ms2*trc)+mc2*cDg**2*gDs*(4*ms2*trs*gDs+2*mc2*trg
     &   *gDs+4*trsgc*cDs+ms2*trsgc+mc2*ms2*trg)+mc2*ms2*(trs+trg+trc)*g
     &   Ds**3*(2*gDs+2*cDs+ms2-mc2))*lVs/(mc2*cDg**2*gDs**2)/4.0_dp+Am
     &   mp

      Ammp = (-2*cDg**3*(4*trs*gDs**2+2*(mc2*(trs+trg)-ms2*trg)*gDs+ms2
     &   *(trsgc+mc2*(trs+trg)))+cDg**2*gDs*(trsgc*(4*(gDs+cDs)+mc2)+mc2
     &   *(ms2*trg-2*(3*trs+trg+2*trc)*gDs))+mc2*cDg*gDs**2*(-4*trs*(gDs
     &   +cDs)-trsgc+mc2*(2*trs+trg)+2*ms2*(-trs+trg+trc))-mc2*(trs+trg+
     &   trc)*gDs**3*(2*gDs+2*cDs+ms2-mc2)+4*trc*cDg**5)*lVc/(cDg**2*gDs
     &   **2)/4.0_dp+Ammp

      Ammp = (cDg**2*(4*trs*gDs**2-2*ms2*trg*gDs+ms2*trsgc)+mc2*trsgc*g
     &   Ds**2-2*trsgc*cDg*gDs*(gDs+cDs))*lRs2/(cDg*gDs**2)/2.0_dp+Ammp

      Ammp = epinv*(cDg**2*(4*trs*gDs**2-2*ms2*trg*gDs+ms2*trsgc)+mc2*t
     &   rsgc*gDs**2-2*trsgc*cDg*gDs*(gDs+cDs))*(2*lRs1+2*lRc1-1)/(cDg*g
     &   Ds**2)/4.0_dp+Ammp

      Ammp = (cDg**2*(4*trs*gDs**2-2*ms2*trg*gDs+ms2*trsgc)+mc2*trsgc*g
     &   Ds**2-2*trsgc*cDg*gDs*(gDs+cDs))*lRc2/(cDg*gDs**2)/2.0_dp+Ammp

      Ammp = 2*LsA*(-mc2**2*trsgc*gDs**4+cDg**3*gDs*(cDs*(4*trs*gDs**2-
     &   2*ms2*trg*gDs+3*ms2*trsgc)+ms2*trsgc*gDs)+cDg**2*gDs**2*(-2*mc2
     &   *(trs+trc)*gDs**2+(-2*trsgc*cDs+mc2*trsgc+mc2*ms2*trg)*gDs-2*tr
     &   sgc*(cDs**2+mc2*ms2))+ms2*cDg**4*(-2*trs*gDs**2+ms2*trg*gDs-ms2
     &   *trsgc)+mc2*cDg*gDs**3*((trsgc+mc2*trg)*gDs+3*trsgc*cDs))/(cDg*
     &   *2*gDs**3)+Ammp

      Ammp = (-mc2*cDg**3*(2*gDs+ms2)*(8*(-6*trsgc+4*mc2*(3*trs-trg+2*t
     &   rc)+3*ms2*trs)*gDs**3+2*mc2*(-9*trsgc+ms2*(21*trs-43*trg+10*trc
     &   )+10*mc2*trs+19*mc2*trg)*gDs**2+9*cDs*(trsgc*(-12*gDs**2+2*(mc-
     &   ms)*(ms+mc)*gDs+mc2*ms2)+mc2*(2*gDs+ms2)*(2*(-trs-trg+trc)*gDs+
     &   ms2*trg))+mc2*ms2*(27*trsgc+10*mc2*trs-2*ms2*(-3*trs+6*trg+8*tr
     &   c)+37*mc2*trg)*gDs+9*mc2*ms2**2*(trsgc+mc2*trg))+mc2*cDg**2*gDs
     &   *(2*gDs+ms2)*(4*(mc2*(7*trs+8*trg+17*trc)-6*ms2*trs)*gDs**3+3*c
     &   Ds*(4*(3*mc2*trc-(2*ms2+mc2)*trs)*gDs**2+6*mc2*(4*trsgc+(ms2+mc
     &   2)*trg)*gDs+3*mc2*ms2*(trsgc+mc2*trg))-2*(mc2*(15*trsgc+ms2*(-6
     &   *trs+31*trg+22*trc))+mc2**2*(2*trs-11*trg+5*trc)-6*ms2**2*(-trs
     &   +trg+trc))*gDs**2+18*mc2*(ms2+mc2)*trsgc*gDs-2*mc2**2*ms2*(9*tr
     &   s-26*trg+9*trc)*gDs+9*mc2**2*ms2*(trsgc+ms2*trg))+2*mc2*cDg**5*
     &   (8*(17*trs+8*trg+5*trc)*gDs**3+(-36*trsgc-4*mc2*(3*trc-8*(trs+t
     &   rg))+2*ms2*(81*trs+27*trg+22*trc))*gDs**2+16*(-trs-trg+trc)*cDs
     &   *gDs*(2*gDs+ms2)+2*ms2*(-18*trsgc+mc2*(8*(trs+trg)-3*trc)+ms2*(
     &   33*trs+6*trg+16*trc))*gDs-9*ms2**2*trsgc+ms2**3*(8*trs-trg+10*t
     &   rc))-mc2*cDg**4*(2*gDs+ms2)*(32*(3*trs+2*trc)*gDs**3-4*(-9*trsg
     &   c+mc2*(7*trs-11*trg+trc)+ms2*(-3*trs+9*trg+8*trc))*gDs**2+18*cD
     &   s*(2*gDs+ms2)*(2*(trc-trs)*gDs+trsgc+ms2*trg)+18*(5*ms2+mc2)*tr
     &   sgc*gDs-2*ms2*(2*ms2*(8*trc-3*(trs+trg))+mc2*(15*(trs-2*trg)+11
     &   *trc))*gDs+9*ms2*(2*ms2+mc2)*trsgc-mc2*ms2**2*(8*trs-19*trg+10*
     &   trc))+mc2*cDg*gDs**2*(2*gDs+ms2)*(12*(mc-ms)*(ms+mc)*(trs+trg+t
     &   rc)*gDs**3+3*cDs*(4*(mc-ms)*(ms+mc)*(trs+trg+trc)*gDs**2+2*mc2*
     &   (-3*trsgc-2*ms2*trs+mc2*(-trs-3*trg+3*trc))*gDs+3*mc2**2*(trsgc
     &   +ms2*trg))+2*(mc2*(6*ms2*(trg+trc)-9*trsgc)+2*mc2**2*(7*trs+2*t
     &   rg+11*trc)-3*ms2**2*(trs+trg+trc))*gDs**2+mc2*(-mc2*(45*trsgc+m
     &   s2*(-12*trs+49*trg+22*trc))+10*mc2**2*trs+6*ms2**2*(-trs+trg+tr
     &   c))*gDs+9*mc2**2*ms2*(trsgc+mc2*trg))-mc2**2*gDs**3*(2*gDs+ms2)
     &   *(-6*(mc-ms)*(ms+mc)*(trs+trg+trc)*gDs**2+cDs*(9*mc2*(trsgc+mc2
     &   *trg)-6*(mc-ms)*(ms+mc)*(trs+trg+trc)*gDs)+(3*mc2*(3*trsgc-2*ms
     &   2*(trs+trg+trc))+3*ms2**2*(trs+trg+trc)+mc2**2*(-7*trs+4*trg-5*
     &   trc))*gDs+9*mc2**2*(trsgc+ms2*trg))+8*cDg**6*gDs*(16*trc*gDs**2
     &   +8*(-trs-trg+trc)*cDs*(2*gDs+ms2)+(2*mc2*(4*(trs+trg)+trc)-ms2*
     &   (3*(trs+trg)+2*trc))*gDs+ms2*(mc2*(4*(trs+trg)+trc)-5*ms2*trc))
     &   +64*(-trs-trg+trc)*cDg**7*gDs*(2*gDs+ms2))/(mc2*cDg**2*(2*cDg+m
     &   c2)*gDs**3*(2*gDs+ms2))/1.2e+1_dp+Ammp

      Ammp = B0cgsf*(-mc2*cDg**3*((2*gDs+ms2)*(-2*(2*trs+trg)*gDs**2+cD
     &   s*(2*(trc-trs)*gDs+ms2*trg)+(2*ms2*trs+3*ms2*trg+mc2*trg-2*ms2*
     &   trc)*gDs+mc2*ms2*trg)+trsgc*(2*gDs*(2*gDs+cDs)+ms2*(2*gDs+cDs)+
     &   ms2**2))+mc2*cDg**2*gDs*(2*(-trsgc+mc2*(trs+trg+trc)-2*ms2*trs)
     &   *gDs**2+(mc2*(2*trg*cDs+trsgc-2*ms2*(trs-trg+trc))-trsgc*(2*cDs
     &   +ms2))*gDs+ms2*(mc2*trg*(cDs+ms2)+trsgc*(cDs+mc2)))+mc2*cDg*gDs
     &   **2*(2*(mc2*(2*trs+trg+2*trc)-2*ms2*trs)*gDs**2+cDs*(2*(mc2*(tr
     &   s+trc)-2*ms2*trs)*gDs+mc2*(trsgc+ms2*trg))-(-2*ms2**2*(-trs+trg
     &   +trc)+mc2*ms2*(-4*trs+5*trg+4*trc)+mc2**2*trg)*gDs+mc2*ms2*(trs
     &   gc+mc2*trg))+mc2*gDs**3*(2*(mc-ms)*(ms+mc)*(trs+trg+trc)*gDs**2
     &   +cDs*(2*(mc-ms)*(ms+mc)*(trs+trg+trc)*gDs-mc2*(trsgc+mc2*trg))-
     &   (mc2*(trsgc-2*ms2*(trs+trg+trc))+ms2**2*(trs+trg+trc)+mc2**2*tr
     &   g)*gDs-mc2**2*(trsgc+ms2*trg))+mc2*cDg**4*(2*gDs+ms2)*(2*(2*trs
     &   +trg)*gDs-trsgc+ms2*(trs+trc))+4*(-trs-trg+trc)*cDg**5*gDs*(gDs
     &   +cDs-mc2)+4*(-trs-trg+trc)*cDg**6*gDs)/(mc2*cDg**2*gDs**3)/4.0d
     &   +0+Ammp

      Ammp = BfunX*(2*mc2*cDg**3*gDs*(trs*(6*gDs**2+2*cDs*(2*gDs+ms2)+3
     &   *ms2*gDs+ms2**2)-(2*gDs+ms2)*((trc-trg)*(gDs+ms2)-2*trg*cDs))+2
     &   *mc2*cDg*gDs**3*(2*(gDs+cDs)+ms2)*(trs*(3*gDs+2*cDs+ms2)+(trg+t
     &   rc)*(gDs-ms2))+4*mc2*cDg**2*gDs**3*(4*trs*(gDs+cDs)-ms2*(-2*trs
     &   +trg+trc))+mc2*(trs+trg+trc)*gDs**4*(2*(gDs+cDs)+ms2)**2-mc2*cD
     &   g**4*(2*gDs+ms2)*(2*(-trs-trg+trc)*gDs+ms2*(trs+trg+trc))+4*(tr
     &   s+trg-trc)*cDg**5*gDs*(gDs+cDs+mc2)+4*(trs+trg-trc)*cDg**6*gDs)
     &   /(mc2*cDg**2*gDs**3)/4.0_dp+Ammp

      Ammp = 3*mc2*tr3s002ft*((trg+trc)*((2*cDg+mc2)*gDs**2+ms2*cDg**2)
     &   +2*trs*cDg*(2*cDg+mc2)*gDs)/(cDg**2*gDs)+Ammp

      Ammp = epinv**2*(cDg**2*(4*trs*gDs**2-2*ms2*trg*gDs+ms2*trsgc)+mc
     &   2*trsgc*gDs**2-2*trsgc*cDg*gDs*(gDs+cDs))/(cDg*gDs**2)/2.0_dp+
     &   Ammp

      Ammp = mc2*tr3s00ft*(-cDg*(2*(trc-2*trs)*gDs**2+cDs*(2*(trc-trs)*
     &   gDs+trsgc+ms2*trg)-(2*trsgc+7*ms2*trg+6*ms2*trc)*gDs+ms2*(trsgc
     &   +mc2*trg))+mc2*gDs*(trg*(gDs+cDs+ms2)+trs*gDs)+cDg**2*(4*trs*gD
     &   s+ms2*(3*trs+trg+2*trc))+trsgc*gDs*(gDs+cDs+mc2))/(cDg**2*gDs)+
     &   Ammp

      Ammp = tr3c00fs*(-mc2*gDs*(mc2*(trg*(ms2*(6*gDs+cDs)+2*gDs*(3*gDs
     &   +cDs)+ms2**2)+ms2*trs*gDs)+trsgc*(2*(2*cDs+mc2)*gDs+ms2*(cDs+mc
     &   2)))+mc2*cDg*((2*gDs+ms2)*(cDs*(2*(trc-trs)*gDs+ms2*trg)-(gDs+m
     &   s2)*(6*trc*gDs-mc2*trg))+trsgc*(ms2*(4*gDs+cDs)+2*gDs*(gDs+cDs)
     &   +ms2**2))+12*cDg**3*gDs*(trc*(gDs+ms2)+(trs+trg)*cDs)-mc2*cDg**
     &   2*(2*gDs+ms2)*(2*(4*trs+3*trg+trc)*gDs-trsgc+ms2*(3*trs+2*(trg+
     &   trc)))+12*(trs+trg)*cDg**4*gDs)/(mc2*gDs**3)+Ammp

      Ammp = 3*tr3c001fs*(4*cDg**3*gDs*(2*trc*gDs+(trs+trg)*cDs+ms2*trc
     &   )-2*mc2*trc*cDg*gDs*(2*gDs+ms2)**2-mc2*(trs+trg)*cDg**2*(2*gDs+
     &   ms2)**2-mc2**2*(trs+trg)*gDs**2*(2*gDs+ms2)+4*(trs+trg)*cDg**4*
     &   gDs)/(mc2*gDs**3)+Ammp

      Ammp = Ammp-3*tr3c002fs*(mc2*trc*cDg**2*(2*gDs+ms2)**2+mc2**2*tr
     &   c*gDs**2*(2*gDs+ms2)+2*mc2**2*(trs+trg)*cDg*gDs*(2*gDs+ms2)-4*c
     &   Dg**3*(trc*cDs+mc2*(trs+trg))*gDs-4*trc*cDg**4*gDs)/(mc2*gDs**3
     &   )

      Ammp = ls*(-8*(trs+trg)*cDg**3*gDs**2+mc2*gDs*(2*gDs+ms2)*(trsgc*
     &   (4*gDs+ms2)-mc2*trg*(2*gDs+ms2))+mc2*cDg*(2*gDs+ms2)**2*(2*(-tr
     &   s+trg+trc)*gDs-trsgc+ms2*trg))/(mc2*gDs*(2*gDs+ms2)**2)/2.0_dp+
     &   Ammp

      Ammp = 3*mc2*tr3s001ft*(mc2*trs*gDs**2+2*cDg*gDs*(trs*gDs+ms2*(tr
     &   g+trc))+ms2*trs*cDg**2)/(cDg**2*gDs)+Ammp

      Ammp = lc*mc2*(mc2*cDg*(ms2*trg-2*(trs-trg+trc)*gDs)+2*cDg**2*(2*
     &   (trg-trs)*gDs+trsgc+ms2*trg)+trsgc*cDg*(2*gDs+mc2)+mc2*(trsgc+m
     &   c2*trg)*gDs)/((2*cDg+mc2)**2*gDs)/2.0_dp+Ammp

      Ammp=Ammp/zmmp

      return
      end
