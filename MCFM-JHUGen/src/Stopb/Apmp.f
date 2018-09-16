      subroutine Aamp_pmp(q,mc,ms,Apmp)
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
      complex(dp):: trc,trg,trs,trsgc,zpmp,Apmp

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

      Apmp = ms2*(trsgc*(gDs*(cDg*(gDs+2*cDs+cDg)-mc2*gDs)-ms2*cDg**2)+
     &   cDg*gDs*(cDg*(ms2*trg-2*(trs+trc)*gDs)+mc2*trg*gDs))*tr3Xs/(cDg
     &   *gDs**2)

      Apmp = (mc2*trsgc*(-mc2*gDs**2+cDg**2*(gDs-ms2)+2*cDg*cDs*gDs)-cD
     &   g*gDs*(trg*(-mc2**2*gDs+2*cDg**2*cDs+cDg**3)+2*mc2*trc*cDg*gDs)
     &   )*tr3Xc/cDg**3+Apmp

      Apmp = Apmp-(cDg**2*(4*trc*gDs**2-2*trsgc*gDs+ms2*trsgc)+mc2*trs
     &   gc*gDs**2-2*cDg*gDs*(mc2*trg*gDs+trsgc*cDs))*tr1Xfs/cDg**2

      Apmp = (-4*trc*cDg*gDs-mc2*trsgc*gDs/cDg+2*mc2*trg*gDs-ms2*trsgc*
     &   cDg/gDs+2*trsgc*cDs+2*trsgc*cDg)*tr1Xfc+Apmp

      Apmp = (cDg**3*(8*mc2*trc*gDs**2+(4*trg*cDs**2+2*mc2*(ms2*(2*trs+
     &   trg+(4-3*sck)*trc)-2*trsgc))*gDs-mc2*ms2*(-3*sck*trsgc+2*trsgc+
     &   2*ms2*trs+3*ms2*trg))+2*mc2*cDg*gDs**2*(ms2*((trg+trc)*(gDs+ms2
     &   )-trs*(gDs+2*cDs+ms2-mc2))+mc2*trsgc)-mc2*cDg**2*gDs*(-2*(ms2*(
     &   trs+trg+trc)-2*mc2*trg)*gDs+4*trsgc*cDs+ms2*(trsgc+ms2*(trg+trc
     &   )+mc2*(4-3*sck)*trg))-mc2*ms2*(trs+trg+trc)*gDs**3*(2*gDs+2*cDs
     &   +ms2-mc2)+6*trg*cDg**4*cDs*gDs+2*trg*cDg**5*gDs)*lVs/(mc2*cDg**
     &   3*gDs)/4.0_dp+Apmp

      Apmp = (2*cDg**3*(2*(trs+trg+trc)*gDs**2+trg*(2*cDs-ms2)*gDs+ms2*
     &   (trsgc+mc2*(trs+trg)))+cDg**2*gDs*(mc2*trc*(2*gDs+ms2)-trsgc*(2
     &   *gDs+4*cDs+mc2))+mc2*cDg*gDs**2*(-2*(-trs+trg+trc)*gDs+4*trs*cD
     &   s+trsgc-mc2*(2*trs+trg)-2*ms2*(-trs+trg+trc))+mc2*(trs+trg+trc)
     &   *gDs**3*(2*gDs+2*cDs+ms2-mc2)+2*trg*cDg**4*gDs)*lVc/(cDg**3*gDs
     &   )/4.0_dp+Apmp

      Apmp = ((mc2*(2*trg*cDg-trsgc)/cDg**2/2.0_dp-2*trc)*gDs-ms2*trsgc
     &   /gDs/2.0_dp+trsgc*cDs/cDg+trsgc)*lRs2+Apmp

      Apmp = epinv*(trsgc*(cDg**2*(gDs*(4*lRs1+4*lRc1-2)-ms2*(2*lRs1+2*
     &   lRc1-2*sck+1))+2*cDg*cDs*gDs*(2*lRs1+2*lRc1-1)+mc2*gDs**2*(-2*l
     &   Rs1-2*lRc1+1))-2*cDg*(2*trc*cDg-mc2*trg)*gDs*(gDs*(2*lRs1+2*lRc
     &   1-1)+ms2*(sck-1)))/(cDg**2*gDs)/4.0_dp+Apmp

      Apmp = ((mc2*(2*trg*cDg-trsgc)/cDg**2/2.0_dp-2*trc)*gDs-ms2*trsgc
     &   /gDs/2.0_dp+trsgc*cDs/cDg+trsgc)*lRc2+Apmp

      Apmp = (cDg**3*(32*mc2*(2*trs+3*trc)*gDs**4+4*mc2*(36*trg*cDs-12*
     &   trsgc+ms2*(16*trs-12*trg+(41-15*sck)*trc))*gDs**3+(-72*trg*cDs*
     &   *3+36*mc2*ms2*(-trs+4*trg+trc)*cDs+2*mc2*ms2*(15*sck*trsgc+ms2*
     &   (11*trs-9*trg-15*sck*trc+13*trc)+10*mc2*trs+19*mc2*trg))*gDs**2
     &   +ms2*(-36*trg*cDs**3+18*mc2*(trsgc+ms2*(-trs+3*trg+trc))*cDs+mc
     &   2*ms2*(15*(sck+2)*trsgc+ms2*(6*(trs+trg)-16*trc)+mc2*(10*trs+37
     &   *trg)))*gDs+9*mc2*ms2**2*(trsgc*(cDs+ms2)+ms2*trg*(cDs+mc2)))-m
     &   c2*cDg**2*gDs*(2*gDs+ms2)*(4*(trs+8*trg+17*trc)*gDs**3+cDs*(36*
     &   (trc-trs)*gDs**2+(54*trsgc+26*ms2*trg+8*ms2*trc)*gDs+9*ms2*(trs
     &   gc+mc2*trg))+2*(-9*trsgc+5*mc2*trs-3*ms2*(-trs+6*trg+3*trc)+12*
     &   mc2*trg)*gDs**2-ms2*(-18*trsgc+mc2*(18*trs+(15*sck-44)*trg+23*t
     &   rc)+3*ms2*(trg+trc))*gDs+9*mc2*ms2*(trsgc+ms2*trg))-mc2*cDg*gDs
     &   **2*(2*gDs+ms2)*(2*(-9*trsgc+2*mc2*(7*trs+2*trg+11*trc)+3*ms2*(
     &   -trs+trg+trc))*gDs**2+3*cDs*(3*mc2*(trsgc+ms2*trg)-2*(3*trsgc+m
     &   c2*(trs+3*trg-3*trc)+2*ms2*trs)*gDs)+(-mc2*(45*trsgc+ms2*(-12*t
     &   rs+49*trg+22*trc))+10*mc2**2*trs+6*ms2**2*(-trs+trg+trc))*gDs+9
     &   *mc2*ms2*(trsgc+mc2*trg))+mc2*gDs**3*(2*gDs+ms2)*(-6*(mc-ms)*(m
     &   s+mc)*(trs+trg+trc)*gDs**2+cDs*(9*mc2*(trsgc+mc2*trg)-6*(mc-ms)
     &   *(ms+mc)*(trs+trg+trc)*gDs)+(3*mc2*(3*trsgc-2*ms2*(trs+trg+trc)
     &   )+3*ms2**2*(trs+trg+trc)+mc2**2*(-7*trs+4*trg-5*trc))*gDs+9*mc2
     &   **2*(trsgc+ms2*trg))+cDg**4*(2*gDs+ms2)*(36*mc2*trg*gDs**2-2*(4
     &   5*trg*cDs**2+mc2*ms2*(17*trs-trg+trc))*gDs+mc2*ms2*(9*trsgc+ms2
     &   *(-8*trs+trg-10*trc)))-72*trg*cDg**5*cDs*gDs*(2*gDs+ms2)-18*trg
     &   *cDg**6*gDs*(2*gDs+ms2))/(mc2*cDg**3*gDs**2*(2*gDs+ms2))/1.2e+1_dp
     &   +Apmp

      Apmp = 2*LsA*(mc2**2*trsgc*gDs**4-cDg**3*gDs*(cDs*(4*trc*gDs**2-2
     &   *trsgc*gDs+3*ms2*trsgc)+ms2*(trsgc+mc2*trg)*gDs)+ms2*cDg**4*(2*
     &   (trs+trc)*gDs**2-(trsgc+ms2*trg)*gDs+ms2*trsgc)+cDg**2*gDs**2*(
     &   2*mc2*trc*gDs**2+mc2*(2*trg*cDs-trsgc)*gDs+2*trsgc*(cDs**2+mc2*
     &   ms2))-mc2*cDg*gDs**3*(mc2*trg*gDs+3*trsgc*cDs))/(cDg**3*gDs**2)
     &   +Apmp

      Apmp = B0cgsf*(cDg**3*(-4*mc2*(trg+2*trc)*gDs**3+(mc2*(2*(2*trg*c
     &   Ds+trsgc)+ms2*(-2*trs-3*trg+trc))-4*trg*cDs**2)*gDs**2+(-4*trg*
     &   cDs**3+2*mc2*ms2*(-trs+2*trg+trc)*cDs+mc2*ms2*(ms2*(2*trs+3*trg
     &   -2*trc)+mc2*trg))*gDs+mc2*ms2*(trsgc*(cDs+ms2)+ms2*trg*(cDs+mc2
     &   )))-mc2*cDg**2*gDs*(4*(trg+2*trc)*gDs**3+cDs*(4*(trc-trs)*gDs**
     &   2-2*trsgc*gDs+ms2*(3*trg+trc)*gDs+ms2*(trsgc+mc2*trg))+(-4*trsg
     &   c+2*mc2*(trs-trg+trc)-ms2*(-2*trs+3*trg+trc))*gDs**2-(mc2*(trsg
     &   c+2*ms2*(trs-trg+trc))+ms2*(ms2*(trg+trc)-trsgc))*gDs+mc2*ms2*(
     &   trsgc+ms2*trg))-cDg**4*(mc2*(-2*trg*gDs**2+ms2*(4*trs*gDs-trsgc
     &   )+ms2**2*(trs+trc))+2*trg*cDs*gDs*(3*gDs+5*cDs))-mc2*cDg*gDs**2
     &   *(2*(-trsgc+2*mc2*(trs+trc)+ms2*(-trs+trg+trc))*gDs**2+cDs*(2*(
     &   -trsgc+mc2*(trs-trg+trc)-2*ms2*trs)*gDs+mc2*(trsgc+ms2*trg))-(m
     &   c2*(2*trsgc+ms2*(-4*trs+7*trg+4*trc))-2*ms2**2*(-trs+trg+trc)+m
     &   c2**2*trg)*gDs+mc2*ms2*(trsgc+mc2*trg))+mc2*gDs**3*(-2*(mc-ms)*
     &   (ms+mc)*(trs+trg+trc)*gDs**2+cDs*(mc2*(trsgc+mc2*trg)-2*(mc-ms)
     &   *(ms+mc)*(trs+trg+trc)*gDs)+(mc2*(trsgc-2*ms2*(trs+trg+trc))+ms
     &   2**2*(trs+trg+trc)+mc2**2*trg)*gDs+mc2**2*(trsgc+ms2*trg))-2*tr
     &   g*cDg**5*gDs*(gDs+4*cDs)-2*trg*cDg**6*gDs)/(mc2*cDg**3*gDs**2)/
     &   4.0_dp+Apmp

      Apmp = tr3c00fs*(cDg*(-mc2*cDs*(8*trg*gDs**2+2*ms2*(-trs+2*trg+tr
     &   c)*gDs+ms2*(trsgc+ms2*trg))-mc2*ms2*(-8*trc*gDs**2+(2*trsgc+mc2
     &   *trg-6*ms2*trc)*gDs+ms2*(trsgc+mc2*trg))+4*trg*cDs**3*gDs)+cDg*
     &   *2*(-4*mc2*trg*gDs**2+2*(5*trg*cDs**2+mc2*ms2*(4*trs+2*trg+trc)
     &   )*gDs+mc2*ms2*(ms2*(3*trs+2*(trg+trc))-trsgc))+mc2*ms2*gDs*(mc2
     &   *(trg*(2*gDs+cDs+ms2)+trs*gDs)+trsgc*(cDs+mc2))+8*trg*cDg**3*cD
     &   s*gDs+2*trg*cDg**4*gDs)/(mc2*cDg*gDs**2)+Apmp

      Apmp = ls*(-cDg*((8*trg*cDs**2+4*mc2*ms2*(trg-trc))*gDs**2+2*ms2*
     &   (2*trg*cDs**2+mc2*(trsgc-ms2*(trs-trg+trc)))*gDs+mc2*ms2**2*(tr
     &   sgc+ms2*trg))-6*trg*cDg**2*cDs*gDs*(2*gDs+ms2)-2*trg*cDg**3*gDs
     &   *(2*gDs+ms2)-mc2*ms2*(trsgc+mc2*trg)*gDs*(2*gDs+ms2))/(mc2*cDg*
     &   (2*gDs+ms2)**2)/2.0_dp+Apmp

      Apmp = tr3s001ft*(3*ms2*(trg+trc)*cDg*(cDg*(cDs-3*gDs)-2*mc2*gDs)
     &   -3*trs*(2*cDg+mc2)*((2*cDg+mc2)*gDs**2+ms2*cDg**2))/cDg**3+Apm
     &   p

      Apmp = Apmp-tr3s00ft*(cDg*(2*(trsgc+mc2*(3*trs+trg-trc))*gDs**2-
     &   cDs*(-2*trsgc*gDs+2*mc2*(-trs-trg+trc)*gDs+mc2*(trsgc+ms2*trg))
     &   +mc2*(4*trsgc+9*ms2*trg+6*ms2*trc)*gDs-mc2*ms2*(trsgc+mc2*trg))
     &   +cDg**2*(-4*(trc-2*trs)*gDs**2-cDs*(4*(trc-trs)*gDs+4*trsgc+5*m
     &   s2*trg+3*ms2*trc)+(2*trsgc+6*mc2*trs+11*ms2*trg+9*ms2*trc)*gDs+
     &   ms2*(-2*trsgc+3*mc2*trs-3*mc2*trg+2*mc2*trc))+mc2**2*gDs*(trg*(
     &   gDs+cDs+ms2)+trs*gDs)+3*cDg**3*(4*trs*gDs+ms2*(2*trs-trg+trc))+
     &   mc2*trsgc*gDs*(gDs+cDs+mc2))/cDg**3

      Apmp = 3*ms2*tr3c001fs*(mc2*(trs+trg)*gDs**2+2*trc*cDg*gDs*(2*gDs
     &   +ms2)+(trs+trg)*cDg**2*(2*gDs+ms2))/(cDg*gDs**2)+Apmp

      Apmp = 3*ms2*tr3c002fs*(mc2*trc*gDs**2+trc*cDg**2*(2*gDs+ms2)+2*m
     &   c2*(trs+trg)*cDg*gDs)/(cDg*gDs**2)+Apmp

      Apmp = BfunX*(cDg**2*gDs**2*((trg+trc)*(4*gDs*(gDs+cDs)+ms2*(gDs-
     &   cDs)-ms2**2)-2*trs*gDs*(2*(gDs+cDs)+ms2))+2*cDg*gDs**3*(2*(gDs+
     &   cDs)+ms2)*(ms2*(-trs+trg+trc)-2*trs*(gDs+cDs))-(trs+trg+trc)*gD
     &   s**4*(2*(gDs+cDs)+ms2)**2-ms2*cDg**3*gDs*((4*trs+5*trg+trc)*gDs
     &   +4*(trs+trg)*cDs+2*ms2*(trs+trg-trc))+ms2*cDg**4*(2*(-trs-trg+t
     &   rc)*gDs+ms2*(trs+trg+trc)))/(cDg**3*gDs**2)/4.0_dp+Apmp

      Apmp = lc*((2*cDg+mc2)*(cDg*(2*(-trs-trg+trc)*gDs+ms2*trg)-mc2*tr
     &   g*gDs)-trsgc*(-mc2*gDs+cDg*(mc2-2*gDs)+4*cDg**2))/(cDg*(2*cDg+m
     &   c2))/2.0_dp+Apmp

      Apmp = epinv**2*((mc2*(2*trg*cDg-trsgc)/cDg**2/2.0_dp-2*trc)*gDs-
     &   ms2*trsgc/gDs/2.0_dp+trsgc*cDs/cDg+trsgc)+Apmp

      Apmp = Apmp-3*tr3s002ft*(2*cDg+mc2)**2*gDs*((trg+trc)*gDs+2*trs*
     &   cDg)/cDg**3

      Apmp=Apmp/zpmp

      return
      end
