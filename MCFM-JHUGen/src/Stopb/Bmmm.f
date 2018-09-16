      subroutine Bamp_mmm(q,mc,ms,Bmmm)
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
      complex(dp):: trc,trg,trs,trsgc,zmmm,Bmmm

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
      zmmm=zb(2,4)*zb(2,3)

      Bmmm = ms*ms2*(cDg*gDs*(mc2*(2*trs*gDs+trsgc-ms2*trg)+cDs*(-10*tr
     &   c*gDs+3*trsgc+3*mc2*trg)+2*trg*cDs**2)+gDs**2*(mc2*(2*trs*gDs+t
     &   rsgc-3*ms2*trg)+cDs*(-4*trc*gDs+trsgc+3*mc2*trg)+4*trg*cDs**2)+
     &   2*cDg**2*cDs*(trsgc-3*trc*gDs))*tr5Xs/((cDs**2-mc2*ms2)*gDs**2)
     &   /2.0_dp

      Bmmm = mc2*ms*(2*mc2*trg*cDs*gDs**2-cDg**2*(cDs*(2*(trg+trc)*gDs-
     &   trsgc+mc2*trg)-2*(mc2*trs+ms2*trg)*gDs+mc2*(trsgc+ms2*trg))+cDg
     &   *gDs*(cDs*(-2*(trg+trc)*gDs+trsgc+mc2*trg)-2*trg*cDs**2+mc2*(ms
     &   2*trg-trsgc))+2*(mc2*trs+ms2*trg)*cDg**3)*tr5Xc/(cDg**2*(cDs**2
     &   -mc2*ms2))/2.0_dp+Bmmm

      Bmmm = ms*ms2*(gDs*(-mc2*(2*trs*gDs+trsgc-3*ms2*trg)-cDs*(-4*trc*
     &   gDs+trsgc+3*mc2*trg)-4*trg*cDs**2)+cDg*cDs*(6*trc*gDs-2*trsgc))
     &   *tr4Xs/((cDs**2-mc2*ms2)*gDs)/2.0_dp+Bmmm

      Bmmm = mc2*ms*(cDg*(cDs*(2*(trg+trc)*gDs-trsgc+mc2*trg)+mc2*(trsg
     &   c+ms2*trg))-2*mc2*trg*cDs*gDs-2*(mc2*trs+ms2*trg)*cDg**2)*tr4Xc
     &   /(cDg*(cDs**2-mc2*ms2))/2.0_dp+Bmmm

      Bmmm = ms*(cDg*(2*(mc2*ms2*(2*trg+trc)-trg*cDs**2)*gDs+mc2*(trsgc
     &   *(cDs-ms2)+ms2*trg*(cDs+mc2)))-2*mc2**2*ms2*trg*gDs-2*(mc2*trs+
     &   ms2*trg)*cDg**2*cDs)*tr3Xs/(cDg*(cDs-mc*ms)*(cDs+mc*ms))/2.0_dp
     &   +Bmmm

      Bmmm = Bmmm-mc2*ms*(gDs*(2*trs*cDs*gDs-4*ms2*trc*gDs+trsgc*cDs+m
     &   s2*trg*cDs+ms2*trsgc+3*mc2*ms2*trg)+2*ms2*cDg*(trsgc-3*trc*gDs)
     &   )*tr3Xc/((cDs**2-mc**2*ms**2)*gDs)/2.0_dp

      Bmmm = ms*((cDs-mc*ms)*(cDs+mc*ms)*gDs**2*(trs*(2*(gDs+cDs)-mc2)-
     &   ms2*(-trs+trg+trc))+cDg**2*(-mc2*ms2*(2*(-trs-trg+trc)*gDs-3*sc
     &   k*trsgc+2*trsgc+ms2*trs+2*ms2*trg)+cDs**2*(-4*trg*gDs-3*sck*trs
     &   gc+2*trsgc+ms2*trs+2*ms2*trg)-2*(trsgc-mc2*trs+ms2*trc)*cDs*gDs
     &   )+cDg*(cDs-mc*ms)*(cDs+mc*ms)*gDs*(2*(2*trs+trg)*gDs-ms2*(trs+t
     &   rg+trc))+4*trs*cDg**3*cDs*gDs)*lVs/(cDg**2*(cDs-mc*ms)*(cDs+mc*
     &   ms)*gDs)/4.0_dp+Bmmm

      Bmmm = (mc2*(mc2*ms2-cDs**2)*gDs**2*(trs*(2*(gDs+cDs)-mc2)-ms2*(-
     &   trs+trg+trc))+cDg**2*(mc2*ms2*(ms2*(2*trsgc+mc2*(trs+trg))-2*(t
     &   rsgc+2*ms2*trg+3*ms2*trc)*gDs)+cDs**2*(2*mc2*trs*gDs+4*ms2*(trg
     &   +trc)*gDs-ms2*(2*trsgc+mc2*(trs+trg)))-2*mc2*ms2*(-trs+trg+trc)
     &   *cDs*gDs)+cDg*(cDs-mc*ms)*(cDs+mc*ms)*gDs*(mc2*(-4*trs*gDs-trsg
     &   c+ms2*(trs+trc))+2*trg*gDs*(gDs+cDs))+4*mc2*ms2*trs*cDg**3*gDs)
     &   *lVc/(ms*cDg**2*(cDs-mc*ms)*(cDs+mc*ms)*gDs)/4.0_dp+Bmmm

      Bmmm = xs*cDs*(cDg*(trsgc-2*(trg+trc)*gDs)+mc2*trg*gDs)*lRcs/(mc*
     &   (xs**2-1)*cDg*gDs)+Bmmm

      Bmmm = (cDg**2*((2*gDs+ms2)*(16*trs*gDs**3+2*ms2*(3*trs+9*trg+47*
     &   trc)*gDs**2-ms2*(18*(trs+2*trg)*cDs+ms2*(3*(trs+trg)-44*trc)+mc
     &   2*(5*trs+23*trg))*gDs-9*ms2**2*trg*(cDs+mc2))-3*ms2*trsgc*(2*(5
     &   *sck+6)*gDs**2+3*cDs*(2*gDs+ms2)+ms2*(5*sck+13)*gDs+3*ms2**2))-
     &   gDs**2*(2*gDs+ms2)*(18*trg*gDs**3+6*(6*trg*cDs+(ms2-mc2)*trs)*g
     &   Ds**2+(18*trg*cDs**2+mc2*(-6*trs*cDs-9*trsgc+ms2*(-6*trs+2*trg+
     &   29*trc))+6*ms2*trs*cDs-5*mc2**2*trs-3*ms2**2*(-trs+trg+trc))*gD
     &   s-9*mc2*(trsgc*(cDs+ms2)+ms2*trg*(cDs+mc2)))+cDg*gDs**2*(2*gDs+
     &   ms2)*(mc2*(ms2*(5*(trg+trc)-11*trs)-6*trs*(4*gDs+3*cDs))+16*trs
     &   *gDs*(gDs+cDs)+12*ms2*(3*trg-trs)*gDs+3*ms2**2*(trs+trg+trc))+m
     &   s2*cDg**3*(2*gDs+ms2)*((20*trc-2*(trs+10*trg))*gDs-9*trsgc+ms2*
     &   (8*trs-trg+10*trc)))/(ms*cDg**2*gDs**2*(2*gDs+ms2))/1.2e+1_dp+Bmm
     &   m

      Bmmm = B0cgsf*(cDg*gDs**2*((2*trs-4*trg)*gDs**3+2*((trs-3*trg)*cD
     &   s+mc2*trs+ms2*(trg-3*trs))*gDs**2-2*(trg*cDs**2+ms2*trs*cDs+mc2
     &   *(2*ms2*trc-trsgc)-ms2**2*(trg+trc))*gDs+mc2*(trsgc*(cDs+ms2)+m
     &   s2*trg*(cDs+mc2)))+cDg**3*(2*trs*gDs**3+(4*ms2*trc-2*mc2*trs)*g
     &   Ds**2-2*ms2*(trs*cDs+2*trg*cDs+trsgc+ms2*(trg-3*trc)+mc2*trg)*g
     &   Ds-ms2*((trsgc+ms2*trg)*cDs+ms2*(trsgc+mc2*trg)))-cDg**2*gDs*(2
     &   *(trg-2*trs)*gDs**3+2*((trg-trs)*cDs+mc2*trs-ms2*(-2*trs+2*trg+
     &   trc))*gDs**2-(mc2*(trsgc-2*trs*cDs)+ms2*(-2*trs*cDs-trsgc+mc2*(
     &   trg-2*trs))+ms2**2*(trg+6*trc))*gDs+ms2*trsgc*(cDs+ms2)+ms2**2*
     &   trg*(cDs+mc2))+gDs**3*(-2*trg*gDs**3-2*(2*trg*cDs+(ms2-mc2)*trs
     &   )*gDs**2+(-2*trg*cDs**2+mc2*(2*trs*cDs+trsgc-ms2*(-2*trs+trg+4*
     &   trc))-2*ms2*trs*cDs+ms2**2*(-trs+trg+trc))*gDs+mc2*(trsgc*(cDs+
     &   ms2)+ms2*trg*(cDs+mc2)))+ms2*cDg**4*(2*(trc-trg)*gDs-trsgc+ms2*
     &   (trs+trc)))/(ms*cDg**2*gDs**2*(gDs+cDg))/4.0_dp+Bmmm

      Bmmm = tr3s00ft*(2*trg*gDs**3+(4*trg*cDs-6*trs*cDg)*gDs**2+(2*trg
     &   *(cDs**2-2*ms2*cDg)-6*trs*cDg*(cDs+cDg)+mc2*(4*trs*cDg-trsgc+2*
     &   ms2*trg+5*ms2*trc))*gDs-mc2*(cDg*(ms2*(trg-3*trs)-2*trs*cDs)+tr
     &   sgc*(cDs+ms2)+ms2*trg*(cDs+mc2)))/(ms*cDg**2)+Bmmm

      Bmmm = BfunX*(cDg**2*gDs*(6*trs*gDs**2-2*ms2*(trs+trg+2*trc)*gDs+
     &   ms2*(2*(trs+trg)*cDs+ms2*(trs+trg-trc)))-cDg*gDs**2*(trs*(-10*g
     &   Ds*(gDs+cDs)+2*ms2*(cDs-2*gDs)+ms2**2)+ms2*(trg+trc)*(4*gDs+2*c
     &   Ds+ms2))+gDs**3*(2*(gDs+cDs)+ms2)*(2*trs*(gDs+cDs)-ms2*(-trs+tr
     &   g+trc))-ms2*cDg**3*(2*trc*gDs+ms2*(trs+trg+trc)))/(ms*cDg**2*gD
     &   s**2)/4.0_dp+Bmmm

      Bmmm = LsB2*ms*(-2*mc2**2*trg*cDs*gDs**2-cDg**2*(mc2*(2*mc2*trs*g
     &   Ds-2*ms2*trc*gDs+ms2*trsgc-mc2*ms2*trg)+2*cDs**2*(2*(trg+trc)*g
     &   Ds-trsgc+mc2*trg)+mc2*(trsgc+ms2*trg)*cDs)+mc2*cDg*gDs*(cDs*(2*
     &   (trg+trc)*gDs-trsgc+mc2*trg)+4*trg*cDs**2+mc2*(trsgc-ms2*trg))+
     &   2*(mc2*trs+ms2*trg)*cDg**3*cDs)/(cDg**2*(cDs**2-mc2*ms2))+Bmmm

      Bmmm = lc*(-2*trg*gDs**2-2*trg*cDs*gDs+mc2*(-2*trs*cDg+trsgc+ms2*
     &   trg))/(ms*(2*cDg+mc2))/2.0_dp+Bmmm

      Bmmm = ms*tr3c00fs*(-12*trc*gDs**2-cDg*(2*(2*trs+trg+2*trc)*gDs-t
     &   rsgc+ms2*(3*trs+2*(trg+trc)))+(2*(trs+2*trg)*cDs+3*trsgc+2*mc2*
     &   trg-7*ms2*trc)*gDs+trsgc*(cDs+ms2)+ms2*trg*(cDs+mc2))/gDs**2+B
     &   mmm

      Bmmm = LsB1*ms*(-cDg*gDs*(2*cDs**2*(2*(trs+2*trc)*gDs-trsgc+ms2*t
     &   rg)+mc2*ms2*(-2*(trs+trc)*gDs+trsgc-ms2*trg)+ms2*cDs*(-4*trc*gD
     &   s+trsgc+3*mc2*trg))+gDs**2*(mc2*cDs*(2*trs*gDs+trsgc-7*ms2*trg)
     &   +4*cDs**2*(-2*trc*gDs+trsgc+mc2*trg)-mc2*ms2*(-4*trc*gDs+3*trsg
     &   c+mc2*trg)+8*trg*cDs**3)+2*ms2*cDg**2*cDs*(3*trc*gDs-trsgc))/((
     &   cDs**2-mc2*ms2)*gDs**2)+Bmmm

      Bmmm = B0csf*ms*(cDg*(mc2*(-2*trs*gDs+trsgc-mc2*trs+ms2*(trg-trs)
     &   )+cDs*(-2*trs*gDs+trsgc+mc2*(-2*trs+trg+trc)+ms2*trc)+2*trc*cDs
     &   **2)+mc2*((trsgc-ms2*trs-mc2*trs+ms2*trg)*gDs+ms2*(ms2+mc2)*trg
     &   )+cDs*((trsgc+mc2*(-2*trs+trg+trc)+ms2*trc)*gDs+2*mc2*ms2*trg)+
     &   cDs**2*(2*trc*gDs-(ms2+mc2)*trg)-2*trg*cDs**3-2*trs*cDg**2*(cDs
     &   +mc2))/((cDs**2-mc2*ms2)*(gDs+cDg))/2.0_dp+Bmmm

      Bmmm = tr3s001ft*(3*mc2*ms2*(trg+trc)*gDs-3*trs*cDg*(2*gDs*(gDs+c
     &   Ds)-mc2*ms2))/(ms*cDg**2)+Bmmm

      Bmmm = epinv*(cDg*(2*(trg+trc)*(mc*ms*(xs**2-1)-4*lp*xs*cDs)*gDs+
     &   trsgc*(4*lp*xs*cDs-mc*ms*(2*sck-1)*(xs**2-1)))+mc2*trg*(4*lp*xs
     &   *cDs+mc*(ms-ms*xs**2))*gDs)/(mc*(xs**2-1)*cDg*gDs)/4.0_dp+Bmmm

      Bmmm = ls*(ms*(2*gDs+ms2)*(-2*trs*gDs+4*trg*gDs+3*ms2*trg)-ms*trs
     &   gc*(4*gDs+ms2))/(2*gDs+ms2)**2/2.0_dp+Bmmm

      Bmmm = ms*tr2fu*cDs*(cDg*(trsgc-2*(trg+trc)*gDs)+mc2*trg*gDs)/(cD
     &   g*gDs)+Bmmm

      Bmmm = Bmmm-3*ms*tr3c002fs*(trc*cDg*(2*gDs+ms2)+mc2*(trs+trg)*gD
     &   s)/gDs**2

      Bmmm = 3*tr3s002ft*(trs*(mc2**2-4*cDg**2)*gDs+mc2*ms2*(trg+trc)*c
     &   Dg)/(ms*cDg**2)+Bmmm

      Bmmm = Bmmm-3*ms*tr3c001fs*(2*gDs+ms2)*(trc*gDs+(trs+trg)*cDg)/g
     &   Ds**2

      Bmmm=Bmmm/zmmm

      return
      end
