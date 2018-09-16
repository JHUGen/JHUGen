      subroutine Aamp_mmm(q,mc,ms,Ammm)
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
      complex(dp):: trc,trg,trs,trsgc,zmmm,Ammm

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

      Ammm = ms*(65536*trs*cDg**2*(cDg+mc2)*gDs**3+ms2*gDs*(cDg**2*(-65
     &   536*trg*gDs**2+2*(trg+2*trc)*cDs*gDs+mc2*ms2*trg)+mc2**2*trg*gD
     &   s**2-mc2*cDg*gDs*((trg+2*trc)*gDs+2*trg*cDs)-ms2*(trg+2*trc)*cD
     &   g**3)+ms2*trsgc*cDg*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2))*tr3X
     &   s/(gDs**2*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2))

      Ammm = mc2*ms*(cDg*(trsgc-(trg+2*trc)*gDs)+mc2*trg*gDs)*tr3Xc/cDg
     &   **2+Ammm

      Ammm = ms*(cDg*(trsgc-2*(trg+trc)*gDs)+mc2*trg*gDs)*tr1Xfs/cDg+A
     &   mmm

      Ammm = ms*(cDg*(trsgc-2*(trg+trc)*gDs)+mc2*trg*gDs)*tr1Xfc/gDs+A
     &   mmm

      Ammm = Ammm-ms*(cDg**4*(131072*trs*gDs**2+2*ms2*(trs+trg)*gDs-ms
     &   2*((2-3*sck)*trsgc+ms2*(trs+2*trg)))+cDg*gDs**3*(trs*(2*(2*cDs-
     &   mc2)*(gDs+cDs)+ms2*(2*cDs+mc2))-ms2*(trg+trc)*(2*cDs-mc2))+mc2*
     &   gDs**4*(trs*(mc2-2*(gDs+cDs))+ms2*(-trs+trg+trc))+cDg**3*gDs*(2
     &   *ms2*((trs+2*trg)*cDs-(trs+65536*trg)*gDs)-2*cDs*(2*(trs+trg)*g
     &   Ds+(3*sck-2)*trsgc)+131072*mc2*trs*gDs+ms2**2*(trs+trg+trc))+cD
     &   g**2*gDs**2*(mc2*(2*(trs+trg)*gDs+(3*sck-2)*trsgc)-2*ms2*(trs*g
     &   Ds+mc2*trg)-2*cDs*(ms2*(2*trs+trg+trc)-2*trs*gDs)+ms2**2*(-trs+
     &   trg+trc)))*lVs/(cDg**2*gDs*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2
     &   ))/4.0_dp

      Ammm = (-cDg**4*(131072*trs*gDs**3+131072*trs*cDs*gDs**2-4*ms2**2
     &   *(trg+trc)*gDs+ms2**2*(2*trsgc+mc2*(trs+trg)))+cDg**3*gDs*(mc2*
     &   (-131072*trs*gDs**2-ms2*(2*trs*gDs+trsgc)+ms2**2*(trs+trc))+131
     &   072*ms2*trg*gDs**2+2*cDs*(ms2*(2*trsgc+mc2*(trs+trg))-4*(16384*
     &   mc2*trs+ms2*(trc-16383*trg))*gDs))+mc2**2*gDs**4*(trs*(mc2-2*(g
     &   Ds+cDs))+ms2*(-trs+trg+trc))+mc2*cDg**2*gDs**2*(2*cDs*(2*trs*gD
     &   s+trsgc-ms2*(2*trs+trc))+ms2*(-2*trs*gDs+4*(trg+trc)*gDs-2*trsg
     &   c+ms2*(-trs+trg+trc)-mc2*trg))+mc2*cDg*gDs**3*(mc2*(-2*trs*gDs-
     &   trsgc+ms2*(trs+trc))+2*trs*cDs*(2*gDs+ms2-mc2)+4*trs*cDs**2-2*m
     &   s2*(trg+trc)*cDs))*lVc/(ms*cDg**2*gDs*(mc2*gDs**2-2*cDg*cDs*gDs
     &   +ms2*cDg**2))/4.0_dp+Ammm

      Ammm = ms*(cDg*(trsgc-2*(trg+trc)*gDs)+mc2*trg*gDs)*lRs2/(cDg*gDs
     &   )/2.0_dp+Ammm

      Ammm = epinv*(ms*trsgc*cDg*(2*lRs1+2*lRc1-2*sck+1)+ms*(mc2*trg-2*
     &   (trg+trc)*cDg)*gDs*(2*lRs1+2*lRc1-1))/(cDg*gDs)/4.0_dp+Ammm

      Ammm = ms*(cDg*(trsgc-2*(trg+trc)*gDs)+mc2*trg*gDs)*lRc2/(cDg*gDs
     &   )/2.0_dp+Ammm

      Ammm = tr3s00ft*(cDg**3*(mc2*ms2*(393218*trs*gDs+ms2*(3*trs+trg+2
     &   *trc))-131072*gDs*(trs*(gDs+cDs)**2+2*ms2**2*trg))-cDg**2*(-131
     &   072*(ms2*trg-mc2*trs)*gDs**3+cDs*(-4*(65536*ms2*trg-65537*mc2*t
     &   rs)*gDs**2+2*mc2*ms2*(3*trs+trg+2*trc)*gDs+mc2*ms2*(trsgc+ms2*t
     &   rg))-131072*(ms2*trg-mc2*trs)*cDs**2*gDs+mc2*ms2*(trsgc-131072*
     &   mc2*trs+131068*ms2*trg-5*ms2*trc)*gDs+mc2*ms2**2*(trsgc+mc2*trg
     &   ))+mc2*cDg*gDs*(2*cDs*((trsgc-4*ms2*trg-5*ms2*trc)*gDs+ms2*(trs
     &   gc+mc2*trg))+mc2*gDs*(2*trs*gDs+ms2*(3*trs+trg+2*trc))+2*(trsgc
     &   +ms2*trg)*cDs**2)+mc2**2*gDs**2*(ms2*(5*trc*gDs-trg*(-4*gDs+cDs
     &   +mc2))-trsgc*(gDs+cDs+ms2))+262144*ms2*trs*cDg**4*gDs)/(ms*cDg*
     &   *2*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2))+Ammm

      Ammm = B0cgsf*(cDg**4*(131072*trs*gDs**4+cDs*(131072*trs*gDs**3-4
     &   *ms2*(trs+trc)*gDs**2+2*ms2*(trsgc-ms2*(trs+trg+trc))*gDs-ms2**
     &   2*(trsgc+ms2*trg))-131072*ms2*trs*gDs**3+2*ms2**2*(trs+trc)*gDs
     &   **2-ms2**2*(trsgc+ms2*(trs+2*trg-3*trc)+2*mc2*trg)*gDs-ms2**3*(
     &   trsgc+mc2*trg))+cDg**3*gDs*(2*cDs*(131072*trs*gDs**3-2*(ms2*(tr
     &   s+32768*trg+trc)-32768*mc2*trs)*gDs**2+ms2*(trsgc+ms2*(trs+2*tr
     &   g-3*trc)+2*mc2*trg)*gDs+ms2**2*(trsgc+mc2*trg))+gDs*(131072*trs
     &   *gDs**3-131072*(ms2*trg-mc2*trs)*gDs**2+2*ms2*(-ms2*trs-131070*
     &   mc2*trs+65536*ms2*trg+mc2*trc)*gDs+ms2**2*(ms2-mc2)*(trs+trg+tr
     &   c))+2*cDs**2*(65536*trs*gDs**2+2*ms2*trg*gDs+ms2*(trsgc+ms2*trg
     &   )))+cDg*gDs**3*(-2*cDs*(2*(mc-ms)*(ms+mc)*trs*gDs**2+(mc2*(trsg
     &   c+ms2*(2*trs-3*trg-4*trc))+ms2**2*(-trs+trg+trc))*gDs+mc2*ms2*(
     &   trsgc+mc2*trg))+mc2*gDs*(mc2*(2*trs*gDs+trsgc-ms2*(2*trs+trg+2*
     &   trc))+ms2*(ms2*(trs+trg+trc)-2*trs*gDs))-2*cDs**2*(mc2*(2*trs*g
     &   Ds+trsgc+ms2*trg)-2*ms2*trs*gDs))+cDg**2*gDs**3*(-131072*(ms2*t
     &   rg-mc2*trs)*gDs**2+cDs*(mc2*(262140*trs*gDs-2*trsgc+6*ms2*trs+4
     &   *ms2*trc)-2*ms2*(2*(65536*trg-trs)*gDs+ms2*(2*trs+trg+trc)))+2*
     &   ms2*(mc2*(2*trs+trc)-ms2*trs)*gDs-131072*(ms2*trg-mc2*trs)*cDs*
     &   *2+ms2*(-2*mc2**2*(65536*trs+trg)+mc2*ms2*(trs+131067*trg-trc)+
     &   ms2**2*(-trs+trg+trc)))+mc2*gDs**4*(2*(mc-ms)*(ms+mc)*trs*gDs**
     &   2+(mc2*(2*trs*cDs+trsgc+ms2*(2*trs-3*trg-4*trc))+ms2*(ms2*(-trs
     &   +trg+trc)-2*trs*cDs))*gDs+mc2*(trsgc*(cDs+ms2)+ms2*trg*(cDs+mc2
     &   )))+ms2**2*cDg**5*((trs+trc)*(2*gDs+ms2)-trsgc))/(ms*cDg**2*gDs
     &   **2*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2))/4.0_dp+Ammm

      Ammm = (cDg**3*gDs*(2*cDs*((2*gDs+ms2)*(1179648*trs*gDs**3-2*ms2*
     &   (3*(trs+5*trg)+29*trc)*gDs**2+ms2*(3*ms2*trs+5*mc2*trs+3*ms2*tr
     &   g+23*mc2*trg-26*ms2*trc)*gDs+9*mc2*ms2**2*trg)+3*ms2*trsgc*(2*(
     &   5*sck+6)*gDs**2+ms2*(5*sck+13)*gDs+3*ms2**2))+gDs*(2*gDs+ms2)*(
     &   1179648*trs*gDs**3+2*ms2*(-3*ms2*trs-1769462*mc2*trs+1179648*ms
     &   2*trg-mc2*trg+10*mc2*trc)*gDs-mc2*ms2*(9*trsgc+ms2*(3*trs+14*tr
     &   g+3*trc))+3*ms2**3*(trs+trg+trc))+18*cDs**2*(2*gDs+ms2)*(65536*
     &   trs*gDs**2+2*ms2*trg*gDs+ms2*(trsgc+ms2*trg)))-ms2*cDg**4*((2*g
     &   Ds+ms2)*(2359296*trs*gDs**3-2*ms2*(3*(trs+5*trg)+29*trc)*gDs**2
     &   +ms2*(3*ms2*trs+5*mc2*trs+3*ms2*trg+23*mc2*trg-26*ms2*trc)*gDs+
     &   9*mc2*ms2**2*trg)+3*ms2*trsgc*(2*(5*sck+6)*gDs**2+ms2*(5*sck+13
     &   )*gDs+3*ms2**2)-9*trsgc*cDs*(4*gDs**2-ms2**2)+cDs*(2*gDs+ms2)**
     &   2*(2*(8*trs-trg+10*trc)*gDs+9*ms2*trg))+cDg**2*gDs**3*((-2*gDs-
     &   ms2)*(1179648*(ms2*trg-mc2*trs)*gDs**2-2*ms2*(mc2*(6*trs+15*trg
     &   +29*trc)-3*ms2*trs)*gDs+ms2*(mc2**2*(1179648*trs+23*trg)-3*ms2*
     &   *2*(-trs+trg+trc)+mc2*ms2*(-3*trs-1179625*trg+3*trc)))-2*cDs*(2
     &   *gDs+ms2)*(3*ms2*(-2*trs*gDs+393216*trg*gDs+ms2*(2*trs+trg+trc)
     &   )-mc2*(1179644*trs*gDs+ms2*(14*trs+4*trg+13*trc)))-3*mc2*ms2*tr
     &   sgc*(2*(5*sck+3)*gDs+ms2*(5*sck+4))-1179648*(ms2*trg-mc2*trs)*c
     &   Ds**2*(2*gDs+ms2))-cDg*gDs**3*(2*gDs+ms2)*(2*cDs*(6*(mc-ms)*(ms
     &   +mc)*trs*gDs**2+(mc2*(9*trsgc+ms2*(6*trs-20*trg-29*trc))+5*mc2*
     &   *2*trs+3*ms2**2*(-trs+trg+trc))*gDs+9*mc2*ms2*(trsgc+mc2*trg))+
     &   6*cDs**2*(mc2*(2*trs*gDs+3*trsgc+3*ms2*trg)-2*ms2*trs*gDs)+mc2*
     &   gDs*(mc2*(ms2*(11*trs+13*(trg+trc))-4*trs*gDs)+6*ms2*trs*gDs-3*
     &   ms2**2*(trs+trg+trc)))+mc2*gDs**4*(2*gDs+ms2)*(6*(mc-ms)*(ms+mc
     &   )*trs*gDs**2+(mc2*(6*trs*cDs+9*trsgc+ms2*(6*trs-20*trg-29*trc))
     &   +3*ms2*(ms2*(-trs+trg+trc)-2*trs*cDs)+5*mc2**2*trs)*gDs+9*mc2*(
     &   trsgc*(cDs+ms2)+ms2*trg*(cDs+mc2)))+ms2**2*cDg**5*(2*gDs+ms2)*(
     &   (8*trs-trg+10*trc)*(2*gDs+ms2)-9*trsgc))/(ms*cDg**2*gDs**2*(2*g
     &   Ds+ms2)*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*cDg**2))/1.2e+1_dp+Ammm

      Ammm = 2*LsA*ms*(-mc2**2*trg*gDs**3-cDg**2*gDs*(2*(trg+trc)*cDs*g
     &   Ds-trsgc*cDs+mc2*ms2*trg)+mc2*cDg*gDs**2*((trg+2*trc)*gDs+trg*c
     &   Ds-trsgc)+ms2*cDg**3*((trg+2*trc)*gDs-trsgc))/(cDg**2*gDs**2)+
     &   Ammm

      Ammm = lc*(cDg**2*(-131072*(ms2*trg-mc2*trs)*gDs**2-131072*(ms2*t
     &   rg-mc2*trs)*cDs*gDs+mc2*ms2*(trsgc+ms2*trg))+mc2**2*(trsgc+ms2*
     &   trg)*gDs**2+131072*trs*cDg**3*gDs*(gDs+cDs)-2*mc2*(trsgc+ms2*tr
     &   g)*cDg*cDs*gDs)/(ms*(2*cDg+mc2)*(mc2*gDs**2-2*cDg*cDs*gDs+ms2*c
     &   Dg**2))/2.0_dp+Ammm

      Ammm = BfunX*(cDg**2*gDs*(4*trs*gDs**2-2*ms2*(trs+trg+2*trc)*gDs+
     &   ms2*(2*(trs+trg)*cDs+ms2*(trs+trg-trc)))+cDg*gDs**2*(trs*(4*gDs
     &   -ms2)*(2*(gDs+cDs)+ms2)-ms2*(trg+trc)*(4*gDs+2*cDs+ms2))+gDs**3
     &   *(2*(gDs+cDs)+ms2)*(2*trs*(gDs+cDs)-ms2*(-trs+trg+trc))-ms2*cDg
     &   **3*(2*trc*gDs+ms2*(trs+trg+trc)))/(ms*cDg**2*gDs**2)/4.0_dp+A
     &   mmm

      Ammm = ms*tr3c00fs*(-8*trc*gDs**2+cDg*(trsgc-(3*trs+2*(trg+trc))*
     &   (2*gDs+ms2))+(2*trg*(cDs+mc2)+3*trsgc-5*ms2*trc)*gDs+trsgc*(cDs
     &   +ms2)+ms2*trg*(cDs+mc2))/gDs**2+Ammm

      Ammm = epinv**2*ms*(cDg*(trsgc-2*(trg+trc)*gDs)+mc2*trg*gDs)/(cDg
     &   *gDs)/2.0_dp+Ammm

      Ammm = Ammm-3*ms*tr3c002fs*(trc*cDg*(2*gDs+ms2)+mc2*(trs+trg)*gD
     &   s)/gDs**2

      Ammm = ls*ms*(ms2*trsgc/(2*gDs+ms2)**2+trg)/2.0_dp+Ammm

      Ammm = 3*mc2*tr3s002ft*(trs*(2*cDg+mc2)*gDs+ms2*(trg+trc)*cDg)/(m
     &   s*cDg**2)+Ammm

      Ammm = 3*mc2*ms*tr3s001ft*((trg+trc)*gDs+trs*cDg)/cDg**2+Ammm

      Ammm = Ammm-3*ms*tr3c001fs*(2*gDs+ms2)*(trc*gDs+(trs+trg)*cDg)/g
     &   Ds**2

      Ammm=Ammm/zmmm

      return
      end
