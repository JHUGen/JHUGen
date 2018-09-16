      subroutine Bamp_ppp(q,mc,ms,Bppp)
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
      complex(dp):: trc,trg,trs,trsgc,zppp,Bppp

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
      zppp=za(2,4)*za(2,3)

      Bppp = mc*ms2*(cDg*gDs*(2*(-(trs+trg)*cDs+mc2*trg+ms2*trc)*gDs-2*
     &   trg*cDs**2+trsgc*(cDs-ms2)+ms2*trg*(cDs+mc2))-2*cDg**2*cDs*((tr
     &   s+trg)*gDs-ms2*trg)+gDs**2*(2*(mc2*trg+ms2*trc)*gDs+trsgc*(cDs-
     &   ms2)-ms2*trg*(cDs+mc2)))*tr5Xs/((cDs-mc*ms)*(cDs+mc*ms)*gDs**2)
     &   /2.0_dp

      Bppp = mc*mc2*(cDg*(cDg*(2*(ms2*trc-5*trs*cDs)*gDs+trg*(4*cDs**2+
     &   3*ms2*cDs-3*mc2*ms2))+gDs*(-6*trs*cDs*gDs+2*trg*cDs**2+3*ms2*tr
     &   g*cDs-mc2*ms2*trg)+2*cDg**2*(ms2*trc-2*trs*cDs))+trsgc*(gDs+cDg
     &   )*(2*cDs*gDs+cDg*(cDs+ms2)))*tr5Xc/(cDg**2*(cDs-mc*ms)*(cDs+mc*
     &   ms))/2.0_dp+Bppp

      Bppp = mc*ms2*(-2*(mc2*trg+ms2*trc)*gDs**2+(ms2*trg*(cDs+mc2)+2*(
     &   trs+trg)*cDg*cDs+trsgc*(ms2-cDs))*gDs-2*ms2*trg*cDg*cDs)*tr4Xs/
     &   ((cDs-mc*ms)*(cDs+mc*ms)*gDs)/2.0_dp+Bppp

      Bppp = mc*mc2*(-cDg*(cDs*(-6*trs*gDs+trsgc+3*ms2*trg)+4*trg*cDs**
     &   2+ms2*(trsgc-3*mc2*trg))-2*trsgc*cDs*gDs+cDg**2*(4*trs*cDs-2*ms
     &   2*trc))*tr4Xc/(cDg*(cDs**2-mc**2*ms**2))/2.0_dp+Bppp

      Bppp = Bppp-mc*(cDg*(mc2*ms2*(-6*trs*gDs+2*trg*gDs+trsgc+3*ms2*t
     &   rg)-2*trg*cDs**2*gDs+ms2*(trsgc+mc2*trg)*cDs)+2*mc2*ms2*trsgc*g
     &   Ds+2*ms2*cDg**2*(trc*cDs-2*mc2*trs))*tr3Xs/(cDg*(cDs**2-mc2*ms2
     &   ))/2.0_dp

      Bppp = mc*(mc2*ms2*(2*cDg*((trs+trg)*gDs-ms2*trg)+(ms2*trg-trsgc)
     &   *gDs)+cDs*gDs*(ms2*(trsgc+mc2*trg)-2*(mc2*trg+ms2*trc)*gDs))*tr
     &   3Xc/((cDs**2-mc2*ms2)*gDs)/2.0_dp+Bppp

      Bppp = (mc2*ms2*(trs+trg+trc)*(mc2*ms2-cDs**2)*gDs**3+cDg**2*gDs*
     &   (mc2*ms2*(4*trc*gDs**2-2*(trsgc+3*mc2*trs+2*mc2*trg)*gDs+ms2*(t
     &   rsgc+mc2*(3*sck-2)*trg))+cDs**2*(2*(2*mc2*(trs+trg)+ms2*trc)*gD
     &   s-ms2*(trsgc+mc2*(3*sck-2)*trg))+2*mc2*ms2*(-trs-trg+trc)*cDs*g
     &   Ds)-ms2*trc*cDg**3*(mc2*ms2-cDs**2)*(ms2-2*(gDs+cDs))+mc2*cDg*(
     &   mc2*ms2-cDs**2)*gDs**2*(2*trg*gDs+2*trsgc+ms2*trs)+2*ms2*trc*cD
     &   g**4*(mc2*ms2-cDs**2))*lVs/(mc*cDg**2*(cDs-mc*ms)*(cDs+mc*ms)*g
     &   Ds**2)/4.0_dp+Bppp

      Bppp = mc*(mc2*ms2*(trs+trg+trc)*(cDs-mc*ms)*(cDs+mc*ms)*gDs**3+c
     &   Dg*(cDs-mc*ms)*(cDs+mc*ms)*gDs**2*(2*trc*gDs**2-2*trg*cDs*gDs+m
     &   s2*(mc2*(trs+trg)-trsgc))+2*ms2*cDg**2*gDs**2*(cDs*(2*trc*gDs-t
     &   rsgc-mc2*trs+ms2*trc)-2*trg*cDs**2+mc2*ms2*(-trs+trg+trc))+ms2*
     &   trc*cDg**3*(mc2*ms2-cDs**2)*(ms2-2*(gDs+cDs))+2*ms2*trc*cDg**4*
     &   (cDs-mc*ms)*(cDs+mc*ms))*lVc/(ms2*cDg**2*(cDs-mc*ms)*(cDs+mc*ms
     &   )*gDs**2)/4.0_dp+Bppp

      Bppp = xs*cDs*(cDg*(ms2*trg-2*(trs+trg)*gDs)+trsgc*gDs)*lRcs/(ms*
     &   (xs**2-1)*cDg*gDs)+Bppp

      Bppp = (mc2*ms2*cDg**2*(-4*(-10*trs+28*trg+trc)*gDs**3-9*cDs*(4*(
     &   2*trg+trc)*gDs**2+2*ms2*trc*gDs-ms2*(trsgc+mc2*trg))+2*(-33*trs
     &   gc+3*ms2*trs+88*mc2*trs-26*ms2*trg+9*mc2*trg-8*ms2*trc+3*mc2*tr
     &   c)*gDs**2+mc2*ms2*(8*trs-15*sck*trg+23*trg-8*trc)*gDs+9*mc2*ms2
     &   *(trsgc+ms2*trg))+mc2*cDg*gDs**2*(4*(trg-8*trc)*gDs**3-2*cDs*(4
     &   *(4*trc-5*trg)*gDs**2+9*ms2*(trsgc+mc2*(3*trg+trc)))-2*(3*mc2-1
     &   6*ms2)*(trg+trc)*gDs**2+36*trg*cDs**2*gDs+2*ms2*(-9*trsgc+mc2*(
     &   17*trs-50*trg+4*trc)+3*ms2*(trs+trg+trc))*gDs+mc2*ms2*(-54*trsg
     &   c+ms2*(3*trs-44*trg-8*trc)+41*mc2*trs))+mc2**2*gDs**2*(2*(trg-8
     &   *trc)*gDs**3-cDs*(4*(4*trc-5*trg)*gDs**2+9*ms2*(trsgc+mc2*trg))
     &   +16*ms2*(trg+trc)*gDs**2+18*trg*cDs**2*gDs+ms2*(-9*trsgc+mc2*(7
     &   *trs-22*trg+5*trc)+3*ms2*(trs+trg+trc))*gDs-9*mc2*ms2*(trsgc+ms
     &   2*trg))+ms2*cDg**3*(4*mc2*(47*trs+9*trg+3*trc)*gDs**2-2*(18*ms2
     &   *trc*cDs+mc2*ms2*(-8*trs+15*sck*trg-23*trg+15*trc)+3*mc2**2*trc
     &   )*gDs+6*(3*ms2*trsgc+mc2*ms2*(3*trg+trc)-mc2**2*trc)*cDs+mc2*ms
     &   2*(27*trsgc+mc2*(-26*trs-17*trg+3*trc)+18*ms2*trg+5*ms2*trc))-2
     &   *ms2*cDg**4*(2*(7*ms2+3*mc2)*trc*gDs+6*(mc-ms)*(ms+mc)*trc*cDs-
     &   9*ms2*trsgc+26*mc2*ms2*trs+17*mc2*ms2*trg-5*ms2**2*trc-6*mc2*ms
     &   2*trc+3*mc2**2*trc)+12*ms2*(ms2-mc2)*trc*cDg**5)/(mc*ms2*cDg**2
     &   *(2*cDg+mc2)*gDs**2)/1.2e+1_dp+Bppp

      Bppp = Bppp-LsB1*mc*(cDs*(-2*(mc2*trg+ms2*trc)*gDs**3+ms2*(trsgc
     &   +mc2*trg)*gDs**2-ms2*cDg*(2*(trs+trg)*cDg-trsgc+ms2*trg)*gDs+2*
     &   ms2**2*trg*cDg**2)+2*cDs**2*gDs*(2*cDg*((trs+trg)*gDs-ms2*trg)+
     &   (ms2*trg-trsgc)*gDs)+ms2*gDs*(cDg*(-2*mc2*trs*gDs+2*ms2*trc*gDs
     &   -ms2*trsgc+mc2*ms2*trg)+mc2*(trsgc-ms2*trg)*gDs))/((cDs**2-mc2*
     &   ms2)*gDs**2)

      Bppp = B0cgsf*(cDg**2*gDs*(-2*mc2*trc*gDs**3-2*(-mc2*trg*cDs+mc2*
     &   ms2*(trg-2*trs)+ms2**2*trc)*gDs**2+ms2*(-2*(ms2+mc2)*trc*cDs+ms
     &   2*trsgc+mc2*(ms2*(2*trs+trg-2*trc)-trsgc)+mc2**2*(4*trs+trg))*g
     &   Ds+ms2**2*(mc2*trg*(cDs+ms2)+trsgc*(cDs+mc2)))-mc2*gDs**3*(2*tr
     &   c*gDs**3+cDs*(2*(trc-trg)*gDs**2+ms2*(trsgc+mc2*trg))-2*ms2*(tr
     &   g+trc)*gDs**2-2*trg*cDs**2*gDs-ms2*(-trsgc+ms2*(trs+trg+trc)-3*
     &   mc2*trg)*gDs+mc2*ms2*(trsgc+ms2*trg))+mc2*cDg*gDs**2*(-4*trc*gD
     &   s**3-cDs*(2*(trc-2*trg)*gDs**2+2*ms2*(2*trg+trc)*gDs+ms2*(trsgc
     &   +mc2*trg))+2*ms2*(trs-trg+trc)*gDs**2+2*trg*cDs**2*gDs+2*ms2*(-
     &   trsgc-(ms2+2*mc2)*(trg-trs))*gDs-mc2*ms2*(trsgc+ms2*trg))+ms2*c
     &   Dg**3*(-2*(mc2*(-trs-trg+trc)+ms2*trc)*gDs**2+cDs*(ms2*(trsgc+m
     &   c2*trg)-2*mc2*trc*gDs)+2*ms2*(trsgc-mc2*trs)*gDs+mc2*ms2*(trsgc
     &   +ms2*trg))+ms2*cDg**4*(2*(ms2-2*mc2)*trc*gDs+2*(ms**2-mc**2)*tr
     &   c*cDs+ms2*(trsgc+mc2*(-3*trs-2*trg+trc)))+2*ms2*(ms2-mc2)*trc*c
     &   Dg**5)/(mc*ms2*cDg**2*gDs**2*(gDs+cDg))/4.0_dp+Bppp

      Bppp = mc*tr3s00ft*(-2*trg*gDs**3-2*(2*trg*cDs+3*ms2*(trg+trc))*g
     &   Ds**2+ms2*cDg*(2*(trs+4*trg+trc)*gDs+2*(2*trg+trc)*cDs+3*trsgc-
     &   4*mc2*trs+5*ms2*trg+3*ms2*trc)+(ms2*(trsgc+mc2*(trs+3*trg))-2*t
     &   rg*cDs**2)*gDs+ms2*(mc2*trg*(cDs+ms2)+trsgc*(cDs+mc2))-6*ms2*tr
     &   s*cDg**2)/(ms2*cDg**2)+Bppp

      Bppp = BfunX*(-mc2*cDg*gDs**2*((trg+trc)*(2*gDs**2-ms2**2)+ms2*tr
     &   s*(2*(gDs+cDs)+ms2))-ms2*cDg**3*(4*trc*gDs**2+2*trc*cDs*(4*gDs+
     &   mc2)+2*mc2*trc*gDs+4*trc*cDs**2+mc2*ms2*(trs+trg-trc))-mc2*gDs*
     &   *3*((trg+trc)*(2*ms2*(2*gDs+cDs)+2*gDs*(gDs+cDs)+ms2**2)+ms2*tr
     &   s*(2*(gDs+cDs)+ms2))-2*ms2*trc*cDg**4*(4*(gDs+cDs)+mc2)-mc2*ms2
     &   *cDg**2*gDs*(ms2*(trs+trg+trc)-2*(trg+trc)*gDs)-4*ms2*trc*cDg**
     &   5)/(mc*ms2*cDg**2*gDs**2)/4.0_dp+Bppp

      Bppp = Bppp-lc*mc*(mc2*(-2*trg*gDs**2-2*trg*cDs*gDs+ms2*trsgc-3*
     &   mc2*ms2*trg)+2*cDg*(2*trc*gDs**2-2*trg*cDs*gDs+ms2*(2*trsgc-5*m
     &   c2*trg+mc2*trc))+4*ms2*(trc-2*trg)*cDg**2)/(ms2*(2*cDg+mc2)**2)
     &   /2.0_dp

      Bppp = B0csf*mc*(cDs*(-2*trc*gDs**2+(trsgc+ms2*(trs+trg-2*trc)+mc
     &   2*trs)*gDs+2*mc2*ms2*trg)+ms2*(-2*trc*gDs**2+(trsgc+mc2*trg-(ms
     &   2+mc2)*trc)*gDs+mc2*(ms2+mc2)*trg)+cDg*(ms2*(-trc*(2*gDs+ms2+mc
     &   2)+trsgc+mc2*trg)+cDs*(-2*trc*gDs+trsgc+ms2*(trs+trg-2*trc)+mc2
     &   *trs)+2*trs*cDs**2)-cDs**2*((ms2+mc2)*trg-2*trs*gDs)-2*trg*cDs*
     &   *3)/((cDs-mc*ms)*(cDs+mc*ms)*(gDs+cDg))/2.0_dp+Bppp

      Bppp = LsB2*mc*(-2*mc2*trsgc*cDs*gDs**2+cDg**2*(4*cDs**2*(-(2*trs
     &   +trc)*gDs+trsgc+ms2*trg)+mc2*ms2*(2*(trs+trc)*gDs-3*trsgc-ms2*t
     &   rg)+cDs*(4*mc2*trs*gDs+ms2*trsgc-7*mc2*ms2*trg)+8*trg*cDs**3)-c
     &   Dg*gDs*(mc2*cDs*(-6*trs*gDs+trsgc+3*ms2*trg)+2*(mc2*trg-trsgc)*
     &   cDs**2+mc2*ms2*(trsgc-mc2*trg))+cDg**3*(-8*trs*cDs**2+2*ms2*trc
     &   *cDs+4*mc2*ms2*trs))/(cDg**2*(cDs**2-mc2*ms2))+Bppp

      Bppp = 3*mc*tr3s001ft*(-2*(trg+trc)*gDs**2+trs*(2*cDg+mc2)*gDs+ms
     &   2*(trg+trc)*cDg)/cDg**2+Bppp

      Bppp = 3*mc*tr3s002ft*((trg+trc)*gDs*(ms2*(2*cDg+mc2)-2*gDs*(gDs+
     &   cDs))+ms2*trs*cDg*(2*cDg+mc2))/(ms2*cDg**2)+Bppp

      Bppp = epinv*(cDg*(2*(trs+trg)*(mc*ms*(xs**2-1)-4*lp*xs*cDs)*gDs+
     &   ms2*trg*(4*lp*xs*cDs-mc*ms*(2*sck-1)*(xs**2-1)))+trsgc*(4*lp*xs
     &   *cDs+mc*(ms-ms*xs**2))*gDs)/(ms*(xs**2-1)*cDg*gDs)/4.0_dp+Bppp

      Bppp = mc*tr2fu*cDs*(cDg*(ms2*trg-2*(trs+trg)*gDs)+trsgc*gDs)/(cD
     &   g*gDs)+Bppp

      Bppp = Bppp-ms2*tr3c00fs*(mc2*(3*trs*gDs+4*trg*gDs+trsgc+ms2*trg
     &   )+cDg*(2*trc*gDs+trsgc-mc2*(2*trs+trg)+3*ms2*trc)+cDs*(-2*trc*g
     &   Ds+trsgc+mc2*trg))/(mc*gDs**2)

      Bppp = Bppp-3*ms2*tr3c001fs*(trc*cDg*(2*gDs+ms2)+mc2*(trs+trg)*g
     &   Ds)/(mc*gDs**2)

      Bppp = Bppp-3*mc*ms2*tr3c002fs*(trc*gDs+(trs+trg)*cDg)/gDs**2

      Bppp = ls*ms2*(-2*trc*gDs+trsgc+mc2*trg)/(mc*(2*gDs+ms2))/2.0_dp+
     &   Bppp

      Bppp=Bppp/zppp

      return
      end
