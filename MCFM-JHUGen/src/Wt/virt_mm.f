      function virt_mm(mQ,ig,is,ie,in,ic,p)
      implicit none
      include 'types.f'
      complex(dp):: virt_mm

************************************************************************
*     Author: F. Tramontano                                            *
*     December, 2008                                                   *
************************************************************************
c---- One-loop -- helicity amplitude for W+c production
c---- hQ=-1  with Q(mu)=Q1(mu)+mQ^2/2/Q.g*g(mu)
c---- hg=-1  with s(mu) as gauge vector
c for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4)) + Qbar(p5)
c For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ Q(p5)
c----
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'b0.f'
      include 'scheme.f'
      include 'scale.f'
      include 'zprods_decl.f'
      integer:: is,ig,ic,ie,in,nu,jpart
      real(dp):: p(mxpart,4),q(mxpart,4)
      real(dp):: dot,taucg,taugs,taucs,msq,mQ
      real(dp):: qsq,tcg,tcs,qsqhat,qgs,qcg,qcs,cgs1
      real(dp):: ddilog,epin,epin2,xlog
      complex(dp):: smm,lomm
      complex(dp):: amm1,amm2,amm3,amm4,amm5
      complex(dp):: fq2m2,fc002,fc00,ffc00,ffc002,fI3me
      complex(dp):: ffcg1,ffcg3,ffcs1,ffcs2,fun1,fun2,fun4
      complex(dp):: fL6m1,fL6m2,fL6m3,L6m1,L6m2,L6m3
      complex(dp):: lnrat,C0fa2m,C0fb2m,I3me
      scheme='dred'
      msq=mQ**2
      xlog=log(musq/msq)
      epin=epinv+xlog
      epin2=epinv**2+epinv*xlog+half*xlog**2
      taugs=+two*dot(p,is,ig)
      taucs=+two*dot(p,ic,is)
      taucg=+two*dot(p,ic,ig)
      qsqhat=taucs+taucg+taugs
      qsq=qsqhat+msq
      tcg=taucg+msq
      tcs=taucs+msq
      qgs=qsqhat-taugs
      qcg=qsqhat-taucg
      qcs=qsqhat-taucs
      cgs1=taucs-msq*taugs/taucg
      fq2m2=-msq/qsq+lnrat(-qsqhat,msq)*msq*qsqhat/qsq**2
      ffcs1=fun4(tcs,qsq,msq)
      ffcg1=fun1(tcg,qsq,msq)
      ffcg3=fun2(qsq,tcg,msq)
      ffcs2=fun2(tcs,qsq,msq)
      fL6m1=L6m1(taucg,taucs,taugs,msq,qsq)
      fL6m2=L6m2(taucg,taucs,taugs,msq,qsq)
      fL6m3=L6m3(taucg,taucs,taugs,msq,qsq)
      fI3me=I3me(msq,taugs,qsq)
      ffc00=fc00(msq,taugs,qsq)
      ffc002=fc002(msq,taugs,qsq)
      do nu=1,4
      do jpart=1,5
      q(jpart,nu)=p(jpart,nu)
      if (jpart==ic) q(ic,nu)=p(ic,nu)-msq/taucg*p(ig,nu)
      enddo
      enddo
      call spinoru(5,q,za,zb)
      amm1=za(ig,ie)*zb(ig,in)/za(is,ic)/zb(ig,ic)**2
      amm2=za(ig,ie)*zb(is,in)/zb(ig,ic)
      amm3=za(is,ie)*zb(is,in)/za(is,ic)/zb(ig,ic)**2
      amm4=za(ie,ic)*zb(is,in)*zb(is,ic)/zb(ig,is)/zb(ig,ic)
      amm5=-amm1-amm3
      smm=  - 1.D0/2.D0*taucs*taucg**(-2)*mQ**2*xn**(-1)*qcs*qcg*
     &    cgs1**(-1)*amm1*fL6m1 - 1.D0/2.D0*taucs*taucg**(-2)*mQ**2*
     &    xn**(-1)*qcs*cgs1**(-1)*amm2*fL6m1 - 1.D0/2.D0*taucs*
     &    taucg**(-2)*mQ**2*xn**(-1)*qcs**2*cgs1**(-1)*amm3*fL6m1 - 1.D0
     &    /2.D0*taucs*taucg**(-1)*taugs**(-1)*xn**(-1)*qcs*amm2*fL6m1
     &     - 1.D0/2.D0*taucs*taucg**(-1)*taugs**(-1)*xn**(-1)*qcs**2*
     &    amm3*fL6m1 - 1.D0/2.D0*taucs*taucg**(-1)*qsqhat*xn**(-1)*amm1
     &    *fL6m1 + taucs*taucg*taugs**(-1)*cf*amm3 + 1.D0/2.D0*taucs*
     &    taucg*taugs**(-1)*xn*amm3*fL6m3 + 1.D0/2.D0*taucs*taugs**(-1)
     &    *mQ**2*xn**(-1)*amm1*fL6m2 + taucs*taugs**(-1)*cf*qcs*amm5*
     &    fq2m2 + 1.D0/2.D0*taucs*taugs**(-1)*cf*qcg*amm5*fq2m2 - 4.D0*
     &    taucs*qsq**(-1)*mQ**2*cf*amm1 - 2.D0*taucs*qsq**(-1)*mQ**2*
     &    xn**(-1)*amm1 + 2.D0*taucs*qsq**(-1)*mQ**2*xn*amm1 + 4.D0*
     &    taucs*qsq**(-1)*cf*cgs1*amm1 + 2.D0*taucs*qsq**(-1)*xn**(-1)*
     &    cgs1*amm1 - 2.D0*taucs*qsq**(-1)*xn*cgs1*amm1 - 4.D0*taucs*cf
     &    *amm1*fq2m2
      smm = smm + 3.D0*taucs*cf*amm1 + taucs*cf*amm3 - 2.D0*taucs*
     &    xn**(-1)*amm1*ffcs2 - 2.D0*taucs*xn**(-1)*amm1*fq2m2 + 2.D0*
     &    taucs*xn*amm1*fq2m2 - taucs*xn*amm1 - 1.D0/2.D0*taucs**2*
     &    taugs**(-1)*cf*amm1 - 1.D0/2.D0*taucs**2*taugs**(-1)*cf*amm3
     &     + 2.D0*taucg**(-1)*taugs*mQ**2*cf*amm3*ffcg3 - 2.D0*
     &    taucg**(-1)*mQ**2*cf*amm2*ffcg3 - taucg**(-1)*mQ**2*xn*amm2*
     &    fL6m3 - 1.D0/2.D0*taucg*taugs**(-1)*tcs*xn**(-1)*amm3*fL6m2
     &     - 1.D0/2.D0*taucg*taugs**(-1)*xn**(-1)*amm2*fL6m2 - taucg*
     &    qsq**(-1)*mQ**2*cf*amm3*ffcg3 + 1.D0/2.D0*taugs**(-1)*qsq*
     &    xn**(-1)*amm2*fL6m2 + 1.D0/2.D0*taugs**(-1)*mQ**2*xn**(-1)*
     &    qcs*amm3*fL6m2 + taugs**(-1)*cf*cgs1*amm2 - taugs**(-1)*
     &    xn**(-1)*cgs1*amm2*fL6m2 + 1.D0/2.D0*taugs**(-1)*xn*qcg*amm2*
     &    fL6m3 + taugs*tcs*xn**(-1)*amm1*ffcs1 + taugs*qsq**(-1)*mQ**2
     &    *cf*amm3*ffcg1 - 2.D0*taugs*qsq**(-1)*mQ**2*cf*amm3*ffcg3 +
     &    taugs*xn**(-1)*amm1*ffcs2 + tcg**(-1)*mQ**4*cf*amm3 -
     &    qsq**(-1)*mQ**2*cf*amm2*ffcg3
      smm = smm - mQ**2*cf*amm1 - 2.D0*mQ**2*cf*amm3*ffcg3 - 3.D0*mQ**2
     &    *cf*amm3 - mQ**2*xn**(-1)*amm1*ffcs2 - mQ**2*xn**(-1)*amm3*
     &    ffcg3 - 1.D0/2.D0*mQ**2*xn*qcs*cgs1**(-1)*amm3*fL6m3 - 1.D0/2.
     &    D0*mQ**2*xn*qcg*cgs1**(-1)*amm1*fL6m3 - 1.D0/2.D0*mQ**2*xn*
     &    cgs1**(-1)*amm2*fL6m3 + 1.D0/2.D0*mQ**2*xn*amm1*fL6m3 - mQ**2
     &    *xn*amm3*fL6m3 - xlog*b0*amm2 - xlog*b0*amm4 + 3.D0/2.D0*epin
     &    *cf*amm2 + 3.D0/2.D0*epin*cf*amm4 + epin*b0*amm2 + epin*b0*
     &    amm4 - 1.D0/2.D0*epin*xn**(-1)*amm2 - 1.D0/2.D0*epin*xn**(-1)
     &    *amm4 + 1.D0/2.D0*epin*xn*amm2 + 1.D0/2.D0*epin*xn*amm4 - 1.D0
     &    /2.D0*epin2*xn**(-1)*amm2 - 1.D0/2.D0*epin2*xn**(-1)*amm4 + 3.
     &    D0/2.D0*epin2*xn*amm2 + 3.D0/2.D0*epin2*xn*amm4 + 1.D0/2.D0*
     &    cf*qcs*amm5*fq2m2 + 5.D0/2.D0*cf*amm2 + 5.D0/2.D0*cf*amm4 - 1.
     &    D0/2.D0*xn**(-1)*pisqo6*amm2 - 1.D0/2.D0*xn**(-1)*pisqo6*amm4
     &     - xn**(-1)*amm2*ffcs2 - xn**(-1)*amm2*fL6m2 - xn**(-1)*amm2
     &     + 1.D0/2.D0*xn*pisqo6*amm2 + 1.D0/2.D0*xn*pisqo6*amm4 -
     &    ddilog(
     & msq**(-1)*tcs)*xn**(-1)*amm2 - ddilog(msq**(-1)*tcs)*xn**(-1)*
     &    amm4 + ddilog(msq**(-1)*tcg)*xn*amm2 + ddilog(msq**(-1)*tcg)*
     &    xn*amm4 - lnrat( - taucs,msq)*taucs*tcs**(-1)*mQ**2*xn**(-1)*
     &    amm1 - 2.D0*lnrat( - taucs,msq)*taucs**2*tcs**(-1)*xn**(-1)*
     &    amm1 + lnrat( - taucs,msq)*tcs**(-1)*mQ**2*xn**(-1)*amm2 -
     &    lnrat( - taucs,msq)*tcs*xn**(-1)*amm3 + lnrat( - taucs,msq)*
     &    epin*xn**(-1)*amm2 + lnrat( - taucs,msq)*epin*xn**(-1)*amm4
     &     - lnrat( - taucs,msq)**2*xn**(-1)*amm2 - lnrat( - taucs,msq)
     &    **2*xn**(-1)*amm4 - 2.D0*lnrat( - taucg,msq)*taucg*tcg**(-1)*
     &    mQ**2*cf*amm1 - lnrat( - taucg,msq)*taucg*tcg**(-1)*mQ**2*
     &    xn**(-1)*amm1 + lnrat( - taucg,msq)*taucg**2*tcg**(-2)*mQ**2*
     &    cf*amm3 - 2.D0*lnrat( - taucg,msq)*tcg**(-1)*mQ**2*cf*qcs*
     &    amm3 - 2.D0*lnrat( - taucg,msq)*tcg**(-1)*mQ**4*cf*amm1 -
     &    lnrat( - taucg,msq)*tcg**(-1)*mQ**4*xn**(-1)*amm1 - lnrat( -
     &    taucg,msq)*epin*xn*amm2 - lnrat( - taucg,msq)*epin*xn*amm4 +
     &    lnrat(
     &  - taucg,msq)**2*xn*amm2 + lnrat( - taucg,msq)**2*xn*amm4 -
     &    lnrat( - taugs,msq)*taucs*taugs**(-1)*cf*amm2 + lnrat( -
     &    taugs,msq)*taucs*cf*amm5 + lnrat( - taugs,msq)*taucs*xn**(-1)
     &    *amm3 + lnrat( - taugs,msq)*taucs**2*taugs**(-1)*cf*amm3 + 1.D
     &    0/2.D0*lnrat( - taugs,msq)*taucs**2*mQ**(-2)*cf*amm5 + 1.D0/2.
     &    D0*lnrat( - taugs,msq)*taugs*cf*amm5 - lnrat( - taugs,msq)*
     &    mQ**2*cf*amm3 - lnrat( - taugs,msq)*mQ**2*cf*amm5 - lnrat( -
     &    taugs,msq)*epin*xn*amm2 - lnrat( - taugs,msq)*epin*xn*amm4 +
     &    1.D0/2.D0*lnrat( - taugs,msq)*xn**(-1)*qcg*amm5 - lnrat( -
     &    taugs,msq)*xn*amm2 + 1.D0/2.D0*lnrat( - taugs,msq)**2*xn*amm2
     &     + 1.D0/2.D0*lnrat( - taugs,msq)**2*xn*amm4 - 3.D0*lnrat( -
     &    qsqhat,msq)*taucs*taucg*taugs**(-1)*qsq**(-1)*qsqhat*cf*amm3
     &     - 2.D0*lnrat( - qsqhat,msq)*taucs*taugs**(-1)*qsq**(-1)*
     &    mQ**2*qsqhat*cf*amm3 + 3.D0*lnrat( - qsqhat,msq)*taucs*
     &    qsq**(-1)*qsqhat*cf*amm3 + 2.D0*lnrat( - qsqhat,msq)*taucs*
     &    qsq**(-1)*qsqhat*cf*amm5
      smm = smm + 3.D0/2.D0*lnrat( - qsqhat,msq)*taucs*qsq**(-1)*qsqhat
     & *xn*amm1 - lnrat( - qsqhat,msq)*taucs**2*taugs**(-1)*qsq**(-1)*
     &    qsqhat*cf*amm3 - 1.D0/2.D0*lnrat( - qsqhat,msq)*taucs**2*
     &    taugs**(-1)*qsq**(-1)*qsqhat*cf*amm5 - 1.D0/2.D0*lnrat( -
     &    qsqhat,msq)*taucs**2*qsq**(-1)*mQ**(-2)*qsqhat*cf*amm5 +
     &    lnrat( - qsqhat,msq)*taucg**(-1)*taugs*qsq**(-1)*mQ**2*qsqhat
     &    *cf*amm3 - lnrat( - qsqhat,msq)*taugs**(-1)*qsq**(-1)*qsqhat*
     &    cf*qcg*amm2 - lnrat( - qsqhat,msq)*taugs**(-1)*qsq**(-1)*
     &    qsqhat*cf*cgs1*amm2 - lnrat( - qsqhat,msq)*taugs*qsq**(-1)*
     &    qsqhat*cf*amm3 - 3.D0/2.D0*lnrat( - qsqhat,msq)*taugs*
     &    qsq**(-1)*qsqhat*cf*amm5 - 1.D0/2.D0*lnrat( - qsqhat,msq)*
     &    taugs*qsq**(-1)*qsqhat*xn*amm1 + 2.D0*lnrat( - qsqhat,msq)*
     &    qsq**(-1)*mQ**2*qsqhat*cf*amm3 - 3.D0*lnrat( - qsqhat,msq)*
     &    qsq**(-1)*mQ**2*qsqhat*cf*amm5 + 3.D0/2.D0*lnrat( - qsqhat,
     &    msq)*qsq**(-1)*mQ**2*qsqhat*xn**(-1)*amm1 - 1.D0/2.D0*lnrat(
     &     - qsqhat,msq)*qsq**(-1)*mQ**2*qsqhat*xn*amm1
      smm = smm + lnrat( - qsqhat,msq)*qsq**(-1)*qsqhat*xn**(-1)*qcg*
     & amm1 + 1.D0/2.D0*lnrat( - qsqhat,msq)*qsq**(-1)*qsqhat*xn**(-1)*
     &    qcg*amm3 + lnrat( - qsqhat,msq)*qsq**(-1)*qsqhat*xn**(-1)*
     &    amm2 - 2.D0*C0fb2m(tcg,msq)*taucs**(-1)*taucg**(-1)*taugs*
     &    mQ**4*xn**(-1)*amm2 + C0fb2m(tcg,msq)*taucs**(-1)*taucg*mQ**2
     &    *xn**(-1)*cgs1*amm3 + C0fb2m(tcg,msq)*taucs**(-1)*mQ**2*
     &    xn**(-1)*qcg*amm2 - C0fb2m(tcg,msq)*taucs**(-1)*mQ**4*
     &    xn**(-1)*amm2 + 2.D0*C0fa2m(tcs,qsq,msq)*taucs**(-1)*
     &    taucg**(-2)*mQ**4*xn**(-1)*qcs**2*amm2 - C0fa2m(tcs,qsq,msq)*
     &    taucs**(-1)*taucg**(-1)*mQ**2*xn**(-1)*qcs**2*amm2 - C0fa2m(
     &    tcs,qsq,msq)*taucs**(-1)*taucg**(-1)*mQ**4*xn**(-1)*qcs*amm2
     &     - C0fa2m(tcs,qsq,msq)*taucs**(-1)*mQ**2*xn**(-1)*qcs*cgs1*
     &    amm3 + C0fa2m(tcs,qsq,msq)*taucs**(-1)*mQ**2*xn**(-1)*qcs*
     &    amm2 - 2.D0*C0fa2m(tcs,qsq,msq)*taucs*mQ**2*xn**(-1)*amm1 +
     &    C0fa2m(tcs,qsq,msq)*taucg**(-1)*taugs*mQ**4*xn**(-1)*amm1 - 2.
     &    D0*C0fa2m(tcs,qsq,msq)*taucg**(-1)*mQ**2*xn**(-1)*qcs*amm2
      smm = smm - C0fa2m(tcs,qsq,msq)*mQ**2*xn**(-1)*amm2 - 18.D0*
     &    ffc002*taucs*taugs**(-1)*cf*qcs*amm5 + 6.D0*ffc002*taucs**2*
     &    taugs**(-1)*mQ**(-2)*cf*qgs*amm5 + 6.D0*ffc002*taucg*
     &    taugs**(-1)*cf*qcg*amm5 + 12.D0*ffc002*mQ**2*cf*amm5 + 24.D0*
     &    ffc00*taucs*taucg*taugs**(-1)*cf*amm3 + 4.D0*ffc00*taucs*
     &    taugs**(-2)*cf*qgs*amm2 + 4.D0*ffc00*taucs*taugs**(-1)*tcs*cf
     &    *amm1 + 16.D0*ffc00*taucs*taugs**(-1)*mQ**2*cf*amm3 + 12.D0*
     &    ffc00*taucs*taugs**(-1)*cf*qcs*amm1 + 8.D0*ffc00*taucs*
     &    taugs**(-1)*cf*amm2 + 12.D0*ffc00*taucs*cf*amm3 - 4.D0*ffc00*
     &    taucs**2*taugs**(-2)*qsqhat*cf*amm3 + 8.D0*ffc00*taucs**2*
     &    taugs**(-1)*cf*amm3 + 12.D0*ffc00*taucg*taugs**(-1)*mQ**2*cf*
     &    amm1 - 8.D0*ffc00*taucg*taugs**(-1)*mQ**2*cf*amm3 - 4.D0*
     &    ffc00*taucg*taugs**(-1)*cf*amm2 - 16.D0*ffc00*taugs**(-1)*
     &    mQ**2*cf*qcs*amm1 - 8.D0*ffc00*taugs**(-1)*mQ**2*cf*amm2 + 2.D
     &    0*ffc00*taugs**(-1)*xn**(-1)*qgs*qcg*amm5 - 4.D0*ffc00*tcg*
     &    xn**(-1)*amm5
      smm = smm - 28.D0*ffc00*mQ**2*cf*amm3 + fI3me*taucs*mQ**2*
     &    xn**(-1)*amm3 - 3.D0*fI3me*taucs*mQ**2*xn*amm3 - fI3me*taucs*
     &    xn*amm2 - fI3me*taucg**(-1)*taugs*mQ**4*xn*amm1 + 2.D0*fI3me*
     &    taucg**(-1)*taugs*mQ**4*xn*amm3 - 2.D0*fI3me*taucg**(-1)*
     &    mQ**2*xn*cgs1*amm2 - fI3me*taucg*mQ**2*xn**(-1)*amm3 - fI3me*
     &    taugs*mQ**2*xn**(-1)*amm1 - 2.D0*fI3me*taugs*mQ**2*xn**(-1)*
     &    amm3 + 2.D0*fI3me*mQ**2*cf*amm2 - fI3me*mQ**2*xn*amm2 + fI3me
     &    *xn*cgs1*amm2


      lomm =-(za(ig,ie)*zb(ig,is)+za(ie,ic)*zb(is,ic))*zb(is,in)
     & /zb(ig,is)/zb(ig,ic)

c--- include finite counterterm to go from DR to MSbar scheme
c--- alphas(DR) = alphas(MSbar) * (1+ (Nc / 6) * alphas(MSbar) / (2*pi))
      smm=smm + lomm * xn/six

c--- include finite counterterm due to FDH scheme
c--- gw = gw * ( 1 - 2 * cf * alphas(MSbar) / (2*pi))
      smm=smm - lomm * cf * two

      virt_mm=smm
      return
      end
