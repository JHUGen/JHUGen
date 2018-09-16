      function virt_pm(mQ,ig,is,ie,in,ic,p)
      implicit none
      include 'types.f'
      complex(dp):: virt_pm

************************************************************************
*     Author: F. Tramontano                                            *
*     December, 2008                                                   *
************************************************************************
c---- One-loop +- helicity amplitude for W+c production
c---- hQ=+1  with Q(mu)=Q1(mu)+mQ^2/2/Q.g*g(mu)
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
      complex(dp):: spm,lopm
      complex(dp):: apm1,apm2,apm3,apm4,apm5
      complex(dp):: fq2m2,fc002,fc00,ffc00,ffc002,fI3me
      complex(dp):: ffcg1,ffcg3,fun1,fun2
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
      ffcg1=fun1(tcg,qsq,msq)
      ffcg3=fun2(qsq,tcg,msq)
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
      apm1=za(ig,ie)*zb(ig,in)/za(is,ic)**2/zb(ig,is)/zb(ig,ic)
      apm2=za(ig,ie)*zb(is,in)*zb(is,ic)/zb(ig,is)
      apm3=za(ig,ie)*zb(is,in)/za(is,ic)/zb(ig,is)
      apm4=za(is,ie)*zb(is,in)/za(is,ic)**2/zb(ig,is)/zb(ig,ic)
      apm5=-apm1-apm4
      spm=  - 3.D0*taucs*taucg**(-1)*taugs*mQ*cf*apm5 - 5.D0/2.D0*taucs
     &    *taucg**(-1)*taugs**2*mQ**(-1)*cf*apm5*fq2m2 - taucs*taucg*
     &    qsq*mQ**(-3)*cf*apm1*fq2m2 - taucs*taucg*qsq*mQ**(-3)*cf*apm5
     &    *fq2m2 - 3.D0/2.D0*taucs*taucg*mQ**(-1)*cf*apm5*fq2m2 - taucs
     &    *taucg*mQ**(-1)*cf*apm5 + taucs*taugs**(-1)*mQ*xn**(-1)*qcs*
     &    apm4*fL6m2 + 1.D0/2.D0*taucs*taugs**(-1)*mQ*xn**(-1)*cgs1*
     &    apm1*fL6m2 + 1.D0/2.D0*taucs*taugs**(-1)*mQ*xn**(-1)*apm3*
     &    fL6m2 - 4.D0*taucs*taugs*mQ**(-1)*cf*apm5*fq2m2 - 3.D0/2.D0*
     &    taucs*mQ*cf*apm5 + 1.D0/2.D0*taucs**2*taucg**(-2)*taugs*mQ*
     &    xn**(-1)*qcs*cgs1**(-1)*apm4*fL6m1 + 1.D0/2.D0*taucs**2*
     &    taucg**(-2)*taugs*mQ*xn**(-1)*qcg*cgs1**(-1)*apm1*fL6m1 + 1.D0
     &    /2.D0*taucs**2*taucg**(-2)*taugs*mQ*xn**(-1)*cgs1**(-1)*apm3*
     &    fL6m1 - 2.D0*taucs**2*taucg**(-1)*taugs*mQ**(-1)*cf*apm5*
     &    fq2m2 + 3.D0/2.D0*taucs**2*taucg**(-1)*mQ*cf*apm5 - 1.D0/2.D0
     &    *taucs**2*taucg**(-1)*mQ*xn**(-1)*apm1*fL6m1 - 1.D0/2.D0*
     &    taucs**2*taugs**(-1)*mQ*xn**(-1)*qcs*cgs1**(-1)*apm4*fL6m2
      spm = spm - 1.D0/2.D0*taucs**2*taugs**(-1)*mQ*xn**(-1)*qcg*
     & cgs1**(-1)*apm1*fL6m2 - 1.D0/2.D0*taucs**2*taugs**(-1)*mQ*
     &    xn**(-1)*cgs1**(-1)*apm3*fL6m2 - 3.D0/2.D0*taucs**2*mQ**(-1)*
     &    cf*apm5*fq2m2 + 1.D0/2.D0*taucs**3*taucg**(-1)*taugs**(-1)*
     &    qsq*mQ**(-1)*cf*apm5*fq2m2 - 1.D0/2.D0*taucs**3*taucg**(-1)*
     &    taugs**(-1)*mQ*cf*apm5*fq2m2 + taucg**(-2)*taugs*tcg*mQ*xn*
     &    apm3*fL6m3 + 2.D0*taucg**(-2)*taugs*qsq*mQ*cf*apm3*ffcg3 - 2.D
     &    0*taucg**(-2)*taugs**2*qsq*mQ*cf*apm4*ffcg3 - 1.D0/2.D0*
     &    taucg**(-1)*taugs*tcg*mQ*xn*apm1*fL6m3 + taucg**(-1)*taugs*
     &    qsq*mQ*xn**(-1)*apm4*ffcg3 - 2.D0*taucg**(-1)*taugs*mQ*cf*qcs
     &    *apm1 + 3.D0*taucg**(-1)*taugs*mQ*cf*qcs*apm4*ffcg3 +
     &    taucg**(-1)*taugs*mQ*cf*apm3*ffcg3 + taucg**(-1)*taugs*mQ**3*
     &    cf*apm5 + taucg**(-1)*taugs**2*mQ*cf*apm1 - taucg**(-1)*
     &    taugs**2*mQ*cf*apm4*ffcg1 - 3.D0/2.D0*taucg**(-1)*taugs**2*mQ
     &    *cf*apm5 - 1.D0/2.D0*taucg**(-1)*taugs**3*mQ**(-1)*cf*apm5*
     &    fq2m2
      spm = spm - 1.D0/2.D0*taucg**(-1)*qsq*mQ*xn**(-1)*apm3*fL6m2 -
     &    taucg**(-1)*xlog*mQ*b0*apm2 + 3.D0/2.D0*taucg**(-1)*mQ*epin*
     &    cf*apm2 + taucg**(-1)*mQ*epin*b0*apm2 - 1.D0/2.D0*taucg**(-1)
     &    *mQ*epin*xn**(-1)*apm2 + 1.D0/2.D0*taucg**(-1)*mQ*epin*xn*
     &    apm2 - 1.D0/2.D0*taucg**(-1)*mQ*epin2*xn**(-1)*apm2 + 3.D0/2.D
     &    0*taucg**(-1)*mQ*epin2*xn*apm2 - taucg**(-1)*mQ*cf*cgs1*apm3
     &     + 5.D0/2.D0*taucg**(-1)*mQ*cf*apm2 + taucg**(-1)*mQ*xn**(-1)
     &    *cgs1*apm3*fL6m2 - 1.D0/2.D0*taucg**(-1)*mQ*xn**(-1)*pisqo6*
     &    apm2 - 1.D0/2.D0*taucg**(-1)*mQ*xn*qcg*apm3*fL6m3 + 1.D0/2.D0
     &    *taucg**(-1)*mQ*xn*pisqo6*apm2 - 1.D0/2.D0*taucg*taugs*
     &    mQ**(-1)*cf*apm5*fq2m2 - taucg*mQ**(-1)*cf*cgs1*apm1 - 1.D0/2.
     &    D0*taugs**(-1)*mQ*xn**(-1)*qcs*cgs1*apm4*fL6m2 - 3.D0/2.D0*
     &    taugs*mQ*cf*apm5 + 1.D0/2.D0*taugs*mQ*xn*qcs*cgs1**(-1)*apm4*
     &    fL6m3 + 1.D0/2.D0*taugs*mQ*xn*qcg*cgs1**(-1)*apm1*fL6m3 + 1.D0
     &    /2.D0*taugs*mQ*xn*cgs1**(-1)*apm3*fL6m3 + 1.D0/2.D0*taugs*mQ*
     &    xn*apm4*fL6m3
      spm = spm - taugs**2*mQ**(-1)*cf*apm5*fq2m2 - ddilog(msq**(-1)*
     &    tcs)*taucg**(-1)*mQ*xn**(-1)*apm2 + ddilog(msq**(-1)*tcg)*
     &    taucg**(-1)*mQ*xn*apm2 + lnrat( - taucs,msq)*taucs*
     &    taucg**(-1)*taugs*tcs**(-1)*mQ**3*xn**(-1)*apm4 + lnrat( -
     &    taucs,msq)*taucs*taucg**(-1)*mQ*xn**(-1)*qcs*apm4 + lnrat( -
     &    taucs,msq)*taucs*taucg**(-1)*mQ*xn**(-1)*apm3 + lnrat( -
     &    taucs,msq)*taucs**2*taucg**(-1)*tcs**(-1)*qsq*mQ*xn**(-1)*
     &    apm1 - lnrat( - taucs,msq)*taucs**2*tcs**(-1)*mQ*xn**(-1)*
     &    apm4 + lnrat( - taucs,msq)*taucg**(-1)*mQ*epin*xn**(-1)*apm2
     &     + lnrat( - taucs,msq)*tcs**(-1)*mQ*xn**(-1)*cgs1*apm3 -
     &    lnrat( - taucs,msq)*tcs**(-1)*mQ**3*xn**(-1)*cgs1*apm4 -
     &    lnrat( - taucs,msq)**2*taucg**(-1)*mQ*xn**(-1)*apm2 - 2.D0*
     &    lnrat( - taucg,msq)*taucg**(-1)*taugs*tcg**(-1)*mQ*cf*qcg*
     &    apm3 - 4.D0*lnrat( - taucg,msq)*taucg**(-1)*taugs*mQ*cf*apm3
     &     + 2.D0*lnrat( - taucg,msq)*taucg**(-1)*taugs**2*tcg**(-1)*mQ
     &    *cf*qcg*apm4
      spm = spm + 4.D0*lnrat( - taucg,msq)*taucg**(-1)*taugs**2*
     & tcg**(-1)*mQ**3*cf*apm4 - lnrat( - taucg,msq)*taucg**(-1)*mQ*
     &    epin*xn*apm2 - lnrat( - taucg,msq)*taucg*tcg**(-1)*mQ*cf*qcg*
     &    apm1 - lnrat( - taucg,msq)*taugs*tcg**(-1)*qsq*mQ*xn**(-1)*
     &    apm4 + 2.D0*lnrat( - taucg,msq)*taugs*tcg**(-1)*mQ**3*cf*apm4
     &     + 3.D0*lnrat( - taucg,msq)*taugs*mQ*cf*apm1 + 2.D0*lnrat( -
     &    taucg,msq)*taugs**2*tcg**(-1)*mQ*cf*apm4 - lnrat( - taucg,msq
     &    )*tcg**(-1)*mQ*cf*qcs*apm3 - 2.D0*lnrat( - taucg,msq)*
     &    tcg**(-1)*mQ**3*cf*cgs1*apm1 - 2.D0*lnrat( - taucg,msq)*
     &    tcg**(-1)*mQ**3*cf*apm3 + lnrat( - taucg,msq)*mQ*xn**(-1)*qcg
     &    *apm1 + lnrat( - taucg,msq)*mQ*xn**(-1)*apm3 + lnrat( - taucg
     &    ,msq)**2*taucg**(-1)*mQ*xn*apm2 + lnrat( - taugs,msq)*taucs*
     &    taucg**(-1)*taugs*mQ*cf*apm1 + 1.D0/2.D0*lnrat( - taugs,msq)*
     &    taucs**2*taucg**(-1)*mQ**(-1)*cf*qcg*apm5 + lnrat( - taugs,
     &    msq)*taucs**2*taucg**(-1)*mQ**(-1)*cf*apm3 - lnrat( - taugs,
     &    msq)*taucs**2*taucg**(-1)*mQ*cf*apm1
      spm = spm + lnrat( - taugs,msq)*taucs**3*taucg**(-1)*mQ**(-1)*cf*
     & apm1 - 1.D0/2.D0*lnrat( - taugs,msq)*taucg**(-1)*taugs*mQ*cf*qcg
     &    *apm5 - lnrat( - taugs,msq)*taucg**(-1)*taugs*mQ*xn**(-1)*qcs
     &    *apm1 - lnrat( - taugs,msq)*taucg**(-1)*taugs*mQ*xn**(-1)*qcs
     &    *apm5 + lnrat( - taugs,msq)*taucg**(-1)*taugs*mQ*xn**(-1)*qcg
     &    *apm1 + lnrat( - taugs,msq)*taucg**(-1)*taugs*mQ*xn*apm3 -
     &    lnrat( - taugs,msq)*taucg**(-1)*taugs*mQ**3*cf*apm1 - lnrat(
     &     - taugs,msq)*taucg**(-1)*mQ*epin*xn*apm2 - 1.D0/2.D0*lnrat(
     &     - taugs,msq)*taugs*mQ**(-1)*xn**(-1)*cgs1*apm5 - 2.D0*lnrat(
     &     - taugs,msq)*taugs*mQ*cf*apm1 - 2.D0*lnrat( - taugs,msq)*
     &    taugs*mQ*cf*apm5 + 1.D0/2.D0*lnrat( - taugs,msq)**2*
     &    taucg**(-1)*mQ*xn*apm2 + 2.D0*lnrat( - qsqhat,msq)*taucs*
     &    taucg**(-2)*taugs*qsq**(-1)*mQ*qsqhat*cf*apm3 - lnrat( -
     &    qsqhat,msq)*taucs*taucg**(-2)*taugs**2*qsq**(-1)*mQ*qsqhat*cf
     &    *apm4 - 1.D0/2.D0*lnrat( - qsqhat,msq)*taucs*taucg**(-1)*
     &    taugs*qsq**(-1)*mQ*qsqhat*cf*apm1
      spm = spm - 5.D0/2.D0*lnrat( - qsqhat,msq)*taucs*taucg**(-1)*
     & taugs*qsq**(-1)*mQ*qsqhat*cf*apm4 - lnrat( - qsqhat,msq)*taucs*
     &    taucg**(-1)*tcs*qsq**(-1)*mQ*qsqhat*xn**(-1)*apm1 + lnrat( -
     &    qsqhat,msq)*taucs*taucg**(-1)*qsq**(-1)*mQ*qsqhat*cf*apm3 -
     &    lnrat( - qsqhat,msq)*taucs*taucg**(-1)*qsq**(-1)*mQ*qsqhat*xn
     &    *qcs*apm1 + 2.D0*lnrat( - qsqhat,msq)*taucs*taucg**(-1)*
     &    qsq**(-1)*mQ**3*qsqhat*cf*apm1 - lnrat( - qsqhat,msq)*taucs*
     &    taucg*qsq**(-1)*mQ**(-1)*qsqhat*cf*apm4 + 3.D0/2.D0*lnrat( -
     &    qsqhat,msq)*taucs*qsq**(-1)*mQ*qsqhat*cf*apm1 - 3.D0/2.D0*
     &    lnrat( - qsqhat,msq)*taucs*qsq**(-1)*mQ*qsqhat*cf*apm4 + 1.D0/
     &    2.D0*lnrat( - qsqhat,msq)*taucs**2*taucg**(-1)*taugs*
     &    qsq**(-1)*mQ**(-1)*qsqhat*cf*apm1 + 1.D0/2.D0*lnrat( - qsqhat
     &    ,msq)*taucs**2*taucg**(-1)*taugs*qsq**(-1)*mQ**(-1)*qsqhat*cf
     &    *apm4 - lnrat( - qsqhat,msq)*taucs**2*taucg**(-1)*qsq**(-1)*
     &    mQ**(-1)*qsqhat*cf*apm3 + 3.D0/2.D0*lnrat( - qsqhat,msq)*
     &    taucs**2*taucg**(-1)*qsq**(-1)*mQ*qsqhat*cf*apm1
      spm = spm + 5.D0/2.D0*lnrat( - qsqhat,msq)*taucs**2*taucg**(-1)*
     & qsq**(-1)*mQ*qsqhat*cf*apm4 - 1.D0/2.D0*lnrat( - qsqhat,msq)*
     &    taucs**3*taucg**(-1)*qsq**(-1)*mQ**(-1)*qsqhat*cf*apm1 + 1.D0/
     &    2.D0*lnrat( - qsqhat,msq)*taucs**3*taucg**(-1)*qsq**(-1)*
     &    mQ**(-1)*qsqhat*cf*apm4 + lnrat( - qsqhat,msq)*taucg**(-2)*
     &    taugs*qsq**(-1)*mQ**3*qsqhat*cf*apm3 + 2.D0*lnrat( - qsqhat,
     &    msq)*taucg**(-2)*taugs**2*qsq**(-1)*mQ*qsqhat*cf*apm3 - 2.D0*
     &    lnrat( - qsqhat,msq)*taucg**(-2)*taugs**2*qsq**(-1)*mQ**3*
     &    qsqhat*cf*apm4 - lnrat( - qsqhat,msq)*taucg**(-2)*taugs**3*
     &    qsq**(-1)*mQ*qsqhat*cf*apm4 - 1.D0/2.D0*lnrat( - qsqhat,msq)*
     &    taucg**(-1)*taugs*tcs*tcg*qsq**(-1)*mQ**(-1)*qsqhat*xn**(-1)*
     &    apm1 + 3.D0*lnrat( - qsqhat,msq)*taucg**(-1)*taugs*qsq**(-1)*
     &    mQ*qsqhat*cf*apm3 - 2.D0*lnrat( - qsqhat,msq)*taucg**(-1)*
     &    taugs*qsq**(-1)*mQ**3*qsqhat*cf*apm4 - 1.D0/2.D0*lnrat( -
     &    qsqhat,msq)*taucg**(-1)*taugs*mQ*qsqhat*xn*apm1 - 2.D0*lnrat(
     &     - qsqhat,msq)*taucg**(-1)*taugs**2*qsq**(-1)*mQ*qsqhat*cf*
     &    apm4
      spm = spm + 2.D0*lnrat( - qsqhat,msq)*taucg**(-1)*qsq**(-1)*mQ**3
     & *qsqhat*cf*apm3 - lnrat( - qsqhat,msq)*taucg**(-1)*mQ*qsqhat*
     &    xn**(-1)*apm3 - 1.D0/2.D0*lnrat( - qsqhat,msq)*taugs*
     &    qsq**(-1)*mQ**(-1)*qsqhat*xn**(-1)*cgs1*apm4 - 3.D0/2.D0*
     &    lnrat( - qsqhat,msq)*taugs*qsq**(-1)*mQ*qsqhat*cf*apm1 - 5.D0/
     &    2.D0*lnrat( - qsqhat,msq)*taugs*qsq**(-1)*mQ*qsqhat*cf*apm4
     &     + lnrat( - qsqhat,msq)*qsq**(-1)*mQ*qsqhat*cf*apm3 + 2.D0*
     &    C0fb2m(tcg,msq)*taucs**(-1)*taucg**(-2)*taugs**2*mQ**5*
     &    xn**(-1)*apm3 + C0fb2m(tcg,msq)*taucs**(-1)*taucg**(-1)*taugs
     &    *qsq*mQ**3*xn**(-1)*apm3 - C0fb2m(tcg,msq)*taucs**(-1)*
     &    taucg**(-1)*taugs*mQ**5*xn**(-1)*qcs*apm4 - C0fb2m(tcg,msq)*
     &    taucs*tcg*mQ*xn**(-1)*apm1 - 3.D0*C0fb2m(tcg,msq)*taucg**(-1)
     &    *taugs*mQ**3*xn**(-1)*apm3 - C0fb2m(tcg,msq)*mQ*xn**(-1)*qcs*
     &    apm3 - C0fb2m(tcg,msq)*mQ**3*xn**(-1)*cgs1*apm1 - 2.D0*
     &    C0fb2m(tcg,msq)*mQ**3*xn**(-1)*apm3 + 2.D0*C0fa2m(tcs,qsq,msq
     &    )*taucs**(-1)*taucg**(-2)*taugs*mQ**3*xn**(-1)*qcs*cgs1*apm3
      spm = spm + C0fa2m(tcs,qsq,msq)*taucs**(-1)*taucg**(-2)*taugs*
     & mQ**5*xn**(-1)*qcs**2*apm4 + C0fa2m(tcs,qsq,msq)*taucs**(-1)*
     &    taucg**(-1)*mQ*xn**(-1)*qcs**2*cgs1*apm3 + C0fa2m(tcs,qsq,msq
     &    )*taucs**(-1)*taucg**(-1)*mQ**3*xn**(-1)*qcs*cgs1*apm3 +
     &    C0fa2m(tcs,qsq,msq)*taucs*taucg**(-1)*tcg*mQ*xn**(-1)*qcs*
     &    apm1 + C0fa2m(tcs,qsq,msq)*taucg**(-1)*mQ**3*xn**(-1)*qcs*
     &    cgs1*apm1 + C0fa2m(tcs,qsq,msq)*taucg**(-1)*mQ**3*xn**(-1)*
     &    qcs*apm3 - 36.D0*ffc002*taucs*taucg**(-1)*taugs*mQ*cf*apm5 -
     &    18.D0*ffc002*taucs*taucg**(-1)*mQ*cf*qgs*apm5 + 18.D0*ffc002*
     &    taucs**2*taucg**(-1)*mQ**(-1)*cf*qgs*apm5 + 36.D0*ffc002*
     &    taucs**2*taucg**(-1)*mQ*cf*apm5 - 6.D0*ffc002*taucs**3*
     &    taucg**(-1)*taugs**(-1)*mQ**(-1)*cf*qgs*apm5 - 12.D0*ffc002*
     &    taucs**3*taucg**(-1)*mQ**(-1)*cf*apm5 + 6.D0*ffc002*
     &    taucg**(-1)*taugs*mQ*cf*qgs*apm5 + 12.D0*ffc002*taucg**(-1)*
     &    taugs*mQ**3*cf*apm5 - 12.D0*ffc00*taucs*taucg**(-1)*taugs*mQ*
     &    cf*apm5
      spm = spm + 8.D0*ffc00*taucs*taucg**(-1)*mQ*cf*apm3 - 20.D0*ffc00
     &    *taucs*taucg**(-1)*mQ**3*cf*apm1 + 8.D0*ffc00*taucs*
     &    taugs**(-1)*mQ*cf*apm3 - 8.D0*ffc00*taucs*mQ*cf*apm1 - 8.D0*
     &    ffc00*taucs*mQ*cf*apm5 + 4.D0*ffc00*taucs**2*taucg**(-1)*
     &    taugs**(-1)*mQ*cf*qcg*apm1 + 8.D0*ffc00*taucs**2*taucg**(-1)*
     &    mQ*cf*apm1 + 28.D0*ffc00*taucs**2*taucg**(-1)*mQ*cf*apm5 - 4.D
     &    0*ffc00*taucs**2*taugs**(-1)*mQ**(-1)*cf*apm3 + 12.D0*ffc00*
     &    taucs**2*taugs**(-1)*mQ*cf*apm1 + 4.D0*ffc00*taucs**2*
     &    mQ**(-1)*cf*apm5 - 4.D0*ffc00*taucs**3*taucg**(-1)*
     &    taugs**(-1)*mQ**(-1)*cf*qgs*apm1 - 4.D0*ffc00*taucs**3*
     &    taucg**(-1)*taugs**(-1)*mQ**(-1)*cf*qgs*apm5 - 4.D0*ffc00*
     &    taucs**3*taucg**(-1)*taugs**(-1)*mQ**(-1)*cf*apm3 + 4.D0*
     &    ffc00*taucg**(-1)*taugs*mQ**3*cf*apm1 - 8.D0*ffc00*
     &    taucg**(-1)*mQ**3*cf*apm3 - 4.D0*ffc00*taugs*mQ*cf*apm5 + 2.D0
     &    *ffc00*mQ**(-1)*xn**(-1)*qgs*cgs1*apm5 + 4.D0*ffc00*mQ*
     &    xn**(-1)*cgs1*apm5
      spm = spm + 4.D0*ffc00*mQ**3*cf*apm1 - 2.D0*fI3me*taucs*
     &    taucg**(-1)*taugs*mQ*cf*apm3 - fI3me*taucs*taucg**(-1)*taugs*
     &    mQ**3*xn**(-1)*apm1 - 2.D0*fI3me*taucs*taugs*mQ*cf*apm4 + 2.D0
     &    *fI3me*taucs**2*taucg**(-1)*taugs*mQ*cf*apm4 + 2.D0*fI3me*
     &    taucg**(-2)*taugs*mQ**3*xn*cgs1*apm3 + fI3me*taucg**(-2)*
     &    taugs**2*mQ**5*xn*apm1 + fI3me*taucg**(-1)*taugs*mQ*xn*cgs1*
     &    apm3 + fI3me*taugs*mQ*xn*cgs1*apm4


      lopm =-mQ/taucg*za(ig,ie)*zb(is,in)*zb(is,ic)/zb(ig,is)

c--- include finite counterterm to go from DR to MSbar scheme
c--- alphas(DR) = alphas(MSbar) * (1+ (Nc / 6) * alphas(MSbar) / (2*pi))
      spm=spm + lopm * xn/six

c--- include finite counterterm due to FDH scheme
c--- gw = gw * ( 1 - 2*cf * alphas(MSbar) / (2*pi))
      spm=spm - lopm * cf * two

      virt_pm=spm

      return
      end
