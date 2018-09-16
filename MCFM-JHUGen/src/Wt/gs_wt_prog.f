      subroutine gs_wt_prog(mq,qwidth,p,ig,is,ie,in,jn,je,jb,ia,gs)
      implicit none
      include 'types.f'
c---- Amplitudes for the proces
c     g(ig)+s(is)-->W^-{e^-(ie)+nbar(in)}+t{W^+[n(jn)+e^+(je)]+b(jb)}+a(ia)
c---- helicities: gs(ha,hg)
c--- label on amplitudes represent gluon helicities for ia,ig
c---   1 = negative helicity, 2 = positive helicity    
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: is,ig,ie,in,ia,je,jn,jb
      real(dp):: p(mxpart,4),tsq,mq,qwidth
      real(dp):: dot,msq,aDt,aDg,aDs,tDg,gDs,sga2,tga2
      complex(dp):: gs(2,2)
      msq=mq**2
      aDt=dot(p,ia,jn)+dot(p,ia,je)+dot(p,ia,jb)
      tDg=dot(p,ig,jn)+dot(p,ig,je)+dot(p,ig,jb)
      aDg=dot(p,ia,ig)
      aDs=dot(p,ia,is)
      gDs=dot(p,ig,is)
      sga2=two*(gDs+aDs+aDg)
      tga2=two*(tDg+aDt+aDg)

c--- additions to cope with off-shell top quark
      tsq=two*(dot(p,jn,je)+dot(p,jn,jb)+dot(p,je,jb))
      tDg=tDg+(tsq-msq)/two
      tga2=tga2+(tsq-msq)

      call spinoru(8,p,za,zb)
      
c--- Do not use the "overall scheme" to include the top width when necessary
c      if (real(tga2+mt**2) > zero) then
c        tga2=tga2+cplx2(zero,mt*twidth)
c        tga2=sqrt(real(tga2**2+(mt*twidth)**2))
c      endif
      
      gs(1,1)=  + tDg**(-1) * ( za(jb,jn)*za(ig,ia)*za(ig,jn)*za(ie,ia)
     &    *zb(is,in)*zb(jn,je) + za(jb,jn)*za(ig,jb)*za(ig,ia)*za(ie,ia
     &    )*zb(jb,je)*zb(is,in) )
      gs(1,1) = gs(1,1) + 1/(zb(ig,ia))*msq*tDg**(-1) * ( za(jb,jn)*za(
     &    ig,ia)*za(ie,ia)*zb(is,in)*zb(ia,je) + za(jb,jn)*za(ig,jn)*
     &    za(ie,ia)*zb(is,in)*zb(jn,je) + za(jb,jn)*za(ig,jb)*za(ie,ia)
     &    *zb(jb,je)*zb(is,in) )
      gs(1,1) = gs(1,1) + 1/(zb(ig,ia))*tga2*aDs**(-1)*tDg**(-1) * (
     &    half*za(jb,jn)*za(ig,ie)*za(ig,jn)*za(is,ia)*zb(ig,is)*zb(
     &    is,in)*zb(jn,je) + half*za(jb,jn)*za(ig,jb)*za(ig,ie)*
     &    za(is,ia)*zb(jb,je)*zb(ig,is)*zb(is,in) )
      gs(1,1) = gs(1,1) + 1/(zb(ig,ia))*tDg**(-1) * ( za(jb,ia)*za(jb,
     &    jn)*za(ig,jn)*za(ie,ia)*zb(jb,ia)*zb(is,in)*zb(jn,je) + za(jb
     &    ,ia)*za(jb,jn)*za(ig,jb)*za(ie,ia)*zb(jb,ia)*zb(jb,je)*zb(is,
     &    in) + za(jb,jn)*za(ig,ia)*za(ig,jn)*za(ie,jn)*zb(ig,jn)*zb(is
     &    ,in)*zb(jn,je) + za(jb,jn)*za(ig,ia)*za(ig,jn)*za(ie,je)*zb(
     &    ig,je)*zb(is,in)*zb(jn,je) + za(jb,jn)*za(ig,ia)*za(ig,jn)*
     &    za(ie,jb)*zb(ig,jb)*zb(is,in)*zb(jn,je) + za(jb,jn)*za(ig,jn)
     &    *za(ie,ia)*za(ia,jn)*zb(is,in)*zb(ia,jn)*zb(jn,je) + za(jb,jn
     &    )*za(ig,jn)*za(ie,ia)*za(ia,je)*zb(is,in)*zb(ia,je)*zb(jn,je)
     &     + za(jb,jn)*za(ig,jb)*za(ig,ia)*za(ie,jn)*zb(jb,je)*zb(ig,jn
     &    )*zb(is,in) + za(jb,jn)*za(ig,jb)*za(ig,ia)*za(ie,je)*zb(jb,
     &    je)*zb(ig,je)*zb(is,in) + za(jb,jn)*za(ig,jb)*za(ig,ia)*za(ie
     &    ,jb)*zb(jb,je)*zb(ig,jb)*zb(is,in) + za(jb,jn)*za(ig,jb)*za(
     &    ie,ia)*za(ia,jn)*zb(jb,je)*zb(is,in)*zb(ia,jn) + za(jb,jn)*
     &    za(ig,jb)*za(ie,ia)*za(ia,je)*zb(jb,je)*zb(is,in)*zb(ia,je) )
      gs(1,1) = gs(1,1) + 1/(zb(ig,ia))/(zb(ig,ia))*msq*tga2*aDs**(-1)*
     & tDg**(-1) * ( half*za(jb,jn)*za(ig,ie)*za(is,ia)*zb(ig,is)*
     &    zb(is,in)*zb(ia,je) )
      gs(1,1) = gs(1,1) + 1/(zb(ig,ia))/(zb(ig,ia))*msq*tDg**(-1) * ( 
     &    za(jb,jn)*za(ig,ia)*za(ie,jn)*zb(ig,jn)*zb(is,in)*zb(ia,je)
     &     + za(jb,jn)*za(ig,ia)*za(ie,je)*zb(ig,je)*zb(is,in)*zb(ia,je
     &    ) + za(jb,jn)*za(ig,ia)*za(ie,jb)*zb(ig,jb)*zb(is,in)*zb(ia,
     &    je) - za(jb,jn)*za(ig,jn)*za(ie,ia)*zb(ig,jn)*zb(is,in)*zb(ia
     &    ,je) - za(jb,jn)*za(ig,je)*za(ie,ia)*zb(ig,je)*zb(is,in)*zb(
     &    ia,je) - za(jb,jn)*za(ig,jb)*za(ie,ia)*zb(ig,jb)*zb(is,in)*
     &    zb(ia,je) )
      gs(1,1) = gs(1,1) + 1/(zb(ig,ia))/(zb(ig,ia))*msq * (  - za(jb,jn
     &    )*za(ig,ie)*zb(ig,je)*zb(is,in) - za(jb,jn)*za(ie,ia)*zb(is,
     &    in)*zb(ia,je) )
      gs(1,1) = gs(1,1) + 1/(zb(ig,ia))/(zb(ig,ia))*sga2**(-1)*tga2*
     & aDs**(-1) * ( za(jb,jn)*za(ig,is)*za(is,ia)*za(ie,jn)*zb(ig,is)*
     &    zb(is,in)*zb(is,ia)*zb(jn,je) + za(jb,jn)*za(ig,is)*za(is,ia)
     &    *za(ie,jb)*zb(jb,je)*zb(ig,is)*zb(is,in)*zb(is,ia) - za(jb,jn
     &    )*za(ig,ia)*za(is,ia)*za(ie,jn)*zb(ig,is)*zb(is,ia)*zb(in,ia)
     &    *zb(jn,je) - za(jb,jn)*za(ig,ia)*za(is,ia)*za(ie,jb)*zb(jb,je
     &    )*zb(ig,is)*zb(is,ia)*zb(in,ia) )
      gs(1,1) = gs(1,1) + 1/(zb(ig,ia))/(zb(ig,ia))*sga2**(-1)*tga2
     &  * ( za(jb,jn)*za(ig,is)*za(ie,jn)*zb(ig,is)*zb(is,in)*zb(jn,je)
     &     + za(jb,jn)*za(ig,is)*za(ie,jb)*zb(jb,je)*zb(ig,is)*zb(is,in
     &    ) - za(jb,jn)*za(ig,ia)*za(ie,jn)*zb(ig,is)*zb(in,ia)*zb(jn,
     &    je) - za(jb,jn)*za(ig,ia)*za(ie,jn)*zb(ig,in)*zb(is,ia)*zb(jn
     &    ,je) - za(jb,jn)*za(ig,ia)*za(ie,jb)*zb(jb,je)*zb(ig,is)*zb(
     &    in,ia) - za(jb,jn)*za(ig,ia)*za(ie,jb)*zb(jb,je)*zb(ig,in)*
     &    zb(is,ia) - za(jb,jn)*za(is,ia)*za(ie,jn)*zb(is,in)*zb(is,ia)
     &    *zb(jn,je) - za(jb,jn)*za(is,ia)*za(ie,jb)*zb(jb,je)*zb(is,in
     &    )*zb(is,ia) )
      gs(1,1) = gs(1,1) + 1/(zb(ig,ia))/(zb(ig,ia))*tga2*aDs**(-1)*
     & tDg**(-1) * ( half*za(jb,jn)*za(ig,jn)*za(is,ia)*za(ie,jn)*
     &    zb(ig,is)*zb(is,in)*zb(ia,jn)*zb(jn,je) + half*za(jb,jn)
     &    *za(ig,jn)*za(is,ia)*za(ie,je)*zb(ig,is)*zb(is,in)*zb(ia,je)*
     &    zb(jn,je) - half*za(jb,jn)*za(ig,jn)*za(is,ia)*za(ie,jb)
     &    *zb(jb,ia)*zb(ig,is)*zb(is,in)*zb(jn,je) + half*za(jb,jn
     &    )*za(ig,jb)*za(is,ia)*za(ie,jn)*zb(jb,je)*zb(ig,is)*zb(is,in)
     &    *zb(ia,jn) + half*za(jb,jn)*za(ig,jb)*za(is,ia)*za(ie,je
     &    )*zb(jb,je)*zb(ig,is)*zb(is,in)*zb(ia,je) - half*za(jb,
     &    jn)*za(ig,jb)*za(is,ia)*za(ie,jb)*zb(jb,ia)*zb(jb,je)*zb(ig,
     &    is)*zb(is,in) )
      gs(1,1) = gs(1,1) + 1/(zb(ig,ia))/(zb(ig,ia))*tDg**(-1) * ( za(jb
     &    ,ia)*za(jb,jn)*za(ig,jn)*za(ie,jn)*zb(jb,ia)*zb(ig,jn)*zb(is,
     &    in)*zb(jn,je) + za(jb,ia)*za(jb,jn)*za(ig,jn)*za(ie,je)*zb(jb
     &    ,ia)*zb(ig,je)*zb(is,in)*zb(jn,je) + za(jb,ia)*za(jb,jn)*za(
     &    ig,jn)*za(ie,jb)*zb(jb,ia)*zb(ig,jb)*zb(is,in)*zb(jn,je) + 
     &    za(jb,ia)*za(jb,jn)*za(ig,jb)*za(ie,jn)*zb(jb,ia)*zb(jb,je)*
     &    zb(ig,jn)*zb(is,in) + za(jb,ia)*za(jb,jn)*za(ig,jb)*za(ie,je)
     &    *zb(jb,ia)*zb(jb,je)*zb(ig,je)*zb(is,in) + za(jb,ia)*za(jb,jn
     &    )*za(ig,jb)*za(ie,jb)*zb(jb,ia)*zb(jb,je)*zb(ig,jb)*zb(is,in)
     &     + za(jb,jn)*za(ig,jn)*za(ie,jn)*za(ia,jn)*zb(ig,jn)*zb(is,in
     &    )*zb(ia,jn)*zb(jn,je) + za(jb,jn)*za(ig,jn)*za(ie,jn)*za(ia,
     &    je)*zb(ig,jn)*zb(is,in)*zb(ia,je)*zb(jn,je) + za(jb,jn)*za(ig
     &    ,jn)*za(ie,je)*za(ia,jn)*zb(ig,je)*zb(is,in)*zb(ia,jn)*zb(jn,
     &    je) + za(jb,jn)*za(ig,jn)*za(ie,je)*za(ia,je)*zb(ig,je)*zb(is
     &    ,in)*zb(ia,je)*zb(jn,je) + za(jb,jn)*za(ig,jn)*za(ie,jb)*za(
     &    ia,jn)*zb(ig,jb)*zb(is,in)*zb(ia,jn)*zb(jn,je) )
      gs(1,1) = gs(1,1) + 1/(zb(ig,ia))/(zb(ig,ia))*tDg**(-1) * ( za(jb
     &    ,jn)*za(ig,jn)*za(ie,jb)*za(ia,je)*zb(ig,jb)*zb(is,in)*zb(ia,
     &    je)*zb(jn,je) + za(jb,jn)*za(ig,jb)*za(ie,jn)*za(ia,jn)*zb(jb
     &    ,je)*zb(ig,jn)*zb(is,in)*zb(ia,jn) + za(jb,jn)*za(ig,jb)*za(
     &    ie,jn)*za(ia,je)*zb(jb,je)*zb(ig,jn)*zb(is,in)*zb(ia,je) + 
     &    za(jb,jn)*za(ig,jb)*za(ie,je)*za(ia,jn)*zb(jb,je)*zb(ig,je)*
     &    zb(is,in)*zb(ia,jn) + za(jb,jn)*za(ig,jb)*za(ie,je)*za(ia,je)
     &    *zb(jb,je)*zb(ig,je)*zb(is,in)*zb(ia,je) + za(jb,jn)*za(ig,jb
     &    )*za(ie,jb)*za(ia,jn)*zb(jb,je)*zb(ig,jb)*zb(is,in)*zb(ia,jn)
     &     + za(jb,jn)*za(ig,jb)*za(ie,jb)*za(ia,je)*zb(jb,je)*zb(ig,jb
     &    )*zb(is,in)*zb(ia,je) )
      gs(1,1) = gs(1,1) + 1/(zb(ig,ia))/(zb(ig,ia)) * (  - za(jb,ia)*
     &    za(jb,jn)*za(ie,jn)*zb(jb,je)*zb(is,in)*zb(ia,jn) - za(jb,ia)
     &    *za(jb,jn)*za(ie,je)*zb(jb,je)*zb(is,in)*zb(ia,je) + za(jb,ia
     &    )*za(jb,jn)*za(ie,jb)*zb(jb,ia)*zb(jb,je)*zb(is,in) - za(jb,
     &    jn)*za(ig,jn)*za(ie,jn)*zb(ig,jn)*zb(is,in)*zb(jn,je) - za(jb
     &    ,jn)*za(ig,jn)*za(ie,je)*zb(ig,je)*zb(is,in)*zb(jn,je) - za(
     &    jb,jn)*za(ig,jn)*za(ie,jb)*zb(ig,jb)*zb(is,in)*zb(jn,je) - 
     &    za(jb,jn)*za(ig,jb)*za(ie,jn)*zb(jb,je)*zb(ig,jn)*zb(is,in)
     &     - za(jb,jn)*za(ig,jb)*za(ie,je)*zb(jb,je)*zb(ig,je)*zb(is,in
     &    ) - za(jb,jn)*za(ig,jb)*za(ie,jb)*zb(jb,je)*zb(ig,jb)*zb(is,
     &    in) + za(jb,jn)*za(ie,jn)*za(ia,jn)*zb(is,in)*zb(ia,jn)*zb(jn
     &    ,je) + za(jb,jn)*za(ie,je)*za(ia,jn)*zb(is,in)*zb(ia,je)*zb(
     &    jn,je) - za(jb,jn)*za(ie,jb)*za(ia,jn)*zb(jb,ia)*zb(is,in)*
     &    zb(jn,je) )
      gs(1,1) = gs(1,1) + 1/(zb(ig,ia)) * (  - za(jb,ia)*za(jb,jn)*za(
     &    ig,ie)*zb(jb,je)*zb(is,in) + za(jb,jn)*za(ig,ie)*za(ia,jn)*
     &    zb(is,in)*zb(jn,je) - za(jb,jn)*za(ig,jn)*za(ie,ia)*zb(is,in)
     &    *zb(jn,je) - za(jb,jn)*za(ig,jb)*za(ie,ia)*zb(jb,je)*zb(is,in
     &    ) )

      gs(1,2)=  + 1/(za(ig,ia))*tDg**(-1) * (  - za(jb,ia)**2*za(jb,jn)
     &    *za(ie,ia)*zb(jb,je)*zb(ig,jb)*zb(is,in) + za(jb,ia)*za(jb,jn
     &    )*za(ie,ia)*za(ia,jn)*zb(ig,jb)*zb(is,in)*zb(jn,je) + za(jb,
     &    ia)*za(jb,jn)*za(ie,ia)*za(ia,jn)*zb(jb,je)*zb(ig,jn)*zb(is,
     &    in) + za(jb,ia)*za(jb,jn)*za(ie,ia)*za(ia,je)*zb(jb,je)*zb(ig
     &    ,je)*zb(is,in) - za(jb,jn)*za(ie,ia)*za(ia,jn)**2*zb(ig,jn)*
     &    zb(is,in)*zb(jn,je) - za(jb,jn)*za(ie,ia)*za(ia,jn)*za(ia,je)
     &    *zb(ig,je)*zb(is,in)*zb(jn,je) )
      gs(1,2) = gs(1,2) + 1/(za(ig,ia))/(zb(ig,ia))*msq*tga2*aDs**(-1)*
     & tDg**(-1) * ( half*za(jb,jn)*za(is,ia)*za(ie,ia)*zb(ig,is)*
     &    zb(ig,je)*zb(is,in) )
      gs(1,2) = gs(1,2) + 1/(za(ig,ia))/(zb(ig,ia))*msq*tDg**(-1) * ( 
     &     - za(jb,ia)*za(jb,jn)*za(ie,ia)*zb(ig,jb)*zb(ig,je)*zb(is,in
     &    ) + za(jb,jn)*za(ie,ia)*za(ia,jn)*zb(ig,jn)*zb(ig,je)*zb(is,
     &    in) + za(jb,jn)*za(ie,ia)*za(ia,je)*zb(ig,je)**2*zb(is,in) )
      gs(1,2) = gs(1,2) + 1/(za(ig,ia))/(zb(ig,ia))*sga2**(-1)*tga2*
     & aDs**(-1) * (  - za(jb,jn)*za(is,ia)**2*za(ie,jn)*zb(ig,is)**2*
     &    zb(is,in)*zb(jn,je) - za(jb,jn)*za(is,ia)**2*za(ie,jb)*zb(jb,
     &    je)*zb(ig,is)**2*zb(is,in) )
      gs(1,2) = gs(1,2) + 1/(za(ig,ia))/(zb(ig,ia))*tga2*aDs**(-1)*
     & tDg**(-1) * ( half*za(jb,ia)*za(jb,jn)*za(is,ia)*za(ie,jn)*
     &    zb(jb,je)*zb(ig,is)*zb(ig,jn)*zb(is,in) + half*za(jb,ia)
     &    *za(jb,jn)*za(is,ia)*za(ie,je)*zb(jb,je)*zb(ig,is)*zb(ig,je)*
     &    zb(is,in) + half*za(jb,ia)*za(jb,jn)*za(is,ia)*za(ie,jb)
     &    *zb(jb,je)*zb(ig,jb)*zb(ig,is)*zb(is,in) - half*za(jb,jn
     &    )*za(is,ia)*za(ie,jn)*za(ia,jn)*zb(ig,is)*zb(ig,jn)*zb(is,in)
     &    *zb(jn,je) - half*za(jb,jn)*za(is,ia)*za(ie,je)*za(ia,jn
     &    )*zb(ig,is)*zb(ig,je)*zb(is,in)*zb(jn,je) - half*za(jb,
     &    jn)*za(is,ia)*za(ie,jb)*za(ia,jn)*zb(ig,jb)*zb(ig,is)*zb(is,
     &    in)*zb(jn,je) )
      gs(1,2) = gs(1,2) + 1/(za(ig,ia))/(zb(ig,ia))*tDg**(-1) * (  - 
     &    za(jb,ia)**2*za(jb,jn)*za(ie,jn)*zb(jb,je)*zb(ig,jb)*zb(ig,jn
     &    )*zb(is,in) - za(jb,ia)**2*za(jb,jn)*za(ie,je)*zb(jb,je)*zb(
     &    ig,jb)*zb(ig,je)*zb(is,in) - za(jb,ia)**2*za(jb,jn)*za(ie,jb)
     &    *zb(jb,je)*zb(ig,jb)**2*zb(is,in) + za(jb,ia)*za(jb,jn)*za(ie
     &    ,jn)*za(ia,jn)*zb(ig,jb)*zb(ig,jn)*zb(is,in)*zb(jn,je) + za(
     &    jb,ia)*za(jb,jn)*za(ie,jn)*za(ia,jn)*zb(jb,je)*zb(ig,jn)**2*
     &    zb(is,in) + za(jb,ia)*za(jb,jn)*za(ie,jn)*za(ia,je)*zb(jb,je)
     &    *zb(ig,jn)*zb(ig,je)*zb(is,in) + za(jb,ia)*za(jb,jn)*za(ie,je
     &    )*za(ia,jn)*zb(ig,jb)*zb(ig,je)*zb(is,in)*zb(jn,je) + za(jb,
     &    ia)*za(jb,jn)*za(ie,je)*za(ia,jn)*zb(jb,je)*zb(ig,jn)*zb(ig,
     &    je)*zb(is,in) + za(jb,ia)*za(jb,jn)*za(ie,je)*za(ia,je)*zb(jb
     &    ,je)*zb(ig,je)**2*zb(is,in) + za(jb,ia)*za(jb,jn)*za(ie,jb)*
     &    za(ia,jn)*zb(ig,jb)**2*zb(is,in)*zb(jn,je) + za(jb,ia)*za(jb,
     &    jn)*za(ie,jb)*za(ia,jn)*zb(jb,je)*zb(ig,jb)*zb(ig,jn)*zb(is,
     &    in) )
      gs(1,2) = gs(1,2) + 1/(za(ig,ia))/(zb(ig,ia))*tDg**(-1) * ( za(jb
     &    ,ia)*za(jb,jn)*za(ie,jb)*za(ia,je)*zb(jb,je)*zb(ig,jb)*zb(ig,
     &    je)*zb(is,in) - za(jb,jn)*za(ie,jn)*za(ia,jn)**2*zb(ig,jn)**2
     &    *zb(is,in)*zb(jn,je) - za(jb,jn)*za(ie,jn)*za(ia,jn)*za(ia,je
     &    )*zb(ig,jn)*zb(ig,je)*zb(is,in)*zb(jn,je) - za(jb,jn)*za(ie,
     &    je)*za(ia,jn)**2*zb(ig,jn)*zb(ig,je)*zb(is,in)*zb(jn,je) - 
     &    za(jb,jn)*za(ie,je)*za(ia,jn)*za(ia,je)*zb(ig,je)**2*zb(is,in
     &    )*zb(jn,je) - za(jb,jn)*za(ie,jb)*za(ia,jn)**2*zb(ig,jb)*zb(
     &    ig,jn)*zb(is,in)*zb(jn,je) - za(jb,jn)*za(ie,jb)*za(ia,jn)*
     &    za(ia,je)*zb(ig,jb)*zb(ig,je)*zb(is,in)*zb(jn,je) )
      gs(1,2) = gs(1,2) + 1/(zb(ig,ia))*sga2**(-1)*tga2*aDs**(-1) * ( 
     &     - za(jb,jn)*za(is,ia)*za(ie,jn)*zb(ig,is)**2*zb(ig,in)*zb(jn
     &    ,je) - za(jb,jn)*za(is,ia)*za(ie,jb)*zb(jb,je)*zb(ig,is)**2*
     &    zb(ig,in) )

      gs(2,1)=  + tga2*aDs**(-1)*tDg**(-1) * ( half*za(jb,jn)*za(
     &    ig,ie)*za(ig,jn)*zb(is,ia)*zb(in,ia)*zb(jn,je) + half*
     &    za(jb,jn)*za(ig,jb)*za(ig,ie)*zb(jb,je)*zb(is,ia)*zb(in,ia) )
      gs(2,1) = gs(2,1) + 1/(za(ig,ia))*tga2*aDs**(-1)*tDg**(-1) * ( 
     &     - half*za(jb,jn)*za(ig,is)*za(ig,ie)*za(ig,jn)*zb(is,in
     &    )*zb(is,ia)*zb(jn,je) - half*za(jb,jn)*za(ig,jb)*za(ig,
     &    is)*za(ig,ie)*zb(jb,je)*zb(is,in)*zb(is,ia) )
      gs(2,1) = gs(2,1) + 1/(za(ig,ia))*tDg**(-1) * (  - za(jb,jn)*za(
     &    ig,ie)*za(ig,jn)**2*zb(is,in)*zb(ia,jn)*zb(jn,je) - za(jb,jn)
     &    *za(ig,ie)*za(ig,jn)*za(ig,je)*zb(is,in)*zb(ia,je)*zb(jn,je)
     &     + za(jb,jn)*za(ig,jb)*za(ig,ie)*za(ig,jn)*zb(jb,ia)*zb(is,in
     &    )*zb(jn,je) - za(jb,jn)*za(ig,jb)*za(ig,ie)*za(ig,jn)*zb(jb,
     &    je)*zb(is,in)*zb(ia,jn) - za(jb,jn)*za(ig,jb)*za(ig,ie)*za(ig
     &    ,je)*zb(jb,je)*zb(is,in)*zb(ia,je) + za(jb,jn)*za(ig,jb)**2*
     &    za(ig,ie)*zb(jb,ia)*zb(jb,je)*zb(is,in) )
      gs(2,1) = gs(2,1) + 1/(za(ig,ia))/(zb(ig,ia))*msq*tga2*aDs**(-1)*
     & tDg**(-1) * (  - half*za(jb,jn)*za(ig,is)*za(ig,ie)*zb(is,
     &    in)*zb(is,ia)*zb(ia,je) )
      gs(2,1) = gs(2,1) + 1/(za(ig,ia))/(zb(ig,ia))*msq*tDg**(-1) * ( 
     &     - za(jb,jn)*za(ig,ie)*za(ig,jn)*zb(is,in)*zb(ia,jn)*zb(ia,je
     &    ) - za(jb,jn)*za(ig,ie)*za(ig,je)*zb(is,in)*zb(ia,je)**2 + 
     &    za(jb,jn)*za(ig,jb)*za(ig,ie)*zb(jb,ia)*zb(is,in)*zb(ia,je) )
      gs(2,1) = gs(2,1) + 1/(za(ig,ia))/(zb(ig,ia))*sga2**(-1)*tga2*
     & aDs**(-1) * (  - za(jb,jn)*za(ig,is)**2*za(ie,jn)*zb(is,in)*zb(
     &    is,ia)**2*zb(jn,je) - za(jb,jn)*za(ig,is)**2*za(ie,jb)*zb(jb,
     &    je)*zb(is,in)*zb(is,ia)**2 )
      gs(2,1) = gs(2,1) + 1/(za(ig,ia))/(zb(ig,ia))*tga2*aDs**(-1)*
     & tDg**(-1) * (  - half*za(jb,jn)*za(ig,is)*za(ig,jn)*za(ie,
     &    jn)*zb(is,in)*zb(is,ia)*zb(ia,jn)*zb(jn,je) - half*za(jb
     &    ,jn)*za(ig,is)*za(ig,jn)*za(ie,je)*zb(is,in)*zb(is,ia)*zb(ia,
     &    je)*zb(jn,je) + half*za(jb,jn)*za(ig,is)*za(ig,jn)*za(ie
     &    ,jb)*zb(jb,ia)*zb(is,in)*zb(is,ia)*zb(jn,je) - half*za(
     &    jb,jn)*za(ig,jb)*za(ig,is)*za(ie,jn)*zb(jb,je)*zb(is,in)*zb(
     &    is,ia)*zb(ia,jn) - half*za(jb,jn)*za(ig,jb)*za(ig,is)*
     &    za(ie,je)*zb(jb,je)*zb(is,in)*zb(is,ia)*zb(ia,je) + half
     &    *za(jb,jn)*za(ig,jb)*za(ig,is)*za(ie,jb)*zb(jb,ia)*zb(jb,je)*
     &    zb(is,in)*zb(is,ia) )
      gs(2,1) = gs(2,1) + 1/(za(ig,ia))/(zb(ig,ia))*tDg**(-1) * (  - 
     &    za(jb,jn)*za(ig,jn)**2*za(ie,jn)*zb(is,in)*zb(ia,jn)**2*zb(jn
     &    ,je) - za(jb,jn)*za(ig,jn)**2*za(ie,je)*zb(is,in)*zb(ia,jn)*
     &    zb(ia,je)*zb(jn,je) + za(jb,jn)*za(ig,jn)**2*za(ie,jb)*zb(jb,
     &    ia)*zb(is,in)*zb(ia,jn)*zb(jn,je) - za(jb,jn)*za(ig,jn)*za(ig
     &    ,je)*za(ie,jn)*zb(is,in)*zb(ia,jn)*zb(ia,je)*zb(jn,je) - za(
     &    jb,jn)*za(ig,jn)*za(ig,je)*za(ie,je)*zb(is,in)*zb(ia,je)**2*
     &    zb(jn,je) + za(jb,jn)*za(ig,jn)*za(ig,je)*za(ie,jb)*zb(jb,ia)
     &    *zb(is,in)*zb(ia,je)*zb(jn,je) + za(jb,jn)*za(ig,jb)*za(ig,jn
     &    )*za(ie,jn)*zb(jb,ia)*zb(is,in)*zb(ia,jn)*zb(jn,je) - za(jb,
     &    jn)*za(ig,jb)*za(ig,jn)*za(ie,jn)*zb(jb,je)*zb(is,in)*zb(ia,
     &    jn)**2 + za(jb,jn)*za(ig,jb)*za(ig,jn)*za(ie,je)*zb(jb,ia)*
     &    zb(is,in)*zb(ia,je)*zb(jn,je) - za(jb,jn)*za(ig,jb)*za(ig,jn)
     &    *za(ie,je)*zb(jb,je)*zb(is,in)*zb(ia,jn)*zb(ia,je) - za(jb,jn
     &    )*za(ig,jb)*za(ig,jn)*za(ie,jb)*zb(jb,ia)**2*zb(is,in)*zb(jn,
     &    je) )
      gs(2,1) = gs(2,1) + 1/(za(ig,ia))/(zb(ig,ia))*tDg**(-1) * ( za(jb
     &    ,jn)*za(ig,jb)*za(ig,jn)*za(ie,jb)*zb(jb,ia)*zb(jb,je)*zb(is,
     &    in)*zb(ia,jn) - za(jb,jn)*za(ig,jb)*za(ig,je)*za(ie,jn)*zb(jb
     &    ,je)*zb(is,in)*zb(ia,jn)*zb(ia,je) - za(jb,jn)*za(ig,jb)*za(
     &    ig,je)*za(ie,je)*zb(jb,je)*zb(is,in)*zb(ia,je)**2 + za(jb,jn)
     &    *za(ig,jb)*za(ig,je)*za(ie,jb)*zb(jb,ia)*zb(jb,je)*zb(is,in)*
     &    zb(ia,je) + za(jb,jn)*za(ig,jb)**2*za(ie,jn)*zb(jb,ia)*zb(jb,
     &    je)*zb(is,in)*zb(ia,jn) + za(jb,jn)*za(ig,jb)**2*za(ie,je)*
     &    zb(jb,ia)*zb(jb,je)*zb(is,in)*zb(ia,je) - za(jb,jn)*za(ig,jb)
     &    **2*za(ie,jb)*zb(jb,ia)**2*zb(jb,je)*zb(is,in) )
      gs(2,1) = gs(2,1) + 1/(zb(ig,ia))*msq*tga2*aDs**(-1)*tDg**(-1)
     &  * ( half*za(jb,jn)*za(ig,ie)*zb(is,ia)*zb(in,ia)*zb(ia,je)
     &     )
      gs(2,1) = gs(2,1) + 1/(zb(ig,ia))*sga2**(-1)*tga2*aDs**(-1) * ( 
     &    za(jb,jn)*za(ig,is)*za(ie,jn)*zb(is,ia)**2*zb(in,ia)*zb(jn,je
     &    ) + za(jb,jn)*za(ig,is)*za(ie,jb)*zb(jb,je)*zb(is,ia)**2*zb(
     &    in,ia) )
      gs(2,1) = gs(2,1) + 1/(zb(ig,ia))*tga2*aDs**(-1)*tDg**(-1) * (
     &    half*za(jb,jn)*za(ig,jn)*za(ie,jn)*zb(is,ia)*zb(in,ia)*zb(
     &    ia,jn)*zb(jn,je) + half*za(jb,jn)*za(ig,jn)*za(ie,je)*
     &    zb(is,ia)*zb(in,ia)*zb(ia,je)*zb(jn,je) - half*za(jb,jn)
     &    *za(ig,jn)*za(ie,jb)*zb(jb,ia)*zb(is,ia)*zb(in,ia)*zb(jn,je)
     &     + half*za(jb,jn)*za(ig,jb)*za(ie,jn)*zb(jb,je)*zb(is,ia
     &    )*zb(in,ia)*zb(ia,jn) + half*za(jb,jn)*za(ig,jb)*za(ie,
     &    je)*zb(jb,je)*zb(is,ia)*zb(in,ia)*zb(ia,je) - half*za(jb
     &    ,jn)*za(ig,jb)*za(ie,jb)*zb(jb,ia)*zb(jb,je)*zb(is,ia)*zb(in,
     &    ia) )

      gs(2,2)=  + sga2**(-1)*tga2*aDs**(-1) * ( za(jb,jn)*za(ie,jn)*zb(
     &    ig,in)*zb(ig,ia)*zb(is,ia)*zb(jn,je) + za(jb,jn)*za(ie,jb)*
     &    zb(jb,je)*zb(ig,in)*zb(ig,ia)*zb(is,ia) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))*msq*tga2*aDs**(-1)*tDg**(-1)
     &  * ( half*za(jb,jn)*za(ie,ia)*zb(ig,je)*zb(is,ia)*zb(in,ia)
     &     )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))*msq*tDg**(-1) * (  - za(jb,jn)*
     &    za(ie,jn)*zb(ig,je)*zb(is,in)*zb(ia,jn) - za(jb,jn)*za(ie,je)
     &    *zb(ig,je)*zb(is,in)*zb(ia,je) + za(jb,jn)*za(ie,jb)*zb(jb,ia
     &    )*zb(ig,je)*zb(is,in) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))*sga2**(-1)*tga2*aDs**(-1) * ( 
     &    za(jb,jn)*za(ig,is)*za(ie,jn)*zb(ig,is)*zb(ig,in)*zb(is,ia)*
     &    zb(jn,je) + za(jb,jn)*za(ig,is)*za(ie,jb)*zb(jb,je)*zb(ig,is)
     &    *zb(ig,in)*zb(is,ia) + za(jb,jn)*za(is,ia)*za(ie,jn)*zb(ig,ia
     &    )*zb(is,in)*zb(is,ia)*zb(jn,je) + za(jb,jn)*za(is,ia)*za(ie,
     &    jb)*zb(jb,je)*zb(ig,ia)*zb(is,in)*zb(is,ia) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))*sga2**(-1)*tga2 * (  - za(jb,jn
     &    )*za(ie,jn)*zb(ig,is)*zb(in,ia)*zb(jn,je) - za(jb,jn)*za(ie,
     &    jn)*zb(ig,in)*zb(is,ia)*zb(jn,je) - za(jb,jn)*za(ie,jb)*zb(jb
     &    ,je)*zb(ig,is)*zb(in,ia) - za(jb,jn)*za(ie,jb)*zb(jb,je)*zb(
     &    ig,in)*zb(is,ia) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))*tga2*aDs**(-1)*tDg**(-1) * ( 
     &    half*za(jb,ia)*za(jb,jn)*za(ie,jn)*zb(jb,je)*zb(ig,jn)*zb(
     &    is,ia)*zb(in,ia) + half*za(jb,ia)*za(jb,jn)*za(ie,je)*
     &    zb(jb,je)*zb(ig,je)*zb(is,ia)*zb(in,ia) + half*za(jb,ia)
     &    *za(jb,jn)*za(ie,jb)*zb(jb,je)*zb(ig,jb)*zb(is,ia)*zb(in,ia)
     &     - half*za(jb,jn)*za(ie,jn)*za(ia,jn)*zb(ig,jn)*zb(is,ia
     &    )*zb(in,ia)*zb(jn,je) - half*za(jb,jn)*za(ie,je)*za(ia,
     &    jn)*zb(ig,je)*zb(is,ia)*zb(in,ia)*zb(jn,je) - half*za(jb
     &    ,jn)*za(ie,jb)*za(ia,jn)*zb(ig,jb)*zb(is,ia)*zb(in,ia)*zb(jn,
     &    je) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))/(za(ig,ia))*msq*tga2*aDs**(-1)*
     & tDg**(-1) * (  - half*za(jb,jn)*za(ig,is)*za(ie,ia)*zb(ig,
     &    je)*zb(is,in)*zb(is,ia) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))/(za(ig,ia))*msq*tDg**(-1) * ( 
     &    za(jb,ia)*za(jb,jn)*za(ig,ie)*zb(jb,ia)*zb(ig,je)*zb(is,in)
     &     - za(jb,ia)*za(jb,jn)*za(ig,ie)*zb(jb,je)*zb(ig,ia)*zb(is,in
     &    ) + za(jb,jn)*za(ig,ie)*za(ia,jn)*zb(ig,ia)*zb(is,in)*zb(jn,
     &    je) + za(jb,jn)*za(ig,ie)*za(ia,jn)*zb(ig,je)*zb(is,in)*zb(ia
     &    ,jn) + za(jb,jn)*za(ig,ie)*za(ia,je)*zb(ig,je)*zb(is,in)*zb(
     &    ia,je) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))/(za(ig,ia))*msq * (  - za(jb,jn
     &    )*za(ig,ie)*zb(ig,je)*zb(is,in) - za(jb,jn)*za(ie,ia)*zb(is,
     &    in)*zb(ia,je) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))/(za(ig,ia))*sga2**(-1)*tga2*
     & aDs**(-1) * ( za(jb,jn)*za(ig,is)*za(is,ia)*za(ie,jn)*zb(ig,is)*
     &    zb(is,in)*zb(is,ia)*zb(jn,je) + za(jb,jn)*za(ig,is)*za(is,ia)
     &    *za(ie,jb)*zb(jb,je)*zb(ig,is)*zb(is,in)*zb(is,ia) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))/(za(ig,ia))*sga2**(-1)*tga2
     &  * ( za(jb,jn)*za(ig,is)*za(ie,jn)*zb(ig,is)*zb(is,in)*zb(jn,je)
     &     + za(jb,jn)*za(ig,is)*za(ie,jb)*zb(jb,je)*zb(ig,is)*zb(is,in
     &    ) - za(jb,jn)*za(is,ia)*za(ie,jn)*zb(is,in)*zb(is,ia)*zb(jn,
     &    je) - za(jb,jn)*za(is,ia)*za(ie,jb)*zb(jb,je)*zb(is,in)*zb(is
     &    ,ia) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))/(za(ig,ia))*tga2*aDs**(-1)*
     & tDg**(-1) * (  - half*za(jb,ia)*za(jb,jn)*za(ig,is)*za(ie,
     &    jn)*zb(jb,je)*zb(ig,jn)*zb(is,in)*zb(is,ia) - half*za(jb
     &    ,ia)*za(jb,jn)*za(ig,is)*za(ie,je)*zb(jb,je)*zb(ig,je)*zb(is,
     &    in)*zb(is,ia) - half*za(jb,ia)*za(jb,jn)*za(ig,is)*za(ie
     &    ,jb)*zb(jb,je)*zb(ig,jb)*zb(is,in)*zb(is,ia) + half*za(
     &    jb,jn)*za(ig,is)*za(ie,jn)*za(ia,jn)*zb(ig,jn)*zb(is,in)*zb(
     &    is,ia)*zb(jn,je) + half*za(jb,jn)*za(ig,is)*za(ie,je)*
     &    za(ia,jn)*zb(ig,je)*zb(is,in)*zb(is,ia)*zb(jn,je) + half
     &    *za(jb,jn)*za(ig,is)*za(ie,jb)*za(ia,jn)*zb(ig,jb)*zb(is,in)*
     &    zb(is,ia)*zb(jn,je) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))/(za(ig,ia))*tDg**(-1) * (  - 
     &    za(jb,ia)*za(jb,jn)*za(ig,ie)*za(ig,jn)*zb(jb,je)*zb(ig,ia)*
     &    zb(ig,jn)*zb(is,in) - za(jb,ia)*za(jb,jn)*za(ig,ie)*za(ig,je)
     &    *zb(jb,je)*zb(ig,ia)*zb(ig,je)*zb(is,in) - za(jb,ia)*za(jb,jn
     &    )*za(ig,jn)*za(ie,jn)*zb(jb,je)*zb(ig,jn)*zb(is,in)*zb(ia,jn)
     &     - za(jb,ia)*za(jb,jn)*za(ig,jn)*za(ie,je)*zb(jb,je)*zb(ig,jn
     &    )*zb(is,in)*zb(ia,je) + za(jb,ia)*za(jb,jn)*za(ig,jn)*za(ie,
     &    jb)*zb(jb,ia)*zb(jb,je)*zb(ig,jn)*zb(is,in) - za(jb,ia)*za(jb
     &    ,jn)*za(ig,je)*za(ie,jn)*zb(jb,je)*zb(ig,je)*zb(is,in)*zb(ia,
     &    jn) - za(jb,ia)*za(jb,jn)*za(ig,je)*za(ie,je)*zb(jb,je)*zb(ig
     &    ,je)*zb(is,in)*zb(ia,je) + za(jb,ia)*za(jb,jn)*za(ig,je)*za(
     &    ie,jb)*zb(jb,ia)*zb(jb,je)*zb(ig,je)*zb(is,in) - za(jb,ia)*
     &    za(jb,jn)*za(ig,jb)*za(ig,ie)*zb(jb,je)*zb(ig,jb)*zb(ig,ia)*
     &    zb(is,in) - za(jb,ia)*za(jb,jn)*za(ig,jb)*za(ie,jn)*zb(jb,je)
     &    *zb(ig,jb)*zb(is,in)*zb(ia,jn) - za(jb,ia)*za(jb,jn)*za(ig,jb
     &    )*za(ie,je)*zb(jb,je)*zb(ig,jb)*zb(is,in)*zb(ia,je) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))/(za(ig,ia))*tDg**(-1) * ( za(jb
     &    ,ia)*za(jb,jn)*za(ig,jb)*za(ie,jb)*zb(jb,ia)*zb(jb,je)*zb(ig,
     &    jb)*zb(is,in) + za(jb,jn)*za(ig,ie)*za(ig,jn)*za(ia,jn)*zb(ig
     &    ,ia)*zb(ig,jn)*zb(is,in)*zb(jn,je) + za(jb,jn)*za(ig,ie)*za(
     &    ig,je)*za(ia,jn)*zb(ig,ia)*zb(ig,je)*zb(is,in)*zb(jn,je) + 
     &    za(jb,jn)*za(ig,jn)*za(ie,jn)*za(ia,jn)*zb(ig,jn)*zb(is,in)*
     &    zb(ia,jn)*zb(jn,je) + za(jb,jn)*za(ig,jn)*za(ie,je)*za(ia,jn)
     &    *zb(ig,jn)*zb(is,in)*zb(ia,je)*zb(jn,je) - za(jb,jn)*za(ig,jn
     &    )*za(ie,jb)*za(ia,jn)*zb(jb,ia)*zb(ig,jn)*zb(is,in)*zb(jn,je)
     &     + za(jb,jn)*za(ig,je)*za(ie,jn)*za(ia,jn)*zb(ig,je)*zb(is,in
     &    )*zb(ia,jn)*zb(jn,je) + za(jb,jn)*za(ig,je)*za(ie,je)*za(ia,
     &    jn)*zb(ig,je)*zb(is,in)*zb(ia,je)*zb(jn,je) - za(jb,jn)*za(ig
     &    ,je)*za(ie,jb)*za(ia,jn)*zb(jb,ia)*zb(ig,je)*zb(is,in)*zb(jn,
     &    je) + za(jb,jn)*za(ig,jb)*za(ig,ie)*za(ia,jn)*zb(ig,jb)*zb(ig
     &    ,ia)*zb(is,in)*zb(jn,je) + za(jb,jn)*za(ig,jb)*za(ie,jn)*za(
     &    ia,jn)*zb(ig,jb)*zb(is,in)*zb(ia,jn)*zb(jn,je) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))/(za(ig,ia))*tDg**(-1) * ( za(jb
     &    ,jn)*za(ig,jb)*za(ie,je)*za(ia,jn)*zb(ig,jb)*zb(is,in)*zb(ia,
     &    je)*zb(jn,je) - za(jb,jn)*za(ig,jb)*za(ie,jb)*za(ia,jn)*zb(jb
     &    ,ia)*zb(ig,jb)*zb(is,in)*zb(jn,je) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))/(za(ig,ia)) * (  - za(jb,ia)*
     &    za(jb,jn)*za(ig,ie)*zb(jb,je)*zb(ig,ia)*zb(is,in) - za(jb,ia)
     &    *za(jb,jn)*za(ie,jn)*zb(jb,je)*zb(is,in)*zb(ia,jn) - za(jb,ia
     &    )*za(jb,jn)*za(ie,je)*zb(jb,je)*zb(is,in)*zb(ia,je) + za(jb,
     &    ia)*za(jb,jn)*za(ie,jb)*zb(jb,ia)*zb(jb,je)*zb(is,in) + za(jb
     &    ,jn)*za(ig,ie)*za(ia,jn)*zb(ig,ia)*zb(is,in)*zb(jn,je) - za(
     &    jb,jn)*za(ig,jn)*za(ie,ia)*zb(ig,ia)*zb(is,in)*zb(jn,je) - 
     &    za(jb,jn)*za(ig,jn)*za(ie,jn)*zb(ig,jn)*zb(is,in)*zb(jn,je)
     &     - za(jb,jn)*za(ig,jn)*za(ie,je)*zb(ig,je)*zb(is,in)*zb(jn,je
     &    ) - za(jb,jn)*za(ig,jn)*za(ie,jb)*zb(ig,jb)*zb(is,in)*zb(jn,
     &    je) - za(jb,jn)*za(ig,jb)*za(ie,ia)*zb(jb,je)*zb(ig,ia)*zb(is
     &    ,in) - za(jb,jn)*za(ig,jb)*za(ie,jn)*zb(jb,je)*zb(ig,jn)*zb(
     &    is,in) - za(jb,jn)*za(ig,jb)*za(ie,je)*zb(jb,je)*zb(ig,je)*
     &    zb(is,in) - za(jb,jn)*za(ig,jb)*za(ie,jb)*zb(jb,je)*zb(ig,jb)
     &    *zb(is,in) + za(jb,jn)*za(ie,jn)*za(ia,jn)*zb(is,in)*zb(ia,jn
     &    )*zb(jn,je) )
      gs(2,2) = gs(2,2) + 1/(za(ig,ia))/(za(ig,ia)) * ( za(jb,jn)*za(ie
     &    ,je)*za(ia,jn)*zb(is,in)*zb(ia,je)*zb(jn,je) - za(jb,jn)*za(
     &    ie,jb)*za(ia,jn)*zb(jb,ia)*zb(is,in)*zb(jn,je) )

c--- Use the "overall scheme" to include the top width when necessary
      if (tga2+mq**2 > zero) then
        tga2=sqrt(tga2**2+(mq*qwidth)**2)
c      if (real(tga2+mt**2) > zero) then
c        tga2=tga2+cplx2(zero,mt*twidth)
      endif
      
      gs(1,1)=gs(1,1)/tga2
      gs(1,2)=gs(1,2)/tga2
      gs(2,1)=gs(2,1)/tga2
      gs(2,2)=gs(2,2)/tga2

      return
      end
