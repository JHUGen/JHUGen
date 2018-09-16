      subroutine gs_wc_dg(p,ig,is,ie,in,jn,je,jb,ja,jt,amp)
      implicit none
      include 'types.f'
c     g(ig)+s(is)-->W^-{e^-(ie)+nbar(in)}+t{W^+[n(jn)+e^+(je)]+b(jb)}+a(ja)
c---- helicities: gs(ht,ha)
c--- labels on amplitudes represent helicities for the heavy quark
c--- and the gluon (ja) respectively
c---   1 = negative helicity, 2 = positive helicity
c--- heavy quark momentum is made massless (ic) with the gluon momentum ig

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: is,ig,ie,in,ja,je,jn,jb,jt,i,j
      real(dp):: p(mxpart,4)
      real(dp):: dot,mtsq,tDa,tDg,bDa,tsq,propd,propt
      complex(dp):: amp(2,2)

      propd=sqrt((s(je,jn)-wmass**2)**2+(wmass*wwidth)**2)
      tsq  =s(je,jn)+s(je,jb)+s(je,ja)+s(jn,jb)+s(jn,ja)+s(jb,ja)
      propt=sqrt((tsq-mt**2)**2+(mt*twidth)**2)

      mtsq=mt**2
      bDa=+dot(p,jb,ja)
      tDa=+dot(p,je,ja)+dot(p,jn,ja)+dot(p,jb,ja)
      tDg=-dot(p,is,ig)-dot(p,ie,ig)-dot(p,in,ig)

      amp(1,1)=  + mtsq*tDg**(-1)*tDa**(-1) * ( one/four*za(jb,jn)*za(
     &    ig,ja)*zb(ig,je)*zb(is,jt)/zb(is,ja) )
      amp(1,1) = amp(1,1) + mtsq*tDa**(-1) * ( half*za(jb,jn)*za(
     &    ig,ja)*zb(is,je)/za(ig,jt)/zb(is,ja) )
      amp(1,1) = amp(1,1) + tDa**(-1) * ( half*za(jb,jn)*za(ja,jt)
     &    *zb(is,jt)*zb(je,jt)/zb(is,ja) )
      amp(1,1) = amp(1,1) + bDa**(-1) * ( half*za(jb,ja)*za(ja,jn)
     &    *zb(je,jt) + half*za(jb,ja)*za(jb,jn)*zb(is,jb)*zb(je,jt
     &    )/zb(is,ja) )

      amp(1,2)=  + mtsq*tDg**(-1)*tDa**(-1) * (  - one/four*za(jb,jn)*
     &    za(ig,is)*zb(ig,je)*zb(ja,jt)/za(is,ja) )
      amp(1,2) = amp(1,2) + mtsq*tDa**(-1) * (  - half*za(jb,jn)*
     &    za(ig,is)*zb(ja,je)/za(ig,jt)/za(is,ja) )
      amp(1,2) = amp(1,2) + tDa**(-1) * (  - half*za(jb,jn)*za(is,
     &    jt)*zb(ja,jt)*zb(je,jt)/za(is,ja) - half*za(jb,jn)*zb(
     &    ja,jt)*zb(ja,je) )
      amp(1,2) = amp(1,2) + bDa**(-1) * (  - half*za(jb,jn)*za(is,
     &    jb)*zb(jb,ja)*zb(je,jt)/za(is,ja) )

      amp(2,1)=  + mt*mtsq*tDg**(-1)*tDa**(-1) * ( one/four*za(jb,jn)*
     &    za(ig,ja)*zb(ig,is)*zb(ig,je)/zb(ig,jt)/zb(is,ja) )
      amp(2,1) = amp(2,1) + mt*tDa**(-1) * ( half*za(jb,jn)*za(ja,
     &    jt)*zb(ig,is)*zb(je,jt)/zb(ig,jt)/zb(is,ja) + half*za(
     &    jb,jn)*za(ja,jt)*zb(is,je)/zb(is,ja) )
      amp(2,1) = amp(2,1) + mt*bDa**(-1) * ( half*za(jb,ja)*za(ja,
     &    jn)*zb(ig,je)/zb(ig,jt) + half*za(jb,ja)*za(jb,jn)*zb(
     &    ig,je)*zb(is,jb)/zb(ig,jt)/zb(is,ja) )

      amp(2,2)=  + mt*mtsq*tDg**(-1)*tDa**(-1) * (  - one/four*za(jb,
     &    jn)*za(ig,is)*zb(ig,ja)*zb(ig,je)/za(is,ja)/zb(ig,jt) )
      amp(2,2) = amp(2,2) + mt*tDa**(-1) * (  - half*za(jb,jn)*za(
     &    is,jt)*zb(ig,ja)*zb(je,jt)/za(is,ja)/zb(ig,jt) - half*
     &    za(jb,jn)*za(is,jt)*zb(ja,je)/za(is,ja) - half*za(jb,jn
     &    )*zb(ig,ja)*zb(ja,je)/zb(ig,jt) )
      amp(2,2) = amp(2,2) + mt*bDa**(-1) * (  - half*za(jb,jn)*za(
     &    is,jb)*zb(jb,ja)*zb(ig,je)/za(is,ja)/zb(ig,jt) )

      do i=1,2
      do j=1,2
      amp(i,j)=amp(i,j)/propd/propt
      enddo
      enddo

      return
      end
