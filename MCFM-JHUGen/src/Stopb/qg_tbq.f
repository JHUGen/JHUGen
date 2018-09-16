      subroutine qg_tbq(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Lowest order t-channel single top, with explicit b-quark         *
*                                                                      *
*     q(p1) + g(p2) -> t(p3) + b(p4) + q'(p5)                          *
*                                                                      *
*     Originally: R. Frederix and F. Tramontano, February 2008         *
*        Adapted: J. Campbell, February 26, 2008                       *
*                                                                      *
************************************************************************
c     u + g  ->  c + s + d  (t-channel single-charm)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'couple.f'
      include 'masses.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'stopscales.f'
      include 'dynamicscale.f'
      include 'nlooprun.f'
      include 'qcdcouple.f'
      integer:: i3,i4
      real(dp):: p(mxpart,4),fac1,fac2,msq_qg,msq_gq,
     & msq_qbarg,msq_gqbar
      real(dp):: dot,msq(-nf:nf,-nf:nf),
     & Wprop15,Wprop25,mq,ma,gsq_H
      real(dp):: sDg,sDu,sDc,sDd,gDu,gDc,gDd,cDd,uDc,lord
      complex(dp)::amps_qg(2,2,2),amps_qbarg(2,2,2),amps_gq(2,2,2),
     & amps_gqbar(2,2,2)
      logical:: nocheck
      integer:: hg,hc,hs,j,k
      real(dp):: alphas
      real(dp):: b1scale,q2scale,q1scale,b2scale
      common/bqscale/b1scale,q2scale,q1scale,b2scale
!$omp threadprivate(/bqscale/)
c--- set this parameter to .false. to check amplitudes vs. squared ME
      parameter(nocheck=.true.)

c--- initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      gsq_H=fourpi*as_H

      if (dynstring == 'DDIS') then
        fac1=aveqg*2._dp*xn**2*Cf*gwsq**2
     &       *fourpi*alphas(b1scale,amz,nlooprun)
        fac2=aveqg*2._dp*xn**2*Cf*gwsq**2
     &       *fourpi*alphas(b2scale,amz,nlooprun)
      else
        fac1=aveqg*2._dp*xn**2*Cf*gsq_H*gwsq**2
        fac2=fac1
      endif


c--- extra correction factor to go from 4 to 5 flavours
c      fac=fac*(1._dp-as_H/2._dp/pi*2._dp/3._dp*log(renscale_H/facscale_H))

c--- set mass of quark and antiquark according to nwz
      if (nwz == +1) then
        mq=mt
        ma=mb
        i3=3
        i4=4
      else
        mq=mb
        ma=mt
        i3=4
        i4=3
      endif

c--- propagator for qg and qbarg amplitudes
      Wprop15=1._dp/(2._dp*dot(p,1,5)-wmass**2)

      call bornwrap(p,1,2,i3,i4,5,mq,ma,amps_qg)
      call bornwrap(p,5,2,i3,i4,1,mq,ma,amps_qbarg)

c--- propagator for gq and gqbar amplitudes
      Wprop25=1._dp/(2._dp*dot(p,2,5)-wmass**2)

      call bornwrap(p,2,1,i3,i4,5,mq,ma,amps_gq)
      call bornwrap(p,5,1,i3,i4,2,mq,ma,amps_gqbar)

c--- now square up amplitudes and add propagators
      msq_qg=0._dp
      msq_qbarg=0._dp
      msq_gq=0._dp
      msq_gqbar=0._dp
      do hg=1,2
      do hc=1,2
      do hs=1,2
      msq_qg=msq_qg+Wprop15**2*abs(amps_qg(hg,hc,hs))**2
      msq_qbarg=msq_qbarg+Wprop15**2*abs(amps_qbarg(hg,hc,hs))**2
      msq_gq=msq_gq+Wprop25**2*abs(amps_gq(hg,hc,hs))**2
      msq_gqbar=msq_gqbar+Wprop25**2*abs(amps_gqbar(hg,hc,hs))**2
      enddo
      enddo
      enddo

c--- fill matrix elements
      do j=1,4
        msq(+j,0)=fac2*Vsum(+j)*msq_qg
        msq(-j,0)=fac2*Vsum(-j)*msq_qbarg
        msq(0,+j)=fac1*Vsum(+j)*msq_gq
        msq(0,-j)=fac1*Vsum(-j)*msq_gqbar
      enddo

c--- normal return is here
      if (nocheck) return

c--- else, check against squared matrix elements (code below)

c     u + g  ->  c + s + d  (t-channel single-charm)

      sDg=dot(p,4,2)
      sDu=dot(p,4,1)
      sDc=dot(p,4,3)
      sDd=dot(p,4,5)
      gDu=dot(p,1,2)
      gDc=dot(p,2,3)
      gDd=dot(p,2,5)
      cDd=dot(p,3,5)
      uDc=dot(p,1,3)
c      gDs=dot(p,4,2)
c      cDs=dot(p,4,3)
c      cDg=dot(p,2,3)
c      uDd=dot(p,1,5)
      lord =
     &  - 2._dp*sDg**(-2)*sDu*cDd*mb**2
     &  - 2._dp*sDg**(-2)*gDu*cDd*mb**2
     &  + 2._dp*sDg**(-1)*sDu*sDc*gDc**(-1)*gDd
     &  + 4._dp*sDg**(-1)*sDu*sDc*gDc**(-1)*cDd
     &  - 2._dp*sDg**(-1)*sDu*sDd
     &  + 2._dp*sDg**(-1)*sDu*cDd
     &  + 2._dp*sDg**(-1)*sDc*gDu*gDc**(-1)*cDd
     &  + 2._dp*sDg**(-1)*gDu*cDd
     &  - 2._dp*sDu*gDc**(-2)*gDd*mt**2
     &  - 2._dp*sDu*gDc**(-2)*cDd*mt**2
     &  + 2._dp*sDu*gDc**(-1)*gDd
     &  + 2._dp*sDu*gDc**(-1)*cDd
     &  - 2._dp*gDc**(-1)*cDd*uDc
     &
      write(6,*)
c      write(6,*) 'amps(2,2,2)=',amps(2,2,2)
c      write(6,*) 'amps(2,2,1)=',amps(2,2,1)
c      write(6,*) 'amps(2,1,2)=',amps(2,1,2)
c      write(6,*) 'amps(2,1,1)=',amps(2,1,1)
c      write(6,*) 'amps(1,2,2)=',amps(1,2,2)
c      write(6,*) 'amps(1,2,1)=',amps(1,2,1)
c      write(6,*) 'amps(1,1,2)=',amps(1,1,2)
c      write(6,*) 'amps(1,1,1)=',amps(1,1,1)
      write(6,*) 'qg = ',-lord*Wprop15**2,msq_qg,-lord*Wprop15**2/msq_qg

c     g + u  ->  c + s + d  (t-channel single-charm)

      sDg=dot(p,4,1)
      sDu=dot(p,4,2)
      sDc=dot(p,4,3)
      sDd=dot(p,4,5)
      gDu=dot(p,2,1)
      gDc=dot(p,1,3)

      gDd=dot(p,1,5)
      cDd=dot(p,3,5)
      uDc=dot(p,2,3)
c      gDs=dot(p,4,1)
c      cDs=dot(p,4,3)
c      cDg=dot(p,1,3)
c      uDd=dot(p,2,5)
      lord =
     &  - 2._dp*sDg**(-2)*sDu*cDd*mb**2
     &  - 2._dp*sDg**(-2)*gDu*cDd*mb**2
     &  + 2._dp*sDg**(-1)*sDu*sDc*gDc**(-1)*gDd
     &  + 4._dp*sDg**(-1)*sDu*sDc*gDc**(-1)*cDd
     &  - 2._dp*sDg**(-1)*sDu*sDd
     &  + 2._dp*sDg**(-1)*sDu*cDd
     &  + 2._dp*sDg**(-1)*sDc*gDu*gDc**(-1)*cDd
     &  + 2._dp*sDg**(-1)*gDu*cDd
     &  - 2._dp*sDu*gDc**(-2)*gDd*mt**2
     &  - 2._dp*sDu*gDc**(-2)*cDd*mt**2
     &  + 2._dp*sDu*gDc**(-1)*gDd
     &  + 2._dp*sDu*gDc**(-1)*cDd
     &  - 2._dp*gDc**(-1)*cDd*uDc
     &
      write(6,*) 'gq = ',-lord*Wprop25**2,msq_gq,-lord*Wprop25**2/msq_gq
c      pause

      return
      end


c--- wrapper to LO amplitude routine bornampsN that allows the
c---  momenta to be permuted according to i1,i2,i5
      subroutine bornwrap(p,i1,i2,i3,i4,i5,mh,ml,amps)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      integer:: i1,i2,i3,i4,i5,j
      real(dp):: p(mxpart,4),q(mxpart,4),mh,ml,dot
      complex(dp):: amps(2,2,2)

      do j=1,4
        q(1,j)=p(i1,j)
        q(2,j)=p(i2,j)
        q(3,j)=p(i3,j)-mh**2/2._dp/dot(p,i2,i3)*p(i2,j)
        q(4,j)=p(i4,j)-ml**2/2._dp/dot(p,i2,i4)*p(i2,j)
        q(5,j)=p(i5,j)
      enddo

c--- set up spinor products
      call spinoru(5,q,za,zb)

c---- calling amps(hg,hc,hs)
      call bornampsN(q,mh,ml,amps)

      return
      end


      subroutine bornampsN(q,mc,ms,amps)
      implicit none
      include 'types.f'
c     u + g  ->  c + s + d  (t-channel single-charm)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      real(dp):: q(mxpart,4),dot,cDs,gDs,cDg,mc,ms
      complex(dp):: trg,trs,trc,trsgc,amps(2,2,2)

      cDg=dot(q,3,2)
      gDs=dot(q,4,2)
      cDs=dot(q,3,4)+mc**2*gDs/2._dp/cDg
     &              +ms**2*cDg/2._dp/gDs

      trg=2._dp*za(5,2)*zb(2,1)
      trs=2._dp*za(5,4)*zb(4,1)+ms**2*za(5,2)*zb(2,1)/gDs
      trc=2._dp*za(5,3)*zb(3,1)+mc**2*za(5,2)*zb(2,1)/cDg
      trsgc=2._dp*zb(1,4)*za(4,2)*zb(2,3)*za(3,5)


      amps(2,2,2)=(mc*(cDg*(-ms**2*trg+2._dp*gDs*(trg + trs))-
     -      gDs*trsgc))/(4._dp*cDg*gDs)

      amps(2,2,1)=-(mc*ms*(-2._dp*cDg*cDs*gDs+gDs**2*mc**2+
     -       cDg**2*ms**2)*trg)/(4._dp*cDg*gDs**2)

      amps(2,1,2)=(gDs**2*mc**2*trsgc-2._dp*cDg*gDs*(gDs*mc**2*trg
     -  +cDs*trsgc)+cDg**2*(4._dp*gDs**2*trc-2._dp*gDs*trsgc +
     -       ms**2*trsgc))/(4._dp*cDg**2*gDs)

      amps(2,1,1)=-(ms*(-2._dp*cDg*cDs*gDs+gDs**2*mc**2+
     -   cDg**2*ms**2)*(2._dp*cDg*gDs*trc-gDs*mc**2*trg-cDg*trsgc))/
     -  (4._dp*cDg**2*gDs**2)

      amps(1,2,2)=-(mc*(-2._dp*cDg*cDs*gDs+gDs**2*mc**2+cDg**2*ms**2)*
     -     (-(cDg*ms**2*trg) + 2._dp*cDg*gDs*trs - gDs*trsgc)
     -     )/(4._dp*cDg**2*gDs**2)

      amps(1,2,1)=(mc*ms*(-2._dp*cDg*cDs*gDs + gDs**2*mc**2 +
     -      cDg**2*ms**2)*trg)/(4._dp*cDg**2*gDs)

      amps(1,1,2)=-(-2._dp*cDg*gDs*(cDs+gDs)*trsgc+gDs**2*mc**2*trsgc+
     -     cDg**2*(-2._dp*gDs*ms**2*trg + 4._dp*gDs**2*trs +
     -        ms**2*trsgc))/(4._dp*cDg*gDs**2)

      amps(1,1,1)=(ms*(-gDs*mc**2*trg+cDg*(2._dp*gDs*(trc+trg)-trsgc)))/
     -  (4._dp*cDg*gDs)


      amps(2,2,2)=amps(2,2,2)/za(2,4)/za(2,3)
      amps(2,2,1)=amps(2,2,1)/za(2,3)**2/zb(4,3)
      amps(2,1,2)=amps(2,1,2)/za(2,4)**2/zb(4,3)
      amps(2,1,1)=amps(2,1,1)/za(2,4)/za(2,3)/zb(4,3)**2
      amps(1,2,2)=amps(1,2,2)/za(4,3)**2/zb(2,4)/zb(2,3)
      amps(1,2,1)=amps(1,2,1)/za(4,3)/zb(2,4)**2
      amps(1,1,2)=amps(1,1,2)/za(4,3)/zb(2,3)**2
      amps(1,1,1)=amps(1,1,1)/zb(2,4)/zb(2,3)

      return
      end

