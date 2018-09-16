      subroutine qqb_wgam_g(p,msq)
      implicit none
      include 'types.f'
      
C-----Author Keith Ellis, September 2002
C----- updated: John Campbell, August 2011 (anomalous couplings)
c----Matrix element for W gam production
C----averaged over initial colours and spins
C For nwz=+1
c     u(-p1)+dbar(-p2)-->W^+(n(p3)+e^+(p4)) + gamma(p5) + g(p6)
C For nwz=-1
c     ubar(-p1)+d(-p2)-->W^-(e^-(p3)+nbar(p4)) + gamma(p5) + g(p6)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'zprods_com.f'
      include 'nwz.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      real(dp):: qbq,qqb,qg,gq,qbg,gqb
      real(dp):: ubdgmsq


      call spinoru(6,p,za,zb)
      fac=V/2._dp*gwsq**2*4._dp*gsq*esq

      if (nwz == -1) then 
      qbq=aveqq*fac*ubdgmsq(1,2,3,4,5,6,za,zb)
      qqb=aveqq*fac*ubdgmsq(2,1,3,4,5,6,za,zb)
      gq=aveqg*fac*ubdgmsq(6,2,3,4,5,1,za,zb)
      qg=aveqg*fac*ubdgmsq(6,1,3,4,5,2,za,zb)
      gqb=aveqg*fac*ubdgmsq(2,6,3,4,5,1,za,zb)
      qbg=aveqg*fac*ubdgmsq(1,6,3,4,5,2,za,zb)

      elseif (nwz == +1) then 
      qbq=aveqq*fac*ubdgmsq(2,1,4,3,5,6,zb,za)
      qqb=aveqq*fac*ubdgmsq(1,2,4,3,5,6,zb,za)
      gq=aveqg*fac*ubdgmsq(2,6,4,3,5,1,zb,za)
      qg=aveqg*fac*ubdgmsq(1,6,4,3,5,2,zb,za)
      gqb=aveqg*fac*ubdgmsq(6,2,4,3,5,1,zb,za)
      qbg=aveqg*fac*ubdgmsq(6,1,4,3,5,2,zb,za)

      endif
     

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
          if ((j > 0) .and. (k < 0)) then
            msq(j,k)=Vsq(j,k)*qqb
          elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=Vsum(k)*gqb
            
          elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=Vsum(k)*gq
          elseif ((j < 0) .and. (k > 0)) then
            msq(j,k)=Vsq(j,k)*qbq
          elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=Vsum(j)*qg
         elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=Vsum(j)*qbg
          endif
      enddo
      enddo
      return
      end

      function ubdgmsq(p1,p2,p3,p4,p5,p6,za,zb)
      implicit none
      include 'types.f'
      real(dp):: ubdgmsq
      
C     Matrix element for 
C     ub(-p1)+d(-p2)=e-(p3)+nu~(p4)+gamma(p5)+g(p6)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'zerowidth.f'
      include 'anomcoup.f'
      include 'xanomcoup.f'
      integer:: p1,p2,p3,p4,p5,p6
      complex(dp):: aLL,aRR,aRL,aLR,prp34,prp345,zazb
      real(dp):: s345,s156,s256,xfac

      zazb(p1,p2,p3,p4)=+za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s156=s(p1,p5)+s(p1,p6)+s(p5,p6)
      s256=s(p2,p5)+s(p2,p6)+s(p5,p6)
      s345=s(p3,p4)+s(p3,p5)+s(p4,p5)
      prp34=s(p3,p4)/cplx2(s(p3,p4)-wmass**2,wmass*wwidth)

    
c--- apply a dipole form factor to anomalous couplings, with power two (only if tevscale > 0)
      if (tevscale > 0._dp) then
        xfac=1._dp/(1._dp+s345/(tevscale*1d3)**2)**2
      else
        xfac=1._dp
      endif
      xdelk_g=xfac*delk_g
      xlambda_g=xfac*lambda_g

      if (zerowidth) then
c--- zerowidth: no final state radiation, so we can set prp345 to zero
        prp345=czip
      else
c--- otherwise, usual Breit-Wigner form
        prp345=s345/cplx2(s345-wmass**2,wmass*wwidth)
      endif

c---  c.f. Eqs.(4.9)-(4.12) of hep-ph/9803250 (multiplied by -i)
c---       for the terms proportional to prp34
      aRR=za(p1,p3)**2/(za(p1,p6)*za(p2,p6)*(s345-s(p3,p4)))*(
     &  Qd*(+zazb(p2,p3,p4,p5)/za(p4,p3)*prp34
     &      +za(p2,p5)*zb(p4,p5)/za(p3,p5)*prp345)/za(p2,p5)
     & +Qu*(+zazb(p1,p3,p4,p5)/za(p3,p4)*prp34
     &      -za(p1,p5)*zb(p4,p5)/za(p3,p5)*prp345)/za(p1,p5))
      
      aLR=Qd*(
     & (-za(p1,p3)*zb(p6,p2)*zazb(p5,p1,p3,p4)
     &  /(zb(p2,p5)*za(p6,p2)*s256)
     &  -zazb(p1,p2,p6,p4)*(za(p3,p4)*za(p1,p5)*zb(p4,p2)
     &                     +za(p3,p5)*za(p1,p6)*zb(p6,p2))
     &  /(zb(p2,p5)*za(p1,p6)*za(p6,p2)*(s345-s(p3,p4))))*prp34/s(p3,p4)
     & +zazb(p1,p2,p6,p4)**2*za(p4,p5)
     &   /(zb(p3,p5)*za(p1,p6)*za(p2,p6)*(s345-s(p3,p4)))*prp345/s345)
     &   +Qu*(
     & (za(p1,p5)*zb(p2,p4)*zazb(p3,p1,p5,p6)
     &  /(zb(p1,p5)*za(p1,p6)*s156)
     & +zazb(p1,p2,p6,p4)*(zazb(p5,p2,p6,p4)*za(p4,p3)
     &                    +za(p5,p3)*s(p2,p6))
     &  /(zb(p1,p5)*za(p1,p6)*za(p6,p2)*(s345-s(p3,p4))))*prp34/s(p3,p4)
     & -zazb(p1,p2,p6,p4)**2*za(p4,p5)
     &   /(zb(p3,p5)*za(p1,p6)*za(p2,p6)*(s345-s(p3,p4)))*prp345/s345)
           
      aLL=zb(p2,p4)**2/(zb(p1,p6)*zb(p2,p6)*(s345-s(p3,p4)))*(
     &  Qu*(+zazb(p5,p3,p4,p1)/zb(p3,p4)*prp34
     &      -zb(p1,p5)*za(p4,p5)/zb(p3,p5)*prp345)/zb(p1,p5)
     & +Qd*(+zazb(p5,p3,p4,p2)/zb(p4,p3)*prp34
     &      +zb(p2,p5)*za(p4,p5)/zb(p3,p5)*prp345)/zb(p2,p5))

      aRL=Qu*(
     & (-zb(p2,p4)*za(p6,p1)*zazb(p3,p2,p4,p5)
     &  /(za(p1,p5)*zb(p6,p1)*s156)
     &  -zazb(p3,p1,p6,p2)*(zb(p4,p3)*zb(p2,p5)*za(p3,p1)
     &                     +zb(p4,p5)*zb(p2,p6)*za(p6,p1))
     &  /(za(p1,p5)*zb(p2,p6)*zb(p6,p1)*(s345-s(p3,p4))))*prp34/s(p3,p4)
     & -zazb(p3,p1,p6,p2)**2*zb(p4,p5)
     &   /(za(p3,p5)*zb(p2,p6)*zb(p1,p6)*(s345-s(p3,p4)))*prp345/s345)
     &   +Qd*(
     & (zb(p2,p5)*za(p1,p3)*zazb(p6,p2,p5,p4)
     &  /(za(p2,p5)*zb(p2,p6)*s256)
     & +zazb(p3,p1,p6,p2)*(zazb(p3,p1,p6,p5)*zb(p3,p4)
     &                    +zb(p5,p4)*s(p1,p6))
     &  /(za(p2,p5)*zb(p2,p6)*zb(p6,p1)*(s345-s(p3,p4))))*prp34/s(p3,p4)
     & +zazb(p3,p1,p6,p2)**2*zb(p4,p5)
     &   /(za(p3,p5)*zb(p2,p6)*zb(p1,p6)*(s345-s(p3,p4)))*prp345/s345)
      
c--- additional anomalous coupling component, Eqs. (11)-(16) of hep-ph/0002138
c---  (multiplied by -i)
c--- note: there is a typo in Eq. (15), <45> should be replaced by [45]
c--- note partial fractioning of W propagators to go beyond zerowidth approx.
      aRR=aRR+(prp34-prp345)*(Qd-Qu)
     & *zb(p4,p5)*(za(p1,p2)*zb(p2,p5)+za(p1,p6)*zb(p6,p5))
     & /(2._dp*s(p3,p4)**2*(s345-s(p3,p4))*za(p1,p6)*za(p6,p2))
     & *(xlambda_g*za(p3,p4)*za(p1,p5)*zb(p4,p5)
     &  +(xdelk_g+xlambda_g)*za(p1,p3)*s(p3,p4))

      aLR=aLR+(prp34-prp345)*(Qd-Qu)
     & *za(p3,p5)*za(p1,p5)
     & /(2._dp*s(p3,p4)**2*(s345-s(p3,p4))*za(p1,p6)*za(p6,p2))
     & *(xlambda_g*za(p3,p5)*zb(p3,p4)
     &    *(za(p1,p2)*zb(p2,p5)+za(p1,p6)*zb(p6,p5))
     &  -(xdelk_g+xlambda_g)
     &    *(za(p1,p2)*zb(p2,p4)+za(p1,p6)*zb(p6,p4))*s(p3,p4))

      aLL=aLL+(prp34-prp345)*(Qu-Qd)
     & *za(p3,p5)*(zb(p2,p1)*za(p1,p5)+zb(p2,p6)*za(p6,p5))
     & /(2._dp*s(p3,p4)**2*(s345-s(p3,p4))*zb(p2,p6)*zb(p6,p1))
     & *(xlambda_g*zb(p4,p3)*zb(p2,p5)*za(p3,p5)
     &  +(xdelk_g+xlambda_g)*zb(p2,p4)*s(p3,p4))

      aRL=aRL+(prp34-prp345)*(Qu-Qd)
     & *zb(p4,p5)*zb(p2,p5)
     & /(2._dp*s(p3,p4)**2*(s345-s(p3,p4))*zb(p2,p6)*zb(p6,p1))
     & *(xlambda_g*zb(p4,p5)*za(p4,p3)
     &    *(zb(p2,p1)*za(p1,p5)+zb(p2,p6)*za(p6,p5))
     &  -(xdelk_g+xlambda_g)
     &    *(zb(p2,p1)*za(p1,p3)+zb(p2,p6)*za(p6,p3))*s(p3,p4))

      ubdgmsq=abs(aLL)**2+abs(aRR)**2+abs(aRL)**2+abs(aLR)**2
      
      return
      end

