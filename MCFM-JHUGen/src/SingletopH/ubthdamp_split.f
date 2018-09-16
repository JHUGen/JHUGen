      subroutine ubthdamp(p,p1,p2,p3,p4,p6,amp)
      implicit none
      include 'types.f'
C     Matrix element
c     u(-p1)+b(p2)->H(p3,p4)+t(p5)+d(p6)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'alpha1.f'
      include 'anomHiggs.f'
      real(dp):: p(mxpart,4),q(mxpart,4),dot,p5Deta,omal
      real(dp):: s16,s25,s345,s34,mw
      integer:: nu,p1,p2,p3,p4,k5,e5,p6,eta,j,k
      complex(dp):: prW,prt,iza,izb,amp(2)
      complex(dp):: zab2,amplower(2),ampmiddle(2),ampscalar(2)
      parameter (k5=5,e5=7)

      prW(s34)=cone/cplx2(s34-wmass**2,zip)
      prt(s34)=cone/cplx2(s34-mt**2,zip)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)

      omal=1d0-alpha1

      mw=wmass
      eta=p2
C   Construct demassified momentum for p5 wrt eta
C   and store in position e5
      p5Deta=dot(p,5,eta)
      do nu=1,4
      do j=1,6
      q(j,nu)=p(j,nu)
      enddo
      q(k5,nu)=p(5,nu)-mt**2*p(eta,nu)/(2d0*p5Deta)
      q(e5,nu)=mt**2*p(eta,nu)/(2d0*p5Deta)
      enddo
      if (p3 == 3) then
        call spinoru(7,q,za,zb)
      else
        call spinoru(7,q,zb,za)
      endif
      s16=s(p1,p6)
      s25=s(p1,p3)+s(p1,p4)+s(p1,p6)+s(p3,p4)+s(p3,p6)+s(p4,p6)
      s345=s(p1,p2)+s(p1,p6)+s(p2,p6)
c      call writeout(q)
c      pause
      Amplower(1)= + prW(s16)*prt(s345)*mt**2*mw**(-1) * ( 1.D0/2.D0*
     &    izb(k5,e5)*zb(p1,p2)*zab2(p6,p1,p2,e5) - 1.D0/2.D0*za(p6,k5)*
     &    zb(p1,p2) )

      Amplower(2)= + prW(s16)*prt(s345)*mt*mw**(-1) * (  - 1.D0/2.D0*
     &    zb(p1,p2)*zab2(p6,p1,p2,k5) )
      Amplower(2) = Amplower(2) + prW(s16)*prt(s345)*mt**3*mw**(-1)
     &  * ( 1.D0/2.D0*iza(k5,e5)*za(p6,e5)*zb(p1,p2) )


      Ampmiddle(1)= + prW(s16)*prW(s25)*mw * (  - za(p6,k5)*zb(p1,p2) )
      Ampmiddle(1) = Ampmiddle(1) + prW(s16)*prW(s25)*mt**2*mw**(-1)
     &  * ( 1.D0/2.D0*izb(k5,e5)*zb(p2,e5)*zab2(p6,p3,p4,p1)*omal )

      Ampmiddle(2)= + prW(s16)*prW(s25)*mt*mw**(-1) * (  - 1.D0/2.D0*
     &    zb(p2,k5)*zab2(p6,p3,p4,p1)*omal )
      Ampmiddle(2) = Ampmiddle(2) + prW(s16)*prW(s25)*mt*mw * ( iza(k5,
     &    e5)*za(p6,e5)*zb(p1,p2) )


      Ampscalar(1)= + prW(s16)*prW(s25)*mt**2*mw**(-1) * ( 1.D0/2.D0*
     &    izb(k5,e5)*zb(p2,e5)*zab2(p6,p3,p4,p1)*alpha1 )

      Ampscalar(2)= + prW(s16)*prW(s25)*mt*mw**(-1) * (  - 1.D0/2.D0*
     &    zb(p2,k5)*zab2(p6,p3,p4,p1)*alpha1 )


c  -- include anomalous Yukawa and EW couplings for the Higgs
      do k=1,2
      Amplower(k)=cttH*Amplower(k)
      Ampmiddle(k)=cWWH*Ampmiddle(k)
      Ampscalar(k)=cWWH*Ampscalar(k)
      Amp(k)=Ampmiddle(k)+Amplower(k)+Ampscalar(k)
      enddo

      return
      end
      
      
      
      
      
      subroutine ubthdamp_dk(p,p1,p2,p3,p4,p6,pmless,amp)
      implicit none
      include 'types.f'
C     Matrix element 
c     u(-p1)+b(p2)->H(p3,p4)+t(p5)+d(p6)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'alpha1.f'
      include 'anomHiggs.f'
      real(dp):: p(mxpart,4),q(mxpart,4),p5Deta,omal
      real(dp):: pmless(4)
      real(dp):: s16,s25,s345,s34,mw
      integer:: nu,p1,p2,p3,p4,k5,e5,p6,j
      complex(dp):: prW,prt,izb,Amp(2)
      complex(dp):: zab2,amplower(2),ampmiddle(2),ampscalar(2)
      parameter (k5=5,e5=7)

      prW(s34)=cone/cplx2(s34-wmass**2,zip)
      prt(s34)=cone/cplx2(s34-mt**2,zip)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      izb(p1,p2)=cone/zb(p1,p2)

      omal=1d0-alpha1

      mw=wmass
C   Construct demassified momentum for p5 wrt pmless
C   and store in position e5
       p5Deta=p(5,4)*pmless(4)-p(5,1)*pmless(1)-p(5,2)*pmless(2)
     & - p(5,3)*pmless(3)
      do nu=1,4
      do j=1,6
      q(j,nu)=p(j,nu)
      enddo
      q(k5,nu)=p(5,nu)-mt**2*pmless(nu)/(2d0*p5Deta)
      q(e5,nu)=mt**2*pmless(nu)/(2d0*p5Deta)
      enddo
      if (p3 == 3) then
        call spinoru(7,q,za,zb)
      else
        call spinoru(7,q,zb,za)
      endif
      s16=s(p1,p6)
      s25=s(p1,p3)+s(p1,p4)+s(p1,p6)+s(p3,p4)+s(p3,p6)+s(p4,p6)
      s345=s(p1,p2)+s(p1,p6)+s(p2,p6)
      Amplower(1)= + prW(s16)*prt(s345)*mt**2*mw**(-1) * ( 1.D0/2.D0*
     &    izb(k5,e5)*zb(p1,p2)*zab2(p6,p1,p2,e5) - 1.D0/2.D0*za(p6,k5)*
     &    zb(p1,p2) )

      Ampmiddle(1)= + prW(s16)*prW(s25)*mw * (  - za(p6,k5)*zb(p1,p2) )
      Ampmiddle(1) = Ampmiddle(1) + prW(s16)*prW(s25)*mt**2*mw**(-1)
     &  * ( 1.D0/2.D0*izb(k5,e5)*zb(p2,e5)*zab2(p6,p3,p4,p1)*omal )

      Ampscalar(1)= + prW(s16)*prW(s25)*mt**2*mw**(-1) * ( 1.D0/2.D0*
     &    izb(k5,e5)*zb(p2,e5)*zab2(p6,p3,p4,p1)*alpha1 )


c  -- include anomalous Yukawa and EW couplings for the Higgs
      Amplower(1)=cttH*Amplower(1)
      Ampmiddle(1)=cWWH*Ampmiddle(1)
      Ampscalar(1)=cWWH*Ampscalar(1)
      Amp(1)=Ampmiddle(1)+Amplower(1)+Ampscalar(1)
      return
      end
