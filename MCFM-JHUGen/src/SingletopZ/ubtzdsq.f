      function ubtzdsq(p,j1,j2,j3,j4,j6)
      implicit none
      include 'types.f'
      real(dp):: ubtzdsq
C     Matrix element squared averaged over initial colors and spins
c     u(-p1)+b(p2)->Z(p3,p4)+t(p5)+d(p6)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      real(dp):: p(mxpart,4),ubtzd
      complex(dp):: amp(2,2)
      integer:: j1,j2,j3,j4,j6

      call ubtzdamp(p,j1,j2,j3,j4,j6,amp)
      ubtzd=abs(amp(1,1))**2+abs(amp(1,2))**2
     &     +abs(amp(2,1))**2+abs(amp(2,2))**2

      ubtzdsq=4d0*xn**2*aveqq*gwsq**2*esq**2*ubtzd

      return
      end




      subroutine ubtzdamp(p,j1,j2,j3,j4,j6,amp)
      implicit none
      include 'types.f'
c--- Author: R. Rontsch (15 Feb 2013)
c--- Amplitudes taken from Campbell, Ellis, Rontsch arXiv:1302.3856 (appendix A)
c--- note: overall factor of 1/(s34-mz^2) extracted in a slightly
c--- different way for numerical stability

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4),q(mxpart,4),dot
      complex(dp):: facuLl,facuRl,facdLl,facdRl
      real(dp):: s25,s16,s34,s234,s235,s345,s346,s134,mtsq
      complex(dp):: amplower(2,2),ampmiddle(2,2),ampextra(2,2),
     &     ampupper(2,2),amp(2,2)
      complex(dp):: prW,zab2,cprop,iprZ
      integer:: eta,p1,p2,p3,p4,pl,pa,j,j1,j2,j3,j4,j6,jb,nu
      prW(s16)=cone/cplx2(s16-wmass**2,zip)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      mtsq=mt**2

      jb=j2   ! momentum to make top massless
      eta=7
      do nu=1,4
      q(1,nu)=p(j1,nu)
      q(2,nu)=p(j2,nu)
      q(3,nu)=p(j3,nu)
      q(4,nu)=p(j4,nu)
      q(5,nu)=p(5,nu)-mtsq*p(jb,nu)/(2d0*dot(p,5,jb))
      q(6,nu)=p(j6,nu)
      q(eta,nu)=mtsq*p(jb,nu)/(2d0*dot(p,5,jb))
      enddo
      s16=2d0*dot(p,j1,j6)
      s25=2d0*dot(p,j2,5)+mt**2
      s345=2d0*dot(p,j1,j2)+2d0*dot(p,j2,j6)+2d0*dot(p,j1,j6)
      s134=2d0*dot(p,j1,j3)+2d0*dot(p,j1,j4)+2d0*dot(p,j3,j4)
      s234=2d0*dot(p,j2,j3)+2d0*dot(p,j2,j4)+2d0*dot(p,j3,j4)
      s346=2d0*dot(p,j6,j3)+2d0*dot(p,j6,j4)+2d0*dot(p,j3,j4)
      s34=2d0*dot(p,j3,j4)
      s235=2d0*dot(p,j1,j6)+2d0*dot(p,j1,j4)+2d0*dot(p,j6,j4)

c--- top production (j3=3,j4=4) as-is
      if (j3==3) then
        call spinoru(7,q,za,zb)
      else
c--- antitop production (j3=4,j4=3) requires c.c.
        call spinoru(7,q,zb,za)
      endif

c--- Implementation of Baur-Zeppenfeld treatment of Z width
      cprop=cplx1(1d0/sqrt((s34-zmass**2)**2+(zmass*zwidth)**2))
      iprZ=cplx1(s34-zmass**2)

      do j=1,2
      if (j == 1) then
         pl=3
         pa=4
         facuLl=cplx1(Qu*q1)*iprZ+cplx1(L(2)*le)*s34
         facuRl=cplx1(Qu*q1)*iprZ+cplx1(R(2)*le)*s34
         facdLl=cplx1(Qd*q1)*iprZ+cplx1(L(1)*le)*s34
         facdRl=cplx1(Qd*q1)*iprZ+cplx1(R(1)*le)*s34
      else
         pl=4
         pa=3
         facuLl=cplx1(Qu*q1)*iprZ+cplx1(L(2)*re)*s34
         facuRl=cplx1(Qu*q1)*iprZ+cplx1(R(2)*re)*s34
         facdLl=cplx1(Qd*q1)*iprZ+cplx1(L(1)*re)*s34
         facdRl=cplx1(Qd*q1)*iprZ+cplx1(R(1)*re)*s34
      endif


      ampupper(j,1)=prW(s25)/s34*( facuLl/s134*za(5,6)*zb(1,pa)*
     &     zab2(pl,1,pa,2) - facdLl/s346*za(pl,6)*zb(1,2)*
     &     (zab2(5,pl,6,pa)))

      ampupper(j,2)=prW(s25)*mt/s34*( facuLl/s134*za(6,eta)*zb(1,pa)/
     &     za(5,eta)*zab2(pl,1,pa,2) + facdLl/s346*za(pl,6)*zb(1,2)
     &     /za(5,eta)*zab2(eta,pl,6,pa))

      amplower(j,1)=prW(s16)/s34*( -facuRl*mtsq/(s345-mtsq)*za(pl,6)
     &     *zb(1,2)*zb(pa,eta)/zb(5,eta) + facdLl/s234*zab2(pl,2,pa,1)*
     &     za(6,5)*zb(2,pa) - facuLl/(s345-mtsq)*zab2(6,1,2,pa)
     &     *za(pl,5)*zb(1,2) )

      amplower(j,2)=prW(s16)*mt/s34*(
     &     facuRl/(s345-mtsq)*za(pl,6)*zb(1,2)* zb(pa,5)
     &     - facdLl/s234*zab2(pl,2,pa,1)*za(6,eta)*zb(2,pa)/
     &     za(5,eta)
     &     +facuLl/(s345-mtsq)*zab2(6,1,2,pa)*za(pl,eta)*
     &     zb(1,2)/za(5,eta) )

      if (j == 1) then
         ampextra(j,1)=prW(s25)*prW(s16)/(2d0*xw*s235)*za(pl,5)*zb(1,pa)
     &        *zab2(6,1,pa,2)*iprZ
         ampextra(j,2)=-prW(s25)*prW(s16)*mt/(2d0*xw*s235)*za(pl,eta)*
     &        zb(1,pa)*zab2(6,1,pa,2)/za(5,eta)*iprZ
      else
         ampextra(j,1)=czip
         ampextra(j,2)=czip
      endif
      ampmiddle(j,1)=prW(s25)*prW(s16)/s34*(-(facuLl-facdLl)*
     &     ( zab2(pl,1,6,pa)*za(6,5)*zb(1,2) + zab2(5,1,6,2)*
     &     za(pl,6)*zb(1,pa) + zab2(6,pl,pa,1)*za(pl,5)*zb(2,pa) )
     &     + mtsq/2d0*za(pl,6)*zb(1,pa)*zb(2,eta)/zb(5,eta)*
     &     (facuLl-facdLl-facuRl+facdRl) )

      ampmiddle(j,2)=prW(s25)*prW(s16)*mt/s34*( (facuLl-facdLl)*(
     &     (za(6,eta)*zb(1,2)*zab2(pl,1,6,pa) + zab2(eta,1,6,2)
     &     *za(pl,6)*zb(1,pa) + zab2(6,pl,pa,1)*za(pl,eta)*zb(2,pa))/
     &     za(5,eta)) + za(pl,6)/2d0*zb(1,pa)*zb(2,5)*
     &     (facuLl-facdLl+facuRl-facdRl))


      enddo
      amp=(amplower+ampupper+ampmiddle+ampextra)*cprop

      return
      end

