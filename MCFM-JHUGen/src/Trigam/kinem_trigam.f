c--- kinematic file for Wga production
c--- note: this is adapted from the point used in Appendix B of 0710.1832

c--- \begin{array}{rrrrrrr}
c---      p_1 & = true( 0.00000000000000, 0.00000000000000, 1021.22119318758,-1021.22119318758
c---      p_2 & = true( 0.00000000000000, 0.00000000000000,-238.714576090637,-238.714576090637
c---      p_3 & = true(-71.5344542606618,-183.877222508616,-3.11006048502754, 197.326337775966
c---      p_4 & = true(-9.92033815503652,-76.1125676676337, 49.0057636944973, 91.0664644166627
c---      p_5 & = true( 32.5059044554765, 245.099246845329,-495.644737899924, 553.889863453468
c---      p_6 & = true( 64.1786550096635, 124.613643938661,-207.896850811885, 250.736037681104
c---      p_7 & = true(-15.2297670494417,-109.723100607740,-124.860731594604, 166.917065951017).\\

c--- momenta are written in Kirill's notation true(E,px,py,pz)
      subroutine set_kinem_trigam(p) 
      include 'constants.f' 
      double precision p(mxpart,4) 
      double precision theta,phi,rho,csig,muk
      double precision ssig,p1true(4),p2true(4),p3true(4)
      double precision p4true(4),p5true(4)
      double precision p1(4),p2(4),p3(4),p4(4),p5(5)
      integer nu,i

      p(:,:)=0d0 
      theta=pi/4d0
      phi=pi/6d0
      rho=pi/3d0
      csig=-7d0/19d0
      muk=6d0

      ssig=sqrt(1d0-csig**2)
      p1true(1)=0.5d0*muk
      p1true(2)=-p1true(1)*sin(theta)
      p1true(3)=-p1true(1)*cos(theta)*sin(phi)
      p1true(4)=-p1true(1)*cos(theta)*cos(phi)

      p2true(1)=muk/3d0
      p2true(2)=p2true(1)
      p2true(3)=0d0
      p2true(4)=0d0

      p3true(1)=-muk/7d0
      p3true(2)=p3true(1)*csig
      p3true(3)=p3true(1)*ssig
      p3true(4)=0d0

      p4true(1)=-0.5d0*muk
      p4true(2)=+p4true(1)*sin(theta)
      p4true(3)=+p4true(1)*cos(theta)*sin(phi)
      p4true(4)=+p4true(1)*cos(theta)*cos(phi)

c--- now construct massless momentum satisfying momentum conservation
c--- (re-use rho,muk here)
      do nu=1,4
      p5true(nu)=-p1true(nu)-p2true(nu)-p3true(nu)-p4true(nu)
      enddo
      muk=p5true(1)**2-p5true(2)**2-p5true(3)**2-p5true(4)**2
      rho=2d0*(p4true(1)*p5true(1)-p4true(2)*p5true(2)
     &        -p4true(3)*p5true(3)-p4true(4)*p5true(4))
      do nu=1,4
      p5true(nu)=-muk/rho*p4true(nu)+p5true(nu)
      enddo
      do nu=1,4
      p4true(nu)=(1d0+muk/rho)*p4true(nu)
      enddo
                  
c--- now form the momenta that will be used in the Kirill routines
      do nu=1,4
      p1(nu)=p1true(nu)
      p2(nu)=p2true(nu)
      p3(nu)=p3true(nu)
      p4(nu)=p4true(nu)
      p5(nu)=p5true(nu)
      enddo

!===== MCFM p 
      p(1,1)=p1(2)
      p(1,2)=p1(3)
      p(1,3)=p1(4)
      p(1,4)=p1(1)
!===== MCFM p 
      p(2,1)=p2(2)
      p(2,2)=p2(3)
      p(2,3)=p2(4)
      p(2,4)=p2(1)
!===== MCFM p 
      p(3,1)=p3(2)
      p(3,2)=p3(3)
      p(3,3)=p3(4)
      p(3,4)=p3(1)
!===== MCFM p 
      p(4,1)=p4(2)
      p(4,2)=p4(3)
      p(4,3)=p4(4)
      p(4,4)=p4(1)
!===== MCFM p 
      p(5,1)=p5(2)
      p(5,2)=p5(3)
      p(5,3)=p5(4)
      p(5,4)=p5(1)

      return 
      end
