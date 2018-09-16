!----- kinematic file for WpWmg production
c--- note: this is the same point as used in Appendix B of 0710.1832

c--- \begin{array}{rrrrrrr}
c---      p_1 & = true( 0.00000000000000, 0.00000000000000, 1021.22119318758,-1021.22119318758
c---      p_2 & = true( 0.00000000000000, 0.00000000000000,-238.714576090637,-238.714576090637
c---      p_3 & = true(-71.5344542606618,-183.877222508616,-3.11006048502754, 197.326337775966
c---      p_4 & = true(-9.92033815503652,-76.1125676676337, 49.0057636944973, 91.0664644166627
c---      p_5 & = true( 32.5059044554765, 245.099246845329,-495.644737899924, 553.889863453468
c---      p_6 & = true( 64.1786550096635, 124.613643938661,-207.896850811885, 250.736037681104
c---      p_7 & = true(-15.2297670494417,-109.723100607740,-124.860731594604, 166.917065951017).\\

c--- momenta are written in Kirill's notation true(E,px,py,pz)

      theta=pi/4._dp
      phi=pi/6._dp
      rho=pi/3._dp
      csig=-7._dp/19._dp

      muk=6._dp
      ssig=sqrt(1._dp-csig**2)

      p1true(1)=-0.5_dp*muk
      p1true(2)=-p1true(1)*sin(theta)
      p1true(3)=-p1true(1)*cos(theta)*sin(phi)
      p1true(4)=-p1true(1)*cos(theta)*cos(phi)

      p4true(1)=muk/3._dp
      p4true(2)=p4true(1)
      p4true(3)=0._dp
      p4true(4)=0._dp

      p3true(1)=muk/7._dp
      p3true(2)=p3true(1)*csig
      p3true(3)=p3true(1)*ssig
      p3true(4)=0._dp

      p2true(1)=-0.5_dp*muk
      p2true(2)=+p2true(1)*sin(theta)
      p2true(3)=+p2true(1)*cos(theta)*sin(phi)
      p2true(4)=+p2true(1)*cos(theta)*cos(phi)

      p5true(1)=muk/6._dp
      p5true(2)=+p5true(1)*cos(rho)*csig
      p5true(3)=+p5true(1)*cos(rho)*ssig
      p5true(4)=+p5true(1)*sin(rho)


      do nu=1,4
      p6true(nu)=-p1true(nu)-p2true(nu)-p3true(nu)-p4true(nu)-p5true(nu)
      enddo

c--- BEGIN: ROTATION TO KEEP BOTH CODES HAPPY
c      tmp=p1true(2)
c      p1true(2)=p1true(4)
c      p1true(3)=-p1true(3)
c      p1true(4)=tmp

c      tmp=p2true(2)
c      p2true(2)=p2true(4)
c      p2true(3)=-p2true(3)
c      p2true(4)=tmp

c      tmp=p3true(2)
c      p3true(2)=p3true(4)
c      p3true(3)=-p3true(3)
c      p3true(4)=tmp

c      tmp=p4true(2)
c      p4true(2)=p4true(4)
c      p4true(3)=-p4true(3)
c      p4true(4)=tmp

c      tmp=p5true(2)
c      p5true(2)=p5true(4)
c      p5true(3)=-p5true(3)
c      p5true(4)=tmp

c      tmp=p6true(2)
c      p6true(2)=p6true(4)
c      p6true(3)=-p6true(3)
c      p6true(4)=tmp

c--- END: ROTATION TO KEEP BOTH CODES HAPPY

c--- another rotation
c      tmp=p1true(3)
c      p1true(3)=+tmp*cos(0.177_dp)+p1true(4)*sin(0.177_dp)     
c      p1true(4)=-tmp*sin(0.177_dp)+p1true(4)*cos(0.177_dp)     

c      tmp=p2true(3)
c      p2true(3)=+tmp*cos(0.177_dp)+p2true(4)*sin(0.177_dp)     
c      p2true(4)=-tmp*sin(0.177_dp)+p2true(4)*cos(0.177_dp)     

c      tmp=p3true(3)
c      p3true(3)=+tmp*cos(0.177_dp)+p3true(4)*sin(0.177_dp)     
c      p3true(4)=-tmp*sin(0.177_dp)+p3true(4)*cos(0.177_dp)     

c      tmp=p4true(3)
c      p4true(3)=+tmp*cos(0.177_dp)+p4true(4)*sin(0.177_dp)     
c      p4true(4)=-tmp*sin(0.177_dp)+p4true(4)*cos(0.177_dp)     

c      tmp=p5true(3)
c      p5true(3)=+tmp*cos(0.177_dp)+p5true(4)*sin(0.177_dp)     
c      p5true(4)=-tmp*sin(0.177_dp)+p5true(4)*cos(0.177_dp)     

c      tmp=p6true(3)
c      p6true(3)=+tmp*cos(0.177_dp)+p6true(4)*sin(0.177_dp)     
c      p6true(4)=-tmp*sin(0.177_dp)+p6true(4)*cos(0.177_dp)     

c      tmp=p7true(3)
c      p7true(3)=+tmp*cos(0.177_dp)+p7true(4)*sin(0.177_dp)     
c      p7true(4)=-tmp*sin(0.177_dp)+p7true(4)*cos(0.177_dp)     

c--- write out momenta
c      do nu=1,4
c      write(6,39) 'p(1,',nu,')=',p1true(nu)
c      write(6,39) 'p(2,',nu,')=',p2true(nu)
c      write(6,39) 'p(3,',nu,')=',p3true(nu)
c      write(6,39) 'p(4,',nu,')=',p4true(nu)
c      write(6,39) 'p(5,',nu,')=',p5true(nu)
c      write(6,39) 'p(6,',nu,')=',p6true(nu)
c      enddo
      
c   39 format(a4,i2,a2,g25.16)   
      
c--- now form the momenta that will be used in the Kirill routines
c      do nu=1,4
c      p1(nu)=p1true(nu)
c      p2(nu)=p2true(nu)
c      p3(nu)=p3true(nu)+p4true(nu)
c      p4(nu)=p5true(nu)+p6true(nu)
c      enddo
      
