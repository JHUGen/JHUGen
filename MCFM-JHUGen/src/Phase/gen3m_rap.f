      subroutine gen3m_rap(r,p,m3,m4,wt,*)
      implicit none
      include 'types.f'
      

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'x1x2.f'
c---- generate phase space for 2-->3 process
c---- with 3 and 4 (masses m3,m4) and 5 massless.
c---- r(mxdim) and 
c---- p1(4),p2(4) are input momenta reversed in sign from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5
c---- with all 2 pi)
      real(dp):: r(mxdim),y5starmin,y5starmax,plstar,plstarsq,
     & Estar,p(mxpart,4),a,E34st,
     & wt,p3(4),p4(4),p5(4),p345(4),p34(4),pstsq,
     & xmin,pt5,ymin,ymax,phi,wt34,
     & dely,sinhy,coshy,y,rtshat,pt2,
     & vs,vsqmax,vsqmin,s34,sinhy5,coshy5,y5,y5max,s34max,s34min,
     & m3,m4,xjac,w,wmax,wmin
c      real(dp):: p3cm(4),beta,costh,sinth
      integer:: j,nu
      include 'energy.f'
      wt=0._dp

      
C---set all vectors to zero
      do nu=1,4
      p345(nu)=0._dp
      do j=1,mxpart
          p(j,nu)=0._dp
      enddo
      enddo 


      xjac=0.5_dp/((twopi)**3*sqrts**2)


C--- generate PT5
      wmin=log(1.e-5_dp)
      wmax=log((sqrts-m3-m4)/2._dp)
      w=wmin+(wmax-wmin)*r(1)
      pt5=exp(w)
      pt2=pt5*pt5
C express in terms of dptsq
      xjac=xjac*pt5**2*(wmax-wmin)
C--generate rapidity
c--- rapidity limited by sqrts=2*pT*coshy
      a=sqrts/(pt5+sqrt(pt2+(m3+m4)**2))
      y5max=log(a+sqrt(a**2-1._dp))
      y5=y5max*(2._dp*r(2)-1._dp)
      sinhy5=sinh(y5)
      coshy5=sqrt(1._dp+sinhy5**2)
      xjac=xjac*2._dp*y5max
        
C--generate phi
      phi=twopi*r(3)
      xjac=xjac*twopi

C in lab frame
      p5(4)=pt5*coshy5
      p5(1)=pt5*cos(phi)
      p5(2)=pt5*sin(phi)
      p5(3)=pt5*sinhy5
        

C  s34=(p_1+p_2-p_5)^2
      s34max=sqrts**2-2._dp*sqrts*pt5
      s34min=(m3+m4)**2
      vsqmax=1._dp/s34min
      vsqmin=1._dp/s34max
      if (vsqmin > vsqmax) then
        write(6,*) 'gen3m:vsqmin',vsqmin
        write(6,*) 'gen3m:vsqmax',vsqmax
      return 1
      endif
      xmin=vsqmin/vsqmax
      vs=(vsqmax-vsqmin)*r(4)+vsqmin
      s34=1/vs
      xjac=xjac*(vsqmax-vsqmin)*s34**2

c--- invariant mass of jets

C plstar is longitudinal momentum of p5 (and -p34) in 34-5 centre of mass
C plstar is obtained by solving s=(p5+p34)^2      
      plstarsq=((sqrts**2-s34)**2-4._dp*pt2*sqrts**2)/(4._dp*sqrts**2)
      if (plstarsq <= 0._dp) then
      write(6,*) 'gen3m:plstarsq,s34,pt2',plstarsq,s34,pt2
      return 1
      endif

      plstar=sqrt(plstarsq)
      Estar=sqrt(plstarsq+pt2)
      y5starmax=log((Estar+plstar)/pt5)
 
      y5starmin=-y5starmax
      ymax=y5-y5starmin
      ymin=y5-y5starmax
      dely=ymax-ymin
      y=ymin+r(5)*dely     
      sinhy=sinh(y)
      coshy=sqrt(1._dp+sinhy**2)
      xjac=xjac*dely
c--- now calculate in the centre of mass 
      pstsq=pt2+(p5(3)*coshy-p5(4)*sinhy)**2
      E34st=sqrt(s34+pstsq)

      rtshat=E34st+sqrt(pstsq)

c--- back in lab frame
      p345(4)=rtshat*coshy
      p345(3)=rtshat*sinhy
            
      do j=1,4
        p34(j)=p345(j)-p5(j)
      enddo
      
      xx(1)=(p345(4)+p345(3))/sqrts
      xx(2)=(p345(4)-p345(3))/sqrts
      
c      if   (xx(1)*xx(2) > 1._dp) then
c      write(6,*) 'gen3m:xx1*xx2,xx(1),xx(2)',xx(1)*xx(2),xx(1),xx(2)
c      endif

      if   ((xx(1) > 1._dp) .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin).or. (xx(2) < xmin)) return 1
      
      xjac=rtshat/E34st*xjac

c--- now make the initial state momenta
        
      
      p(1,4)=-xx(1)*sqrts/2._dp
      p(1,3)=p(1,4)

      p(2,4)=-xx(2)*sqrts/2._dp
      p(2,3)=-p(2,4)
      
c--- decay s34 lump into two particles with mass m3 and m4
      call phi3m(r(6),r(7),p34,p3,p4,m3,m4,wt34,*99)

c      costh=2._dp*r(6)-1._dp
c      sinth=sqrt(1._dp-costh**2)
c      phi=2._dp*pi*r(7)
c      beta=sqrt(1._dp-(m3+m4)**2/s34)
c      xjac=4._dp*pi*xjac/8._dp

c      p3cm(4)=sqrt(s34+m3**2-m4**2)/2._dp
c      p3cm(1)=p3cm(4)*beta*sinth*cos(phi)
c      p3cm(2)=p3cm(4)*beta*sinth*sin(phi)
c      p3cm(3)=p3cm(4)*beta*costh

      
c--- boost into lab frame    
c      call boost(sqrt(s34),p34,p3cm,p3)

      do j=1,4
      p(3,j)=p3(j)
      p(4,j)=p4(j)
      p(5,j)=p5(j)
      enddo

      wt=xjac*wt34
      return
 99   continue
      wt=0._dp
      return 1
      end
      
      
      
      
      
      
      
      
      
      
      
