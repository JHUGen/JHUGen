      subroutine strings_dk(p,q,h3,h5,docc)
! input: momenta of decay products is q, momenta prior to decay is p,
! helicity of 5 is h5
! computes currents for given q and helicity combination
! stores as common block as will be re-used.
! add more currents as needed
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'currentdecl.f'
      double precision q(mxpart,4)
      integer h1,h2,h3,h4,h5,h6
      integer m,fi,nu,ro,si,om,i,j
      parameter(m=4)
      double complex gam0(4,4),gam1(m,4,4),gam2(m,m,4,4),
     & gam3(m,m,m,4,4),gam4(m,m,m,m,4,4),gam5(m,m,m,m,m,4,4),
     & f1(4),f2(4),f3(4),f4(4),f5(4),f6(4)
      double complex p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),
     & q1(4),q2(4),q3(4),q4(4),q5(4),q6(4),q7(4),q8(4)
      double precision p(mxpart,4)
      common/gamodd/gam1,gam3,gam5
      common/gameven/gam0,gam2,gam4
      parameter(h1=-1,h2=-1,h6=-1)
      logical doJ61,doJ52,docc


c -- p(q1) p(q2) e-(q3) e+(q4) nu(q5) mu+(q6) b(q7) d(q8)
      q1(:)=dcmplx(q(1,:))    
      q2(:)=dcmplx(q(2,:))
      q3(:)=dcmplx(q(3,:))
      q4(:)=dcmplx(q(4,:))
      q5(:)=dcmplx(q(5,:))
      q6(:)=dcmplx(q(6,:))
      q7(:)=dcmplx(q(7,:))
      q8(:)=dcmplx(q(8,:))

c -- p(p1) p(p2) e-(p3) e+(p4) t(p5) d(p6)
      p1(:)=dcmplx(p(1,:))    
      p2(:)=dcmplx(p(2,:))
      p3(:)=dcmplx(p(3,:))
      p4(:)=dcmplx(p(4,:))
      p5(:)=dcmplx(p(5,:))
      p6(:)=dcmplx(p(6,:))

c--- anti-lepton helicity h4 set by h3
      h4=h3


c -- decide which strings to compute - for the moment, calculate everything
c -- can modify later to speed up
c      doJ61=.false.
c      doJ52=.false. 
c--- first call, calculate everything 
c      if ((h3 .eq. -1) .and. (h5 .eq. -1)) then
c        doJ61=.true.
c	doJ52=.true.
c      endif        
c      if ((h3 .eq. -1) .and. (h5 .eq. +1)) then
c	doJ52=.true.
c      endif        
      doJ61=.true.
      doJ52=.true.


      if (doJ61) then
      call ubarspinor0(p6,h6,f6)
      call uspinor0(p1,h1,f1)
      endif
      if (doJ52) then
      call UbKlSt(p5,mt,q6,h5,f5)	! Auxiliary vector q6
      call uspinor0(p2,h2,f2)
      endif
      call ubarspinor0(p3,h3,f3)
      call uspinor0(p4,h4,f4)

c      J34x0=czip
c      J61x0=czip      
      if (doJ52) then
      J52x0=czip
      do i=1,4
      do j=1,4
c      J34x0=J34x0+f3(i)*gam0(i,j)*f4(j)
      if (abs(gam0(i,j)) .gt. 1d-8) then  
      J52x0=J52x0+f5(i)*gam0(i,j)*f2(j)      
c      J61x0=J61x0+f6(i)*gam0(i,j)*f1(j)
      endif
      enddo
      enddo
      endif
      
      do fi=1,4
      J34x1(fi)=czip
      if (doJ52) J52x1(fi)=czip
      if (doJ61) J61x1(fi)=czip
      do i=1,4
      do j=1,4
      if (abs(gam1(fi,i,j)) .gt. 1d-8) then  
      J34x1(fi)=J34x1(fi)+f3(i)*gam1(fi,i,j)*f4(j)  
      if (doJ52) J52x1(fi)=J52x1(fi)+f5(i)*gam1(fi,i,j)*f2(j)       
      if (doJ61) J61x1(fi)=J61x1(fi)+f6(i)*gam1(fi,i,j)*f1(j)    
      endif   
      enddo
      enddo
      enddo

      if (doJ52) then
      do fi=1,4
      do nu=1,4
      J52x2(fi,nu)=czip
      do i=1,4
      do j=1,4
      if (abs(gam2(fi,nu,i,j)) .gt. 1d-8) then  
      J52x2(fi,nu)=J52x2(fi,nu)
     & +f5(i)*gam2(fi,nu,i,j)*f2(j)      
      endif
      enddo
      enddo
      enddo
      enddo
      endif
      
c      write(6,*) 'J34',j34
c      write(6,*) 'J52',j52


      do fi=1,4
      do nu=1,4
      do ro=1,4
      J34x3(fi,nu,ro)=czip
      if (doJ52) J52x3(fi,nu,ro)=czip
      if (doJ61) J61x3(fi,nu,ro)=czip
      do i=1,4
      do j=1,4
      if (abs(gam3(fi,nu,ro,i,j)) .gt. 1d-8) then  
      J34x3(fi,nu,ro)=J34x3(fi,nu,ro)
     & +f3(i)*gam3(fi,nu,ro,i,j)*f4(j)      
      if (doJ52) then
      J52x3(fi,nu,ro)=J52x3(fi,nu,ro)
     & +f5(i)*gam3(fi,nu,ro,i,j)*f2(j) 
      endif
      if (doJ61) then     
      J61x3(fi,nu,ro)=J61x3(fi,nu,ro)
     & +f6(i)*gam3(fi,nu,ro,i,j)*f1(j)
      endif
      endif      
      enddo
      enddo
      enddo
      enddo
      enddo


      if (doJ52) then
      do fi=1,4
      do nu=1,4
      do ro=1,4
      do si=1,4
      J52x4(fi,nu,ro,si)=czip
      do i=1,4
      do j=1,4
      if (abs(gam4(fi,nu,ro,si,i,j)) .gt. 1d-8) then  
      J52x4(fi,nu,ro,si)=J52x4(fi,nu,ro,si)
     & +f5(i)*gam4(fi,nu,ro,si,i,j)*f2(j)
      endif      
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      endif

      do fi=1,4
      do nu=1,4
      do ro=1,4
      do si=1,4
      do om=1,4
      if (doJ52) J52x5(fi,nu,ro,si,om)=czip
      if (doJ61) J61x5(fi,nu,ro,si,om)=czip
      do i=1,4
      do j=1,4
      if (abs(gam5(fi,nu,ro,si,om,i,j)) .gt. 1d-8) then  
      if (doJ52) then
      J52x5(fi,nu,ro,si,om)=J52x5(fi,nu,ro,si,om)
     & +f5(i)*gam5(fi,nu,ro,si,om,i,j)*f2(j)      
      endif
      if (doJ61) then
      J61x5(fi,nu,ro,si,om)=J61x5(fi,nu,ro,si,om)
     & +f6(i)*gam5(fi,nu,ro,si,om,i,j)*f1(j)   
      endif
      endif   
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

c--- perform complex conjugation, if required
      if (docc) then
        if (doJ52) then
	  J52x0=-h5*dconjg(J52x0)
	  do fi=1,4
	  J52x1(fi)=-h5*dconjg(J52x1(fi))
	    do nu=1,4
	    J52x2(fi,nu)=-h5*dconjg(J52x2(fi,nu))
	      do ro=1,4
	      J52x3(fi,nu,ro)=-h5*dconjg(J52x3(fi,nu,ro))
	        do si=1,4
                J52x4(fi,nu,ro,si)=-h5*dconjg(J52x4(fi,nu,ro,si))
                  do om=1,4
                  J52x5(fi,nu,ro,si,om)=-h5*
     &                              dconjg(J52x5(fi,nu,ro,si,om))
                  enddo
               enddo
            enddo
         enddo
      enddo
        endif
        if (doJ61) then
	  do fi=1,4
	  J61x1(fi)=dconjg(J61x1(fi))
	    do nu=1,4
	      do ro=1,4
	      J61x3(fi,nu,ro)=dconjg(J61x3(fi,nu,ro))
                do si=1,4
	          do om=1,4
	          J61x5(fi,nu,ro,si,om)=dconjg(J61x5(fi,nu,ro,si,om))
	          enddo
	        enddo
	      enddo
            enddo
          enddo
        endif

        do fi=1,4
         J34x1(fi)=dconjg(J34x1(fi))
         do nu=1,4
            do ro=1,4
               J34x3(fi,nu,ro)=dconjg(J34x3(fi,nu,ro))
            enddo
         enddo
      enddo   

      endif

            end


