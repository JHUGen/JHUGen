      subroutine stringsh(q,h5,docc)
! input: momenta q, helicity of 5 is h5,
!        flag docc to turn on c.c. - appropriate for t~ calculation
! computes currents for given q and helicity combination
! stores as common block as will be re-used.
! add more currents as needed
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'currentdecl.f'
      double precision q(mxpart,4)
      integer h1,h2,h5,h6
      integer m,fi,nu,ro,si,i,j
      parameter(m=4)
      double complex gam0(4,4),gam1(m,4,4),gam2(m,m,4,4),
     & gam3(m,m,m,4,4),gam4(m,m,m,m,4,4),gam5(m,m,m,m,m,4,4),
     & f1(4),f2(4),f5(4),f6(4)
      double complex p1(4),p2(4),p5(4),p6(4)
c     p3(4),p4(4)
      common/gamodd/gam1,gam3,gam5
      common/gameven/gam0,gam2,gam4
      parameter(h1=-1,h2=-1,h6=-1)
      logical doJ61,doJ52,docc

      p1(:)=dcmplx(q(1,:))    
      p2(:)=dcmplx(q(2,:))
c      p3(:)=dcmplx(q(3,:))
c      p4(:)=dcmplx(q(4,:))
      p5(:)=dcmplx(q(5,:))
      p6(:)=dcmplx(q(6,:))

c--- decide whick strings to comput; note, this is tied
c--- to order of array in qq_tchan_htq_v.f
c--- it relies on h5 = -1,1

      doJ52=.true. 
c--- first call, calculate everything 
      if     (h5 .eq. -1) then
        doJ61=.true.
      elseif (h5 .eq. +1) then
        doJ61=.false.
      else
        write(6,*) 'stringsh: unknown h5=',h5
      stop
      endif        

      if (doJ61) then
      call ubarspinor0(p6,h6,f6)
      call uspinor0(p1,h1,f1)
      endif
      if (doJ52) then
      call UbKlSt(p5,mt,p2,h5,f5)      ! Auxiliary vector p2
      call uspinor0(p2,h2,f2)
      endif
      
      if (doJ52) then
      J52x0=czip
      do i=1,4
      do j=1,4
      if (abs(gam0(i,j)) .gt. 1d-8) then  
      J52x0=J52x0+f5(i)*gam0(i,j)*f2(j)      
      endif
      enddo
      enddo
      endif
      
      do fi=1,4
      if (doJ52) J52x1(fi)=czip
      if (doJ61) J61x1(fi)=czip
      do i=1,4
      do j=1,4
      if (abs(gam1(fi,i,j)) .gt. 1d-8) then  
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
      
c      write(6,*) 'J52',j52


      do fi=1,4
      do nu=1,4
      do ro=1,4
      if (doJ52) J52x3(fi,nu,ro)=czip
      if (doJ61) J61x3(fi,nu,ro)=czip
      do i=1,4
      do j=1,4
      if (abs(gam3(fi,nu,ro,i,j)) .gt. 1d-8) then  
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
              enddo
            enddo
          enddo
        endif
      endif      

      end


c      double precision function vec(p1,ro)
c      implicit none
c      double precision p1(4)
c      integer ro
c      vec=p1(ro)
c      return
c      end
