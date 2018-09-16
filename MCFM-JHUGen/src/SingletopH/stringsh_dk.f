      subroutine stringsh_dk(q,pmless,h5,docc)
      implicit none
      include 'types.f'
! input: momenta q, helicity of 5 is h5,
!        flag docc to turn on c.c. - appropriate for t~ calculation
! input momentum pmless is momentum used to make pt massless
! computes currents for given q and helicity combination
! stores as common block as will be re-used.
! add more currents as needed
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'currentdecl.f'
      real(dp):: q(mxpart,4)
      complex(dp):: pmless(4)
      integer:: h1,h2,h5,h6
      integer:: m,fi,nu,ro,si,i,j
      parameter(m=4)
      complex(dp):: gam0(4,4),gam1(m,4,4),gam2(m,m,4,4),
     & gam3(m,m,m,4,4),gam4(m,m,m,m,4,4),gam5(m,m,m,m,m,4,4),
     & f1(4),f2(4),f5(4),f6(4)
      complex(dp):: p1(4),p2(4),p5(4),p6(4)
      common/gamodd/gam1,gam3,gam5
      common/gameven/gam0,gam2,gam4
      parameter(h1=-1,h2=-1,h6=-1)
      logical:: doJ61,doJ52,docc

      p1(:)=cmplx(q(1,:),kind=dp)    
      p2(:)=cmplx(q(2,:),kind=dp)
c      p3(:)=cmplx(q(3,:),kind=dp)
c      p4(:)=cmplx(q(4,:),kind=dp)
      p5(:)=cmplx(q(5,:),kind=dp)
      p6(:)=cmplx(q(6,:),kind=dp)
      

c--- decide whick strings to comput; note, this is tied
c--- to order of array in qq_tchan_htq_v.f
c--- it relies on h5 = -1,1

      doJ52=.true. 
c--- first call, calculate everything 
      if     (h5 == -1) then
        doJ61=.true.
      elseif (h5 == +1) then
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
      call UbKlSt(p5,mt,pmless,h5,f5)      ! Auxiliary vector pmless
      call uspinor0(p2,h2,f2)
      endif
      
      if (doJ52) then
      J52x0=czip
      do i=1,4
      do j=1,4
      if (abs(gam0(i,j)) > 1d-8) then  
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
      if (abs(gam1(fi,i,j)) > 1d-8) then  
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
      if (abs(gam2(fi,nu,i,j)) > 1d-8) then  
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
      if (abs(gam3(fi,nu,ro,i,j)) > 1d-8) then  
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
      if (abs(gam4(fi,nu,ro,si,i,j)) > 1d-8) then  
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
        J52x0=-h5*conjg(J52x0)
        do fi=1,4
        J52x1(fi)=-h5*conjg(J52x1(fi))
          do nu=1,4
          J52x2(fi,nu)=-h5*conjg(J52x2(fi,nu))
            do ro=1,4
            J52x3(fi,nu,ro)=-h5*conjg(J52x3(fi,nu,ro))
              do si=1,4
              J52x4(fi,nu,ro,si)=-h5*conjg(J52x4(fi,nu,ro,si))
              enddo
            enddo
          enddo
        enddo
        endif
        if (doJ61) then
        do fi=1,4
        J61x1(fi)=conjg(J61x1(fi))
          do nu=1,4
            do ro=1,4
            J61x3(fi,nu,ro)=conjg(J61x3(fi,nu,ro))
            enddo
          enddo
        enddo
        endif
      endif      

      end


c      function vec(p1,ro)
c      implicit none
c      include 'types.f'
c      real(dp):: vec
c      
c      real(dp):: p1(4)
c      integer:: ro
c      vec=p1(ro)
c      return
c      end
