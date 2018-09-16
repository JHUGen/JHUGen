      subroutine strings(q,h3,h5,docc)
      implicit none
      include 'types.f'
! input: momenta q, helicity of 5 is h5
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
      integer:: h1,h2,h3,h4,h5,h6
      integer:: m,fi,nu,ro,si,om,i,j
      parameter(m=4)
      complex(dp):: gam0(4,4),gam1(m,4,4),gam2(m,m,4,4),
     & gam3(m,m,m,4,4),gam4(m,m,m,m,4,4),gam5(m,m,m,m,m,4,4),
     & f1(4),f2(4),f3(4),f4(4),f5(4),f6(4)
      complex(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      real(dp):: q1(4),q2(4),q3(4),q4(4),q5(4),
     & q6(4)
c ,g(4)
c      data g/-1,-1,-1,1/
      common/gamodd/gam1,gam3,gam5
      common/gameven/gam0,gam2,gam4
      parameter(h1=-1,h2=-1,h6=-1)
      logical:: doJ61,doJ52,docc
c      common/current/J34x0,J34x1,J34x3,J34x5,
c     & J52x0,J52x1,J52x3,J52x5,
c     & J61x0,J61x1,J61x3,J61x5

      data q1/0.0000000000000000D+00, 0.0000000000000000D+00,
     &  -0.1058011987414947D+02,-0.1058011987414947D+02/


      data q2/0.0000000000000000D+00, 0.0000000000000000D+00,
     &   0.6850815232034704D+04,-0.6850815232034704D+04/

      data q3/0.5355493692468755D+02, -0.8706764175737293D+02,
     & -0.1729151062441010D+04 ,0.1732169824887580D+04/

      data q4/-0.2302108345379800D+01,-0.1976659016226675D+02,
     &  -0.6411873545470962D+02 ,0.6713590712716950D+02/

      data q5/0.1345849548854817D+02 ,0.2847171395914329D+02,
     &  -0.4611339062715607D+04, 0.4614765944327343D+04/

      data q6/-0.6471132406785591D+02 ,0.7836251796049639D+02,
     &  -0.4356262515492281D+03, 0.4473236755667594D+03/

c      data hltmp/-1,1,-1,-1,-1,-1,-1/
C---statement function
c     vec(p1,ro)=p1(ro)
C---end statement function

c      p1(:)=cmplx(q1(:),kind=dp)
c      p2(:)=cmplx(q2(:),kind=dp)
c      p3(:)=cmplx(q3(:),kind=dp)
c      p4(:)=cmplx(q4(:),kind=dp)
c      p5(:)=cmplx(q5(:),kind=dp)
c      p6(:)=cmplx(q6(:),kind=dp)

      p1(:)=cmplx(q(1,:),kind=dp)
      p2(:)=cmplx(q(2,:),kind=dp)
      p3(:)=cmplx(q(3,:),kind=dp)
      p4(:)=cmplx(q(4,:),kind=dp)
      p5(:)=cmplx(q(5,:),kind=dp)
      p6(:)=cmplx(q(6,:),kind=dp)

c      p346=p3+p4+p6
c      s34=2d0*(p3(4)*p4(4)-p3(1)*p4(1)-p3(2)*p4(2)-p3(3)*p4(3))
c      s346=+p346(4)**2-p346(1)**2-p346(2)**2-p346(3)**2

c--- anti-lepton helicity h4 set by h3
      h4=h3

c--- decide whick strings to comput; note, this is tied
c--- to order of array in qq_tchan_ztq_v.f
c--- it relies on (h3,h5) = (-1,-1), (+1,-1), (-1,+1), (+1,+1)

      doJ61=.false.
      doJ52=.false.
c--- first call, calculate everything
      if ((h3 == -1) .and. (h5 == -1)) then
        doJ61=.true.
	doJ52=.true.
      endif
      if ((h3 == -1) .and. (h5 == +1)) then
	doJ52=.true.
      endif

      if (doJ61) then
      call ubarspinor0(p6,h6,f6)
      call uspinor0(p1,h1,f1)
      endif
      if (doJ52) then
      call UbKlSt(p5,mt,p2,h5,f5)	! Auxiliary vector p2
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
      if (abs(gam0(i,j)) > 1d-8) then
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
      if (abs(gam1(fi,i,j)) > 1d-8) then
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
      if (abs(gam2(fi,nu,i,j)) > 1d-8) then
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
      if (abs(gam3(fi,nu,ro,i,j)) > 1d-8) then
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

      do fi=1,4
      do nu=1,4
      do ro=1,4
      do si=1,4
      do om=1,4
      if (doJ52) J52x5(fi,nu,ro,si,om)=czip
      if (doJ61) J61x5(fi,nu,ro,si,om)=czip
      do i=1,4
      do j=1,4
      if (abs(gam5(fi,nu,ro,si,om,i,j)) > 1d-8) then
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
	  J52x0=-h5*conjg(J52x0)
	  do fi=1,4
	  J52x1(fi)=-h5*conjg(J52x1(fi))
	    do nu=1,4
	    J52x2(fi,nu)=-h5*conjg(J52x2(fi,nu))
	      do ro=1,4
	      J52x3(fi,nu,ro)=-h5*conjg(J52x3(fi,nu,ro))
	        do si=1,4
                J52x4(fi,nu,ro,si)=-h5*conjg(J52x4(fi,nu,ro,si))
                  do om=1,4
                  J52x5(fi,nu,ro,si,om)=-h5*
     &                              conjg(J52x5(fi,nu,ro,si,om))
                  enddo
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
                do si=1,4
	          do om=1,4
	          J61x5(fi,nu,ro,si,om)=conjg(J61x5(fi,nu,ro,si,om))
	          enddo
	        enddo
	      enddo
            enddo
          enddo
        endif
c      do fi=1,4
c         J34x1(fi)=conjg(J34x1(fi))
c         do nu=1,4
c            do ro=1,4
c               J34x3(fi,nu,ro)=conjg(J34x3(fi,nu,ro))
c            enddo
c         enddo
c      enddo
      endif


      end

      function vec(p1,ro)
      implicit none
      include 'types.f'
      real(dp):: vec

      real(dp):: p1(4)
      integer:: ro
      vec=p1(ro)
      return
      end
