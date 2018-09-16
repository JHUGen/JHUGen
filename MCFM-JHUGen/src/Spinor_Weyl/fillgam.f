      subroutine fillgam
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'swapxz.f'
      integer:: m,mu,nu,ro,si,om,i,j
      parameter(m=4)
      complex(dp):: gam0(4,4),gam1(m,4,4),gam2(m,m,4,4),
     & gam3(m,m,m,4,4),gam4(m,m,m,m,4,4),gam5(m,m,m,m,m,4,4)

      common/gamodd/gam1,gam3,gam5
      common/gameven/gam0,gam2,gam4

      do i=1,4
      do j=1,4
      do mu=1,4
      gam1(mu,i,j)=czip
      enddo
      gam0(i,j)=czip
      enddo
      enddo

c--- gam0 = identity matrix
      gam0(1,1)=+cone
      gam0(2,2)=+cone
      gam0(3,3)=+cone
      gam0(4,4)=+cone

      if (swapxz) then
        gam1(4,1,3)=+cone
        gam1(4,2,4)=+cone
        gam1(4,3,1)=+cone
        gam1(4,4,2)=+cone

        gam1(3,1,4)=-cone
        gam1(3,2,3)=-cone
        gam1(3,3,2)=+cone
        gam1(3,4,1)=+cone

        gam1(2,1,4)=-im
        gam1(2,2,3)=+im
        gam1(2,3,2)=+im
        gam1(2,4,1)=-im

        gam1(1,1,3)=-cone
        gam1(1,2,4)=+cone
        gam1(1,3,1)=+cone
        gam1(1,4,2)=-cone
      else      
        gam1(4,1,3)=+cone
        gam1(4,2,4)=+cone
        gam1(4,3,1)=+cone
        gam1(4,4,2)=+cone

        gam1(1,1,4)=-cone
        gam1(1,2,3)=-cone
        gam1(1,3,2)=+cone
        gam1(1,4,1)=+cone

        gam1(2,1,4)=+im
        gam1(2,2,3)=-im
        gam1(2,3,2)=-im
        gam1(2,4,1)=+im

        gam1(3,1,3)=-cone
        gam1(3,2,4)=+cone
        gam1(3,3,1)=+cone
        gam1(3,4,2)=-cone
      endif

      do mu=1,4
      do nu=1,4
      do i=1,4
      do j=1,4
      gam2(mu,nu,i,j)=
     & +gam1(mu,i,1)*gam1(nu,1,j)
     & +gam1(mu,i,2)*gam1(nu,2,j)
     & +gam1(mu,i,3)*gam1(nu,3,j)
     & +gam1(mu,i,4)*gam1(nu,4,j)
      enddo
      enddo
      enddo
      enddo
     
     
      do mu=1,4
      do nu=1,4
      do ro=1,4
      do i=1,4
      do j=1,4
      gam3(mu,nu,ro,i,j)=
     & +gam1(mu,i,1)*gam2(nu,ro,1,j)
     & +gam1(mu,i,2)*gam2(nu,ro,2,j)
     & +gam1(mu,i,3)*gam2(nu,ro,3,j)
     & +gam1(mu,i,4)*gam2(nu,ro,4,j)
      enddo
      enddo
      enddo
      enddo
      enddo


      do mu=1,4
      do nu=1,4
      do ro=1,4
      do si=1,4
      do i=1,4
      do j=1,4
      gam4(mu,nu,ro,si,i,j)=
     & +gam1(mu,i,1)*gam3(nu,ro,si,1,j)
     & +gam1(mu,i,2)*gam3(nu,ro,si,2,j)
     & +gam1(mu,i,3)*gam3(nu,ro,si,3,j)
     & +gam1(mu,i,4)*gam3(nu,ro,si,4,j)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      do mu=1,4
      do nu=1,4
      do ro=1,4
      do si=1,4
      do om=1,4
      do i=1,4
      do j=1,4
      gam5(mu,nu,ro,si,om,i,j)=
     & +gam1(mu,i,1)*gam4(nu,ro,si,om,1,j)
     & +gam1(mu,i,2)*gam4(nu,ro,si,om,2,j)
     & +gam1(mu,i,3)*gam4(nu,ro,si,om,3,j)
     & +gam1(mu,i,4)*gam4(nu,ro,si,om,4,j)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      end
