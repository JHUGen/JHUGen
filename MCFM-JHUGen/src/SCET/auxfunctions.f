      function splits2(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      
      real(dp), intent(in) :: z
      real(dp) :: splits2
      
      real(dp) :: fLi2mz,Li2_TR
      !complex(8) :: DiLog,mzcom
      
      !mzcom = -z - ci*0.0_dp  
      
      !fLi2mz = real(dilog(mzcom))
      fLi2mz = Li2_TR(-z)
      
      if(z.eq.1.0_dp) then
       splits2 = 0.0_dp
      ! write(*,*)"splits2:",splits2,"z:",z
       return  
      endif
      
      splits2 = -2*fLi2mz-2*log(1+z)*log(z)-pisq/6.0_dp
      
      !write(*,*)"splits2:",splits2,"z:",z
      
      return      
      end function splits2
      
      
      
      function splits3(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      
      real(dp), intent(in) :: z
      real(dp) :: splits3,splits2
      
      real(dp) :: fLi3z,fLi3omz,fLi3omz2,fLi3ooopz,Li3_TR
      !complex(16) :: xli3,omz,zcom,omz2,ooopz
      
      !zcom = z + 0.0_dp*ci
      !omz = 1.0_dp - zcom
      !omz2 = 1.0_dp - zcom**2
      !ooopz = 1.0_dp/(1.0_dp+zcom)
      
      !fLi3z = real(xLi3_TR(zcom))
      !fLi3omz = real(xLi3_TR(omz))
      !fLi3omz2 = real(xLi3_TR(omz2))
      !fLi3ooopz = real(xLi3_TR(ooopz))
      
      fLi3z = Li3_TR(z)
      fLi3omz = Li3_TR(1.0_dp-z)
      fLi3omz2 = Li3_TR(1.0_dp-z**2)
      fLi3ooopz = Li3_TR(1.0_dp/(1.0_dp+z))
      
      
      if(z.eq.1.0_dp) then
        splits3 = 0.0_dp
      !  write(*,*)"splits3:",splits3,"z:",z
        return
      endif
      
      splits3 = 2*fLi3omz - fLi3z + 4*fLi3ooopz - fLi3omz2
     &        + pisq/3.0_dp*log(1+z) - 2/3.0_dp*log(1+z)**3
     &        - 5*zeta3/2.0_dp + pisq/6.0_dp*log(z)
     &        + splits2(z)*log((1-z)/z) - log(z)**3/4.0_dp
      
      !write(*,*)"splits3:",splits3,"z:",z
      
      return      
      end function splits3
      
      
      function splitT3(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      
      real(dp), intent(in) :: z
      real(dp) :: splitT3
      
      real(dp) :: fLi2z,fLi2omz,fLi3omz,Li2_TR,Li3_TR
      !complex(16) :: xli3,omz,zcom
      !complex(8) :: Dilog
      
      !zcom = z + 0.0_dp*ci
      !omz = 1.0_dp - zcom
      
      !fLi2z = real(Dilog(zcom))
      !fLi2omz = real(Dilog(omz))
      !fLi3omz = real(xLi3_TR(omz))
      
      fLi2z = Li2_TR(z)
      fLi2omz = Li2_TR(1.0_dp - z)
      fLi3omz = Li3_TR(1.0_dp - z)
      
      if(z.eq.1.0_dp) then
       splitT3 = 0.0_dp 
      ! write(*,*)"splitT3:",splitT3,"z:",z
       return
      endif
      
      splitT3 = fLi3omz - fLi2omz*log(1-z)
     &        - log(z)*(fLi2z + log(1-z)**2/2.0_dp
     &                 + 5/12.0_dp*log(z)**2 - pisq/3.0_dp)
      
      !write(*,*)"splitT3:",splitT3,"z:",z
      
      return      
      end function splitT3
      
      
      
      function splitU3(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      
      real(dp), intent(in) :: z
      real(dp) :: splitU3
      
      real(dp) :: fLi2z,fLi2omz,fLi3z,fLi3omz,Li2_TR,Li3_TR
      !complex(16) :: xli3,omz,zcom
      !complex(8) :: Dilog
      
      !zcom = z + 0.0_dp*ci
      !omz = 1.0_dp - zcom
      
      
      !fLi2z = real(Dilog(zcom))
      !fLi2omz = real(Dilog(omz))
      !fLi3z = real(xLi3_TR(zcom))
      !fLi3omz = real(xLi3_TR(omz))
      
      fLi2z = Li2_TR(z)
      fLi2omz = Li2_TR(1.0_dp - z)
      fLi3z = Li3_TR(z)
      fLi3omz = Li3_TR(1.0_dp - z)
      
      if(z.eq.1.0_dp) then
        splitU3 = 0.0_dp
      !  write(*,*)"splitU3:",splitU3,"z:",z
        return
      endif  
      
      splitU3 = -4*fLi3omz + fLi3z - zeta3
     &         -log(1-z)*(fLi2z-pisq/6.0_dp)
     &         +2*fLi2omz*log(z) - log(z)**3/4.0_dp
      
      !write(*,*)"splitU3:",splitU3,"z:",z
      
      return      
      end function splitU3
      
      
      
      function splitV3(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      
      real(dp), intent(in) :: z
      real(dp) :: splitV3
      
      real(dp) :: fLi3z,fLi3omz,Li3_TR
      !complex(16) :: xli3,omz,zcom
      
      
      !zcom = z + 0.0_dp*ci
      !omz = 1.0_dp - zcom
      
      !fLi3z = real(xLi3_TR(zcom))
      !fLi3omz = real(xLi3_TR(omz))
      
      fLi3z = Li3_TR(z)
      fLi3omz = Li3_TR(1.0_dp - z)
      
      if(z.eq.1.0_dp) then
        splitV3 = 0.0_dp
      !  write(*,*)"splitV3:",splitV3,"z:",z
        return
      endif
      
      splitV3 =
     & -4*fLi3omz-5*fLi3z+5*zeta3+log(1-z)*log(z)**2/2.0_dp
     & -(2*log(1-z)**2+11/12.0_dp*log(z)**2-13*pisq/6.0_dp)*log(z)
      
      
      !write(*,*)"splitV3:",splitV3,"z:",z
      !pause
      
      return      
      end function splitV3
      
      
      
      
      subroutine F1G1(f,g,fg)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      real(dp),intent(in) :: f(-1:1), g(-1:1)
      real(dp),intent(out) :: fg(-1:3)
      
      
      fg(-1) = f(-1)*g(-1)-pisq/6.0_dp*f(0)*g(0)
     &        - pisq**2/360.0_dp*f(1)*g(1)
     &        + zeta3*f(1)*g(0) + zeta3*f(0)*g(1)
      
      fg(0) = f(-1)*g(0)+f(0)*g(-1)
     &       - pisq/6.0_dp*f(1)*g(0) - pisq/6.0_dp*f(0)*g(1)
     &       + 2*zeta3*f(1)*g(1)
      
      fg(1) = 2*f(0)*g(0)+f(-1)*g(1)+f(1)*g(-1)-pisq/3.0_dp*f(1)*g(1)
      
      fg(2) = 1.5_dp*(f(1)*g(0)+f(0)*g(1))
      
      fg(3) = f(1)*g(1)
      
      return      
      end subroutine F1G1
      
      
      
!## Li2      
      function Li2_TR(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp) :: z, Li2_TR, Li2
c      complex(8) :: XSPENZ
      include 'cplx.h'
c      Li2_TR = real(XSPENZ(cplx1(z)))

c---- Comparison with the normal MCFM Li2 function
c      write(6,*) 'TR: Li2',Li2_TR
      Li2_TR = Li2(z)
c      write(6,*) 'us: Li2',Li2_TR

      end
    
    
!## Li3    
      function Li3_TR(z)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp) :: z,Li3_TR,Li3
c      complex(8) :: xli3
      include 'cplx.h'
c      Li3_TR = real(xLi3(cplx1(z)))

c---- Comparison with the normal MCFM Li3 function
c      write(6,*) 'TR: Li3',Li3_TR
      Li3_TR = Li3(z)
c      write(6,*) 'us: Li3',Li3_TR

      end
      
      
      
      
