      subroutine writeout(p)
      implicit none
      include 'constants.f'
      integer j,n
      double precision p(mxpart,4),dot,sum(4)
      write(6,*) 'In writeout'
      do n=1,8
      do j=1,4
      write(6,'(a8,i3,a,i3,a3,g24.15,a2)')
     & '      p(',n,',',j,') =',p(n,j),'d0'
      enddo
      enddo
      write(6,'(a3,4e24.16)') 'p1',p(1,1),p(1,2),p(1,3),p(1,4)
      write(6,'(a3,4e24.16)') 'p2',p(2,1),p(2,2),p(2,3),p(2,4)
      write(6,'(a3,4e24.16)') 'p3',p(3,1),p(3,2),p(3,3),p(3,4)
      write(6,'(a3,4e24.16)') 'p4',p(4,1),p(4,2),p(4,3),p(4,4)
      write(6,'(a3,4e24.16)') 'p5',p(5,1),p(5,2),p(5,3),p(5,4)
      write(6,'(a3,4e24.16)') 'p6',p(6,1),p(6,2),p(6,3),p(6,4)
      write(6,'(a3,4e24.16)') 'p7',p(7,1),p(7,2),p(7,3),p(7,4)
      write(6,'(a3,4e24.16)') 'p8',p(8,1),p(8,2),p(8,3),p(8,4)
      write(6,'(a3,4e24.16)') 'p9',p(9,1),p(9,2),p(9,3),p(9,4)
      write(6,'(a3,4e24.16)') 'p10',p(10,1),p(10,2),p(10,3),p(10,4)
      write(6,'(a3,4e24.16)') 'p11',p(11,1),p(11,2),p(11,3),p(11,4)
      write(6,'(a3,4e24.16)') 'p12',p(12,1),p(12,2),p(12,3),p(12,4)

      write(6,*) 's12',2d0*dot(p,1,2)
      write(6,*) 's56',2d0*dot(p,5,6)
      write(6,*) 'sqrt(s34)',dsqrt(2d0*dot(p,3,4))
      write(6,*) 'sqrt(s56)',dsqrt(2d0*dot(p,5,6))
      write(6,*) 'sqrt(s67)',dsqrt(2d0*dot(p,6,7))
      write(6,*) 'sqrt(s78)',dsqrt(2d0*dot(p,7,8))
      write(6,*) 'sqrt(s9,10)',dsqrt(2d0*dot(p,9,10))
      write(6,*) 'sqrt(s11,12)',dsqrt(2d0*dot(p,11,12))
      write(6,*) 'sqrt(s345)', 
     .   dsqrt(2d0*dot(p,3,4)+2d0*dot(p,3,5)+2d0*dot(p,4,5))
      write(6,*) 'sqrt(s346)', 
     .   dsqrt(2d0*dot(p,3,4)+2d0*dot(p,3,6)+2d0*dot(p,4,6))
      write(6,*) 'sqrt(s347)', 
     .   dsqrt(2d0*dot(p,3,4)+2d0*dot(p,3,7)+2d0*dot(p,4,7))
      write(6,*) 'sqrt(s567)', 
     .   dsqrt(2d0*dot(p,5,6)+2d0*dot(p,5,7)+2d0*dot(p,6,7))
      write(6,*) 'sqrt(s348)', 
     .   dsqrt(2d0*dot(p,3,4)+2d0*dot(p,3,8)+2d0*dot(p,4,8))
      write(6,*) 'sqrt(s678)', 
     .   dsqrt(2d0*dot(p,6,7)+2d0*dot(p,6,8)+2d0*dot(p,7,8))
      write(6,*) 'sqrt(s349)', 
     .   dsqrt(2d0*dot(p,3,4)+2d0*dot(p,3,9)+2d0*dot(p,4,9))
      write(6,*) 'sqrt(s789)', 
     .   dsqrt(2d0*dot(p,7,8)+2d0*dot(p,7,9)+2d0*dot(p,8,9))
      write(6,*) 'sqrt(s3456)',
     .   dsqrt(2d0*dot(p,3,4)+2d0*dot(p,3,5)+2d0*dot(p,3,6)
     .                                 +2d0*dot(p,4,5)+2d0*dot(p,4,6)
     .                                                +2d0*dot(p,5,6))
      write(6,*) 'sqrt(s3457)',
     .   dsqrt(2d0*dot(p,3,4)+2d0*dot(p,3,5)+2d0*dot(p,3,7)
     .                                 +2d0*dot(p,4,5)+2d0*dot(p,4,7)
     .                                                +2d0*dot(p,5,7))
      write(6,*) 'sqrt(s3458)',
     .   dsqrt(2d0*dot(p,3,4)+2d0*dot(p,3,5)+2d0*dot(p,3,8)
     .                                 +2d0*dot(p,4,5)+2d0*dot(p,4,8)
     .                                                +2d0*dot(p,5,8))
      write(6,*) 'sqrt(s5678)',
     .   dsqrt(2d0*dot(p,5,6)+2d0*dot(p,5,7)+2d0*dot(p,5,8)
     .                                 +2d0*dot(p,6,7)+2d0*dot(p,6,8)
     .                                                +2d0*dot(p,7,8))
      write(6,*) 'sqrt(s3459)',
     .   dsqrt(2d0*dot(p,3,4)+2d0*dot(p,3,5)+2d0*dot(p,3,9)
     .                                 +2d0*dot(p,4,5)+2d0*dot(p,4,9)
     .                                                +2d0*dot(p,5,9))

      write(6,*) 'sqrt(s6789)',
     .   dsqrt(2d0*dot(p,6,7)+2d0*dot(p,6,8)+2d0*dot(p,6,9)
     .                                 +2d0*dot(p,7,8)+2d0*dot(p,7,9)
     .                                                +2d0*dot(p,8,9))

      do j=1,4
      sum(j)=p(1,j)+p(2,j)
      enddo
      do n=3,mxpart
      do j=1,4
      sum(j)=sum(j)+p(n,j)
      enddo
      enddo

      write(6,*) '     psum1',sum(1)
      write(6,*) '     psum2',sum(2)
      write(6,*) '     psum3',sum(3)
      write(6,*) '     psum4',sum(4)
c      do j=1,4
c      sum(j)=-p(1,j)-p(2,j)
c      enddo
c      do n=3,mxpart
c      do j=1,4
c      sum(j)=sum(j)+p(n,j)
c      enddo
c      enddo

c      write(6,*) '     msum1',sum(1)
c      write(6,*) '     msum2',sum(2)
c      write(6,*) '     msum3',sum(3)
c      write(6,*) '     msum4',sum(4)
      write(6,'(a7,e24.16)') 'p1Dp1',dot(p,1,1)
      write(6,'(a7,e24.16)') 'p2Dp2',dot(p,2,2)
      write(6,'(a7,e24.16)') 'p3Dp3',dot(p,3,3)
      write(6,'(a7,e24.16)') 'p4Dp4',dot(p,4,4)
      write(6,'(a7,e24.16)') 'p5Dp5',dot(p,5,5)
      write(6,'(a7,e24.16)') 'p6Dp6',dot(p,6,6)
      write(6,'(a7,e24.16)') 'p7Dp7',dot(p,7,7)
      write(6,'(a7,e24.16)') 'p8Dp8',dot(p,8,8)
      write(6,'(a7,e24.16)') 'p9Dp9',dot(p,9,9)
      write(6,'(a7,e24.16)') 'p10Dp10',dot(p,10,10)
      write(6,'(a7,e24.16)') 'p11Dp11',dot(p,11,11)
      write(6,*)

      call flush(6)
c      pause

      return
      end
