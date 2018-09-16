      subroutine kingen(p)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      integer:: iout
      real(dp):: p(mxpart,4),s(mxpart,mxpart)
      call dotem(10,p,s)
      write(6,*) 'p1',p(1,1),p(1,2),p(1,3),p(1,4)
      write(6,*) 'p2',p(2,1),p(2,2),p(2,3),p(2,4)
      write(6,*) 'p3',p(3,1),p(3,2),p(3,3),p(3,4)
      write(6,*) 'p4',p(4,1),p(4,2),p(4,3),p(4,4)
      write(6,*) 'p5',p(5,1),p(5,2),p(5,3),p(5,4)
      write(6,*) 'p6',p(6,1),p(6,2),p(6,3),p(6,4)
      write(6,*) 'p7',p(7,1),p(7,2),p(7,3),p(7,4)
      write(6,*) 'p8',p(8,1),p(8,2),p(8,3),p(8,4)
      write(6,*) 'p9',p(9,1),p(9,2),p(9,3),p(9,4)
      write(6,*) 'p0',p(10,1),p(10,2),p(10,3),p(10,4)

      open(unit=7,file='kin.blo',status='unknown')
      do iout=7,7
      write(iout,*)'P ninput'
      write(iout,*)'Id,Numer,p1Dp1,',0
      write(iout,*)'Al,Numer,p1Dp2,',half*s(1,2)
      write(iout,*)'Al,Numer,p1Dp3,',half*s(1,3)
      write(iout,*)'Al,Numer,p1Dp4,',half*s(1,4)
      write(iout,*)'Al,Numer,p1Dp5,',half*s(1,5)
      write(iout,*)'Al,Numer,p1Dp6,',half*s(1,6)
      write(iout,*)'Al,Numer,p1Dp7,',half*s(1,7)
      write(iout,*)'Al,Numer,p1Dp8,',half*s(1,8)
      write(iout,*)'Al,Numer,p1Dp9,',half*s(1,9)
      write(iout,*)'Al,Numer,p1Dp0,',half*s(1,10)
      write(iout,*)'*yep'
      write(iout,*)'Id,Numer,p2Dp2,',0
      write(iout,*)'Al,Numer,p2Dp3,',half*s(2,3)
      write(iout,*)'Al,Numer,p2Dp4,',half*s(2,4)
      write(iout,*)'Al,Numer,p2Dp5,',half*s(2,5)
      write(iout,*)'Al,Numer,p2Dp6,',half*s(2,6)
      write(iout,*)'Al,Numer,p2Dp7,',half*s(2,7)
      write(iout,*)'Al,Numer,p2Dp8,',half*s(2,8)
      write(iout,*)'Al,Numer,p2Dp9,',half*s(2,9)
      write(iout,*)'Al,Numer,p2Dp0,',half*s(2,10)
      write(iout,*)'*yep'
      write(iout,*)'Id,Numer,p3Dp3,',0
      write(iout,*)'Al,Numer,p3Dp4,',half*s(3,4)
      write(iout,*)'Al,Numer,p3Dp5,',half*s(3,5)
      write(iout,*)'Al,Numer,p3Dp6,',half*s(3,6)
      write(iout,*)'Al,Numer,p3Dp7,',half*s(3,7)
      write(iout,*)'Al,Numer,p3Dp8,',half*s(3,8)
      write(iout,*)'Al,Numer,p3Dp9,',half*s(3,9)
      write(iout,*)'Al,Numer,p3Dp0,',half*s(3,10)
      write(iout,*)'*yep'
      write(iout,*)'Id,Numer,p4Dp4,',0
      write(iout,*)'Al,Numer,p4Dp5,',half*s(4,5)
      write(iout,*)'Al,Numer,p4Dp6,',half*s(4,6)
      write(iout,*)'Al,Numer,p4Dp7,',half*s(4,7)
      write(iout,*)'Al,Numer,p4Dp8,',half*s(4,8)
      write(iout,*)'Al,Numer,p4Dp9,',half*s(4,9)
      write(iout,*)'Al,Numer,p4Dp0,',half*s(4,10)
      write(iout,*)'*yep'
      write(iout,*)'Id,Numer,p5Dp5,',0
      write(iout,*)'Al,Numer,p5Dp6,',half*s(5,6)
      write(iout,*)'Al,Numer,p5Dp7,',half*s(5,7)
      write(iout,*)'Al,Numer,p5Dp8,',half*s(5,8)
      write(iout,*)'Al,Numer,p5Dp9,',half*s(5,9)
      write(iout,*)'Al,Numer,p5Dp0,',half*s(5,10)
      write(iout,*)'*yep'
      write(iout,*)'Id,Numer,p6Dp6,',0
      write(iout,*)'Al,Numer,p6Dp7,',half*s(6,7)
      write(iout,*)'Al,Numer,p6Dp8,',half*s(6,8)
      write(iout,*)'Al,Numer,p6Dp9,',half*s(6,9)
      write(iout,*)'Al,Numer,p6Dp0,',half*s(6,10)
      write(iout,*)'*yep'
      write(iout,*)'Id,Numer,p7Dp7,',0
      write(iout,*)'Al,Numer,p7Dp8,',half*s(7,8)
      write(iout,*)'Al,Numer,p7Dp9,',half*s(7,9)
      write(iout,*)'Al,Numer,p7Dp0,',half*s(7,10)
      write(iout,*)'*yep'
      write(iout,*)'Id,Numer,p8Dp8,',0
      write(iout,*)'Al,Numer,p8Dp9,',half*s(8,9)
      write(iout,*)'Al,Numer,p8Dp0,',half*s(8,10)
      write(iout,*)'*yep'
      write(iout,*)'Id,Numer,p9Dp9,',0
      write(iout,*)'Al,Numer,p9Dp0,',half*s(9,10)
      write(iout,*)'*yep'
      write(iout,*)'Id,Numer,p0Dp0,',0
      write(iout,*)'P input'
      write(iout,*)'End'
      enddo
C q1=p2+p3+p4+p5
C r1=p2+p3+p6+p7
C q2=p2+p4+p5
C r2=p2+p6+p7
C t1=p3+p4+t5
C t3=p4+p5
C t4=p3+p5
C t5=p3+t4

      close(unit=7)
      return
      end
