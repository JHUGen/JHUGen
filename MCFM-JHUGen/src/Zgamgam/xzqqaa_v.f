**************************************************************** 
*   Color ordered virtual amplitudes for:
*   0 -> q(-p1) + qb(-p4) + a(p2) + a(p3) + lb(p5) + l(p6)
****************************************************************
      subroutine zaa_a6vh(j1,j2,j3,j4,j5,j6,za,zb,a6vh)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'ipsgen.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a6g,a6vQLslc,a6virtLL
      complex(dp):: a6vh(4,16)
      character*9 st1,st2,st3
      character*14 st4,st5 
      integer:: i,j,k
c-----helicity stamp
      st1='q+qb-g+g+'
      st2='q+qb-g+g-'
      st3='q+qb-g-g+'
      st4='q+qb-g+g+lb-l+'
      st5='q+qb-g+g-lb-l+'
c-----initialize a60h
      do i=1,4
      do j=1,8
      a6vh(i,j)=czip
      enddo
      enddo
c-----QQ contributions
      if (ipsgen .ne. 2) then
      a6vh(1,1) = a6g(st1,j1,j4,j2,j3,j5,j6,za,zb)
     &           +a6g(st1,j1,j4,j3,j2,j5,j6,za,zb)
      a6vh(1,9) = a6g(st1,j1,j4,j2,j3,j5,j6,zb,za)
     &           +a6g(st1,j1,j4,j3,j2,j5,j6,zb,za)
      a6vh(1,2) = a6g(st2,j1,j4,j2,j3,j5,j6,za,zb)
     &           +a6g(st3,j1,j4,j3,j2,j5,j6,za,zb)
      a6vh(1,10)= a6g(st2,j1,j4,j2,j3,j5,j6,zb,za)
     &           +a6g(st3,j1,j4,j3,j2,j5,j6,zb,za)
      a6vh(1,3) = a6g(st3,j1,j4,j2,j3,j5,j6,za,zb)
     &           +a6g(st2,j1,j4,j3,j2,j5,j6,za,zb)
      a6vh(1,11)= a6g(st3,j1,j4,j2,j3,j5,j6,zb,za)
     &           +a6g(st2,j1,j4,j3,j2,j5,j6,zb,za)
      a6vh(1,4) = a6g(st1,j4,j1,j2,j3,j6,j5,zb,za)
     &           +a6g(st1,j4,j1,j3,j2,j6,j5,zb,za)
      a6vh(1,12)= a6g(st1,j4,j1,j2,j3,j6,j5,za,zb)
     &           +a6g(st1,j4,j1,j3,j2,j6,j5,za,zb)
      a6vh(1,5) = a6g(st1,j1,j4,j2,j3,j6,j5,za,zb)
     &           +a6g(st1,j1,j4,j3,j2,j6,j5,za,zb)
      a6vh(1,13)= a6g(st1,j1,j4,j2,j3,j6,j5,zb,za)
     &           +a6g(st1,j1,j4,j3,j2,j6,j5,zb,za)
      a6vh(1,6) = a6g(st2,j1,j4,j2,j3,j6,j5,za,zb)
     &           +a6g(st3,j1,j4,j3,j2,j6,j5,za,zb)
      a6vh(1,14)= a6g(st2,j1,j4,j2,j3,j6,j5,zb,za)
     &           +a6g(st3,j1,j4,j3,j2,j6,j5,zb,za)
      a6vh(1,7) = a6g(st3,j1,j4,j2,j3,j6,j5,za,zb)
     &           +a6g(st2,j1,j4,j3,j2,j6,j5,za,zb)
      a6vh(1,15)= a6g(st3,j1,j4,j2,j3,j6,j5,zb,za)
     &           +a6g(st2,j1,j4,j3,j2,j6,j5,zb,za)
      a6vh(1,8) = a6g(st1,j4,j1,j2,j3,j5,j6,zb,za)
     &           +a6g(st1,j4,j1,j3,j2,j5,j6,zb,za)
      a6vh(1,16)= a6g(st1,j4,j1,j2,j3,j5,j6,za,zb)
     &           +a6g(st1,j4,j1,j3,j2,j5,j6,za,zb)
      endif
c-----LL contributions
      if (ipsgen .ne. 3) then
      a6vh(2,1) = a6virtLL(st1,j6,j5,j2,j3,j4,j1,za,zb)
     &           +a6virtLL(st1,j6,j5,j3,j2,j4,j1,za,zb)
      a6vh(2,9) = a6virtLL(st1,j6,j5,j2,j3,j4,j1,zb,za)
     &           +a6virtLL(st1,j6,j5,j3,j2,j4,j1,zb,za)
      a6vh(2,2) = a6virtLL(st2,j6,j5,j2,j3,j4,j1,za,zb)
     &           +a6virtLL(st3,j6,j5,j3,j2,j4,j1,za,zb)
      a6vh(2,10)= a6virtLL(st2,j6,j5,j2,j3,j4,j1,zb,za)
     &           +a6virtLL(st3,j6,j5,j3,j2,j4,j1,zb,za)
      a6vh(2,3) = a6virtLL(st3,j6,j5,j2,j3,j4,j1,za,zb)
     &           +a6virtLL(st2,j6,j5,j3,j2,j4,j1,za,zb)
      a6vh(2,11)= a6virtLL(st3,j6,j5,j2,j3,j4,j1,zb,za)
     &           +a6virtLL(st2,j6,j5,j3,j2,j4,j1,zb,za)
      a6vh(2,4) = a6virtLL(st1,j5,j6,j2,j3,j1,j4,zb,za)
     &           +a6virtLL(st1,j5,j6,j3,j2,j1,j4,zb,za)
      a6vh(2,12)= a6virtLL(st1,j5,j6,j2,j3,j1,j4,za,zb)
     &           +a6virtLL(st1,j5,j6,j3,j2,j1,j4,za,zb)
      a6vh(2,5) = a6virtLL(st1,j5,j6,j2,j3,j4,j1,za,zb)
     &           +a6virtLL(st1,j5,j6,j3,j2,j4,j1,za,zb)
      a6vh(2,13)= a6virtLL(st1,j5,j6,j2,j3,j4,j1,zb,za)
     &           +a6virtLL(st1,j5,j6,j3,j2,j4,j1,zb,za)
      a6vh(2,6) = a6virtLL(st2,j5,j6,j2,j3,j4,j1,za,zb)
     &           +a6virtLL(st3,j5,j6,j3,j2,j4,j1,za,zb)
      a6vh(2,14)= a6virtLL(st2,j5,j6,j2,j3,j4,j1,zb,za)
     &           +a6virtLL(st3,j5,j6,j3,j2,j4,j1,zb,za)
      a6vh(2,7) = a6virtLL(st3,j5,j6,j2,j3,j4,j1,za,zb)
     &           +a6virtLL(st2,j5,j6,j3,j2,j4,j1,za,zb)
      a6vh(2,15)= a6virtLL(st3,j5,j6,j2,j3,j4,j1,zb,za)
     &           +a6virtLL(st2,j5,j6,j3,j2,j4,j1,zb,za)
      a6vh(2,8) = a6virtLL(st1,j6,j5,j2,j3,j1,j4,zb,za)
     &           +a6virtLL(st1,j6,j5,j3,j2,j1,j4,zb,za)
      a6vh(2,16)= a6virtLL(st1,j6,j5,j2,j3,j1,j4,za,zb)
     &           +a6virtLL(st1,j6,j5,j3,j2,j1,j4,za,zb)
      endif
c-----QL contributions
      if ((ipsgen == 2) .or. (ipsgen == 3)) then
      a6vh(3,1) = -a6vQLslc(st4,j1,j4,j2,j3,j5,j6,za,zb)
      a6vh(3,9) = -a6vQLslc(st4,j1,j4,j2,j3,j5,j6,zb,za)
      a6vh(3,2) = -a6vQLslc(st5,j1,j4,j2,j3,j5,j6,za,zb)
      a6vh(3,10)= -a6vQLslc(st5,j1,j4,j2,j3,j5,j6,zb,za)
      a6vh(3,3) = -a6vQLslc(st5,j4,j1,j2,j3,j6,j5,zb,za)
      a6vh(3,11)= -a6vQLslc(st5,j4,j1,j2,j3,j6,j5,za,zb)
      a6vh(3,4) = -a6vQLslc(st4,j4,j1,j2,j3,j6,j5,zb,za)
      a6vh(3,12)= -a6vQLslc(st4,j4,j1,j2,j3,j6,j5,za,zb)
      a6vh(3,5) = -a6vQLslc(st4,j1,j4,j2,j3,j6,j5,za,zb)
      a6vh(3,13)= -a6vQLslc(st4,j1,j4,j2,j3,j6,j5,zb,za)
      a6vh(3,6) = -a6vQLslc(st5,j1,j4,j2,j3,j6,j5,za,zb)
      a6vh(3,14)= -a6vQLslc(st5,j1,j4,j2,j3,j6,j5,zb,za)
      a6vh(3,7) = -a6vQLslc(st5,j4,j1,j2,j3,j5,j6,zb,za)
      a6vh(3,15)= -a6vQLslc(st5,j4,j1,j2,j3,j5,j6,za,zb)
      a6vh(3,8) = -a6vQLslc(st4,j4,j1,j2,j3,j5,j6,zb,za)
      a6vh(3,16)= -a6vQLslc(st4,j4,j1,j2,j3,j5,j6,za,zb)
      endif
      
      if ((ipsgen == 3) .or. (ipsgen == 4)) then
      a6vh(4,1) = -a6vQLslc(st4,j1,j4,j3,j2,j5,j6,za,zb)
      a6vh(4,9) = -a6vQLslc(st4,j1,j4,j3,j2,j5,j6,zb,za)
      a6vh(4,2) = -a6vQLslc(st5,j4,j1,j3,j2,j6,j5,zb,za)
      a6vh(4,10)= -a6vQLslc(st5,j4,j1,j3,j2,j6,j5,zb,za)
      a6vh(4,3) = -a6vQLslc(st5,j1,j4,j3,j2,j5,j6,za,zb)
      a6vh(4,11)= -a6vQLslc(st5,j1,j4,j3,j2,j5,j6,zb,za)
      a6vh(4,4) = -a6vQLslc(st4,j4,j1,j3,j2,j6,j5,zb,za)
      a6vh(4,12)= -a6vQLslc(st4,j4,j1,j3,j2,j6,j5,za,zb)
      a6vh(4,5) = -a6vQLslc(st4,j1,j4,j3,j2,j6,j5,za,zb)
      a6vh(4,13)= -a6vQLslc(st4,j1,j4,j3,j2,j6,j5,zb,za)
      a6vh(4,6) = -a6vQLslc(st5,j4,j1,j3,j2,j5,j6,zb,za)
      a6vh(4,14)= -a6vQLslc(st5,j4,j1,j3,j2,j5,j6,za,zb)
      a6vh(4,7) = -a6vQLslc(st5,j1,j4,j3,j2,j6,j5,za,zb)
      a6vh(4,15)= -a6vQLslc(st5,j1,j4,j3,j2,j6,j5,zb,za)
      a6vh(4,8) = -a6vQLslc(st4,j4,j1,j3,j2,j5,j6,zb,za)
      a6vh(4,16)= -a6vQLslc(st4,j4,j1,j3,j2,j5,j6,za,zb)
      endif
      
c-----
c      do i=1,4
c      do j=1,8
c         write(*,*) a6vh(i,j) 
c      enddo
c      enddo
c      stop
c-----done here  

      return
      end
