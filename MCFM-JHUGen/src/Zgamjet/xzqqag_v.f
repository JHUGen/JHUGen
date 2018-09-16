**************************************************************** 
*   Color ordered virtual amplitudes for:
*   0 -> q(-p1) + qb(-p4) + g(p2) + a(p3) + lb(p5) + l(p6)
****************************************************************
      subroutine zaj_a6vh(j1,j2,j3,j4,j5,j6,za,zb,a6vh)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a6g,a64v,a64ax
      complex(dp):: a6vQLlc,a6vQLslc,a6vQLfloop
      complex(dp):: a6vh(5,16)
      character*9 st1,st2,st3,st4,st5
      character*14 st6,st7 
      integer:: i,j,k
c-----helicity stamp
      st1='q+g+qb-g+'
      st2='q+g+qb-g-'
      st3='q+qb-g+g+'
      st4='q+qb-g+g-'
      st5='q+qb-g-g+'
      st6='q+qb-g+g+lb-l+'
      st7='q+qb-g+g-lb-l+'
c-----initialize a60h
      do i=1,5
      do j=1,16
      a6vh(i,j)=czip
      enddo
      enddo
c-----QQ contributions
c-----photon is coming from quark line
c-----contribution from: a6g,a6v,a6ax
c-----1:a6g
c-----(1q+2g+3a+4qb-5lb-6l+)
      a6vh(1,1) = -xn*a6g(st1,j1,j2,j4,j3,j5,j6,za,zb)
     &            -(1._dp/xn)*( a6g(st3,j1,j4,j2,j3,j5,j6,za,zb)
     &                       +a6g(st3,j1,j4,j3,j2,j5,j6,za,zb) )
c-----(1q+2g+3a-4qb-5lb-6l+)
      a6vh(1,2) = -xn*a6g(st2,j1,j2,j4,j3,j5,j6,za,zb)
     &            -(1._dp/xn)*( a6g(st4,j1,j4,j2,j3,j5,j6,za,zb)
     &                       +a6g(st5,j1,j4,j3,j2,j5,j6,za,zb) )
c-----(1q+2g-3a+4qb-5lb-6l+)
      a6vh(1,3) = -xn*a6g(st2,j4,j2,j1,j3,j6,j5,zb,za)
     &            -(1._dp/xn)*( a6g(st5,j1,j4,j2,j3,j5,j6,za,zb)
     &                       +a6g(st4,j1,j4,j3,j2,j5,j6,za,zb) )
c-----(1q+2g-3a-4qb-5lb-6l+)
      a6vh(1,4) = -xn*a6g(st1,j4,j2,j1,j3,j6,j5,zb,za)
     &            -(1._dp/xn)*( a6g(st3,j4,j1,j2,j3,j6,j5,zb,za)
     &                       +a6g(st3,j4,j1,j3,j2,j6,j5,zb,za) )
c-----(1q+2g+3a+4qb-5lb+6l-)
      a6vh(1,5) = -xn*a6g(st1,j1,j2,j4,j3,j6,j5,za,zb)
     &            -(1._dp/xn)*( a6g(st3,j1,j4,j2,j3,j6,j5,za,zb)
     &                       +a6g(st3,j1,j4,j3,j2,j6,j5,za,zb) )
c-----(1q+g+3a-4qb-lb+6l-)
      a6vh(1,6) = -xn*a6g(st2,j1,j2,j4,j3,j6,j5,za,zb)
     &            -(1._dp/xn)*( a6g(st4,j1,j4,j2,j3,j6,j5,za,zb)
     &                       +a6g(st5,j1,j4,j3,j2,j6,j5,za,zb) )
c-----(1q+2g-3a+4qb-5lb+6l-)
      a6vh(1,7) = -xn*a6g(st2,j4,j2,j1,j3,j5,j6,zb,za)
     &            -(1._dp/xn)*( a6g(st5,j1,j4,j2,j3,j6,j5,za,zb)
     &                       +a6g(st4,j1,j4,j3,j2,j6,j5,za,zb) )
c-----(1q+2g-3a-4qb-5lb+6l-)
      a6vh(1,8) = -xn*a6g(st1,j4,j2,j1,j3,j5,j6,zb,za)
     &            -(1._dp/xn)*( a6g(st3,j4,j1,j2,j3,j5,j6,zb,za)
     &                       +a6g(st3,j4,j1,j3,j2,j5,j6,zb,za) )
c-----(1q-2g-3a-4qb+5lb+6l-)
      a6vh(1,9)  = -xn*a6g(st1,j1,j2,j4,j3,j5,j6,zb,za)
     &             -(1._dp/xn)*( a6g(st3,j1,j4,j2,j3,j5,j6,zb,za)
     &                        +a6g(st3,j1,j4,j3,j2,j5,j6,zb,za) )
c-----(1q-2g-3a+4qb+5lb+6l-)
      a6vh(1,10) = -xn*a6g(st2,j1,j2,j4,j3,j5,j6,zb,za)
     &             -(1._dp/xn)*( a6g(st4,j1,j4,j2,j3,j5,j6,zb,za)
     &                        +a6g(st5,j1,j4,j3,j2,j5,j6,zb,za) )
c-----(1q-g+3a-4qb+lb+6l-)
      a6vh(1,11) = -xn*a6g(st2,j4,j2,j1,j3,j6,j5,za,zb)
     &             -(1._dp/xn)*( a6g(st5,j1,j4,j2,j3,j5,j6,zb,za)
     &                        +a6g(st4,j1,j4,j3,j2,j5,j6,zb,za) )
c-----(1q-2g+3a+4qb+5lb+6l-)
      a6vh(1,12) = -xn*a6g(st1,j4,j2,j1,j3,j6,j5,za,zb)
     &             -(1._dp/xn)*( a6g(st3,j4,j1,j2,j3,j6,j5,za,zb)
     &                        +a6g(st3,j4,j1,j3,j2,j6,j5,za,zb) )
c-----(1q-2g-3a-4qb+5lb-6l+)
      a6vh(1,13) = -xn*a6g(st1,j1,j2,j4,j3,j6,j5,zb,za)
     &             -(1._dp/xn)*( a6g(st3,j1,j4,j2,j3,j6,j5,zb,za)
     &                        +a6g(st3,j1,j4,j3,j2,j6,j5,zb,za) )
c-----(1q-2g-3a+4qb+5lb-6l+)
      a6vh(1,14) = -xn*a6g(st2,j1,j2,j4,j3,j6,j5,zb,za)
     &             -(1._dp/xn)*( a6g(st4,j1,j4,j2,j3,j6,j5,zb,za)
     &                        +a6g(st5,j1,j4,j3,j2,j6,j5,zb,za) )
c-----(1q-2g+3a-4qb+5lb-6l+)
      a6vh(1,15) = -xn*a6g(st2,j4,j2,j1,j3,j5,j6,za,zb)
     &             -(1._dp/xn)*( a6g(st5,j1,j4,j2,j3,j6,j5,zb,za)
     &                        +a6g(st4,j1,j4,j3,j2,j6,j5,zb,za) )
c-----(1q-2g+3a+4qb+5lb-6l+)
      a6vh(1,16) = -xn*a6g(st1,j4,j2,j1,j3,j5,j6,za,zb)
     &             -(1._dp/xn)*( a6g(st3,j4,j1,j2,j3,j5,j6,za,zb)
     &                        +a6g(st3,j4,j1,j3,j2,j5,j6,za,zb) )
c-----2:a6v
c-----(1q+2g+3a+4qb-5lb-6l+)
      a6vh(2,1)  =  two*a64v(st3,j1,j4,j2,j3,j5,j6,za,zb)
      a6vh(2,9)  =  two*a64v(st3,j1,j4,j2,j3,j5,j6,zb,za)
c-----(1q+2g+3a-4qb-5lb-6l+)
      a6vh(2,2)  =  two*a64v(st4,j1,j4,j2,j3,j5,j6,za,zb)
      a6vh(2,10) =  two*a64v(st4,j1,j4,j2,j3,j5,j6,zb,za)
c-----(1q+2g-3a+4qb-5lb-6l+)
      a6vh(2,3)  =  two*a64v(st5,j1,j4,j2,j3,j5,j6,za,zb)
      a6vh(2,11) =  two*a64v(st5,j1,j4,j2,j3,j5,j6,zb,za)
c-----(1q+2g-3a-4qb-5lb-6l+)
      a6vh(2,4)  =  two*a64v(st3,j4,j1,j2,j3,j6,j5,zb,za)
      a6vh(2,12) =  two*a64v(st3,j4,j1,j2,j3,j6,j5,za,zb)
c-----(1q+2g+3a+4qb-5lb+6l-)
      a6vh(2,5)  =  two*a64v(st3,j1,j4,j2,j3,j6,j5,za,zb)
      a6vh(2,13) =  two*a64v(st3,j1,j4,j2,j3,j6,j5,zb,za)
c-----(1q+2g+3a-4qb-5lb+6l-)
      a6vh(2,6)  =  two*a64v(st4,j1,j4,j2,j3,j6,j5,za,zb)
      a6vh(2,14) =  two*a64v(st4,j1,j4,j2,j3,j6,j5,zb,za)
c-----(1q+2g-3a+4qb-5lb+6l-)
      a6vh(2,7)  =  two*a64v(st5,j1,j4,j2,j3,j6,j5,za,zb)
      a6vh(2,15) =  two*a64v(st5,j1,j4,j2,j3,j6,j5,zb,za)
c-----(1q+2g-3a-4qb-5lb+6l-)
      a6vh(2,8)  =  two*a64v(st3,j4,j1,j2,j3,j5,j6,zb,za)
      a6vh(2,16) =  two*a64v(st3,j4,j1,j2,j3,j5,j6,za,zb)
c-----3:a6ax.
c-----(1q+2g+3a+4qb-5lb-6l+)
      a6vh(3,1)  =  a64ax(st3,j1,j4,j2,j3,j5,j6,za,zb)
     &             +a64ax(st3,j1,j4,j3,j2,j5,j6,za,zb) 
      a6vh(3,9)  =  a64ax(st3,j1,j4,j2,j3,j5,j6,zb,za)
     &             +a64ax(st3,j1,j4,j3,j2,j5,j6,zb,za) 
c-----(1q+2g+3a-4qb-5lb-6l+)
      a6vh(3,2)  =  a64ax(st4,j1,j4,j2,j3,j5,j6,za,zb)
     &             +a64ax(st5,j1,j4,j3,j2,j5,j6,za,zb) 
      a6vh(3,10) =  a64ax(st4,j1,j4,j2,j3,j5,j6,zb,za)
     &             +a64ax(st5,j1,j4,j3,j2,j5,j6,zb,za) 
c-----(1q+2g-3a+4qb-5lb-6l+)
      a6vh(3,3)  =  a64ax(st5,j1,j4,j2,j3,j5,j6,za,zb)
     &             +a64ax(st4,j1,j4,j3,j2,j5,j6,za,zb) 
      a6vh(3,11) =  a64ax(st5,j1,j4,j2,j3,j5,j6,zb,za)
     &             +a64ax(st4,j1,j4,j3,j2,j5,j6,zb,za) 
c-----(1q+2g-3a-4qb-5lb-6l+)
      a6vh(3,4)  =  a64ax(st3,j4,j1,j2,j3,j6,j5,zb,za)
     &             +a64ax(st3,j4,j1,j3,j2,j6,j5,zb,za) 
      a6vh(3,12) =  a64ax(st3,j4,j1,j2,j3,j6,j5,za,zb)
     &             +a64ax(st3,j4,j1,j3,j2,j6,j5,za,zb) 
c-----(1q+2g+3a+4qb-5lb+6l-)
      a6vh(3,5)  =  a64ax(st3,j1,j4,j2,j3,j6,j5,za,zb)
     &             +a64ax(st3,j1,j4,j3,j2,j6,j5,za,zb) 
      a6vh(3,13) =  a64ax(st3,j1,j4,j2,j3,j6,j5,zb,za)
     &             +a64ax(st3,j1,j4,j3,j2,j6,j5,zb,za) 
c-----(1q+2g+3a-4qb-5lb+6l-)
      a6vh(3,6)  =  a64ax(st4,j1,j4,j2,j3,j6,j5,za,zb)
     &             +a64ax(st5,j1,j4,j3,j2,j6,j5,za,zb) 
      a6vh(3,14) =  a64ax(st4,j1,j4,j2,j3,j6,j5,zb,za)
     &             +a64ax(st5,j1,j4,j3,j2,j6,j5,zb,za) 
c-----(1q+2g-3a+4qb-5lb+6l-)
      a6vh(3,7)  =  a64ax(st5,j1,j4,j2,j3,j6,j5,za,zb)
     &             +a64ax(st4,j1,j4,j3,j2,j6,j5,za,zb)
      a6vh(3,15) =  a64ax(st5,j1,j4,j2,j3,j6,j5,zb,za)
     &             +a64ax(st4,j1,j4,j3,j2,j6,j5,zb,za)
c-----(1q+2g-3a-4qb-5lb+6l-)
      a6vh(3,8)  =  a64ax(st3,j4,j1,j2,j3,j5,j6,zb,za)
     &             +a64ax(st3,j4,j1,j3,j2,j5,j6,zb,za) 
      a6vh(3,16) =  a64ax(st3,j4,j1,j2,j3,j5,j6,za,zb)
     &             +a64ax(st3,j4,j1,j3,j2,j5,j6,za,zb) 

c-----QL contributions
c-----photon is coming from lepton line
c-----4:a6lc+a6slc
      a6vh(4,1) = +xn*a6vQLlc(st6,j1,j4,j2,j3,j5,j6,za,zb)
     &            +(1._dp/xn)*a6vQLslc(st6,j1,j4,j2,j3,j5,j6,za,zb)
      a6vh(4,9) = +xn*a6vQLlc(st6,j1,j4,j2,j3,j5,j6,zb,za)
     &            +(1._dp/xn)*a6vQLslc(st6,j1,j4,j2,j3,j5,j6,zb,za)
c-----
      a6vh(4,2) = +xn*a6vQLlc(st7,j1,j4,j2,j3,j5,j6,za,zb)
     &            +(1._dp/xn)*a6vQLslc(st7,j1,j4,j2,j3,j5,j6,za,zb)
      a6vh(4,10)= +xn*a6vQLlc(st7,j1,j4,j2,j3,j5,j6,zb,za)
     &            +(1._dp/xn)*a6vQLslc(st7,j1,j4,j2,j3,j5,j6,zb,za)
c-----
      a6vh(4,3) = +xn*a6vQLlc(st7,j4,j1,j2,j3,j6,j5,zb,za)
     &            +(1._dp/xn)*a6vQLslc(st7,j4,j1,j2,j3,j6,j5,zb,za)
      a6vh(4,11)= +xn*a6vQLlc(st7,j4,j1,j2,j3,j6,j5,za,zb)
     &            +(1._dp/xn)*a6vQLslc(st7,j4,j1,j2,j3,j6,j5,za,zb)
c-----
      a6vh(4,4) = +xn*a6vQLlc(st6,j4,j1,j2,j3,j6,j5,zb,za)
     &            +(1._dp/xn)*a6vQLslc(st6,j4,j1,j2,j3,j6,j5,zb,za)
      a6vh(4,12)= +xn*a6vQLlc(st6,j4,j1,j2,j3,j6,j5,za,zb)
     &            +(1._dp/xn)*a6vQLslc(st6,j4,j1,j2,j3,j6,j5,za,zb)
c-----
      a6vh(4,5) = +xn*a6vQLlc(st6,j1,j4,j2,j3,j6,j5,za,zb)
     &            +(1._dp/xn)*a6vQLslc(st6,j1,j4,j2,j3,j6,j5,za,zb)
      a6vh(4,13)= +xn*a6vQLlc(st6,j1,j4,j2,j3,j6,j5,zb,za)
     &            +(1._dp/xn)*a6vQLslc(st6,j1,j4,j2,j3,j6,j5,zb,za)
c-----
      a6vh(4,6) = +xn*a6vQLlc(st7,j1,j4,j2,j3,j6,j5,za,zb)
     &            +(1._dp/xn)*a6vQLslc(st7,j1,j4,j2,j3,j6,j5,za,zb)
      a6vh(4,14)= +xn*a6vQLlc(st7,j1,j4,j2,j3,j6,j5,zb,za)
     &            +(1._dp/xn)*a6vQLslc(st7,j1,j4,j2,j3,j6,j5,zb,za)
c-----
      a6vh(4,7) = +xn*a6vQLlc(st7,j4,j1,j2,j3,j5,j6,zb,za)
     &            +(1._dp/xn)*a6vQLslc(st7,j4,j1,j2,j3,j5,j6,zb,za)
      a6vh(4,15)= +xn*a6vQLlc(st7,j4,j1,j2,j3,j5,j6,za,zb)
     &            +(1._dp/xn)*a6vQLslc(st7,j4,j1,j2,j3,j5,j6,za,zb)
c-----
      a6vh(4,8) = +xn*a6vQLlc(st6,j4,j1,j2,j3,j5,j6,zb,za)
     &            +(1._dp/xn)*a6vQLslc(st6,j4,j1,j2,j3,j5,j6,zb,za)
      a6vh(4,16)= +xn*a6vQLlc(st6,j4,j1,j2,j3,j5,j6,za,zb)
     &            +(1._dp/xn)*a6vQLslc(st6,j4,j1,j2,j3,j5,j6,za,zb)
c-----5:a6floop
      a6vh(5,1) = a6vQLfloop(st6,j1,j4,j2,j3,j5,j6,za,zb)
      a6vh(5,9) = a6vQLfloop(st6,j1,j4,j2,j3,j5,j6,zb,za)
c-----
      a6vh(5,2) = a6vQLfloop(st7,j1,j4,j2,j3,j5,j6,za,zb)
      a6vh(5,10)= a6vQLfloop(st7,j1,j4,j2,j3,j5,j6,zb,za)
c-----
      a6vh(5,3) = a6vQLfloop(st7,j4,j1,j2,j3,j6,j5,zb,za)
      a6vh(5,11)= a6vQLfloop(st7,j4,j1,j2,j3,j6,j5,za,zb)
c-----
      a6vh(5,4) = a6vQLfloop(st6,j4,j1,j2,j3,j6,j5,zb,za)
      a6vh(5,12)= a6vQLfloop(st6,j4,j1,j2,j3,j6,j5,za,zb)
c-----
      a6vh(5,5) = a6vQLfloop(st6,j1,j4,j2,j3,j6,j5,za,zb)
      a6vh(5,13)= a6vQLfloop(st6,j1,j4,j2,j3,j6,j5,zb,za)
c-----
      a6vh(5,6) = a6vQLfloop(st7,j1,j4,j2,j3,j6,j5,za,zb)
      a6vh(5,14)= a6vQLfloop(st7,j1,j4,j2,j3,j6,j5,zb,za)
c-----
      a6vh(5,7) = a6vQLfloop(st7,j4,j1,j2,j3,j5,j6,zb,za)
      a6vh(5,15)= a6vQLfloop(st7,j4,j1,j2,j3,j5,j6,za,zb)
c-----
      a6vh(5,8) = a6vQLfloop(st6,j4,j1,j2,j3,j5,j6,zb,za)
      a6vh(5,16)= a6vQLfloop(st6,j4,j1,j2,j3,j5,j6,za,zb)
c-----
c      do i=1,5
c      do j=1,8
c         write(*,*) a6vh(i,j) 
c      enddo
c      enddo
c      stop
c-----donehere  
      return
      end
