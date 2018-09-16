      subroutine LRcalc(k1,k2,k3,k4,k5,k6,za,zb,A,xpp,xmp,xpm,xmm)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: k1,k2,k3,k4,k5,k6,l3,l4,l5,l6,h3,h5
      complex(dp):: A(6),Xpp(2,2),Xmp(2,2),Xpm(2,2),Xmm(2,2),
     & zab2,zba2
      real(dp):: s12,s34,s56
C---statementfunctions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      zba2(k1,k2,k3,k4)=zb(k1,k2)*za(k2,k4)+zb(k1,k3)*za(k3,k4)
C---end statement functions 

      s12=s(k1,k2)
      s34=s(k3,k4)
      s56=s(k5,k6)
      do h3=1,2
         if (h3 ==1) then
             l3=k3
             l4=k4
         elseif (h3 ==2) then
             l3=k4
             l4=k3
         endif
      do h5=1,2
         if (h5 ==1) then
             l5=k5
             l6=k6
         elseif (h5 ==2) then
             l5=k6
             l6=k5
         endif

c      Xmp(h3,h5)=(
c     & A(2)*s12*za(l3,l5)*zb(l4,l6)*zab2(k1,k3,k4,k2)*izab2(k2,k3,k4,k1)
c     & -(A(3)+A(4))*za(k1,l3)*za(k1,l5)*zb(k2,l4)*zb(k2,l6)
c     & -A(5)*za(l3,l5)*zb(k2,l4)*zb(k2,l6)*zab2(k1,k3,k4,k2)*izb(k1,k2)
c     & +(A(5)+A(6))*za(k1,l3)*za(k1,l5)*zb(l4,l6)*zab2(k1,k3,k4,k2)
c     & *iza(k1,k2))/s34/s56/s12
c      write(6,*) 'old,h3,h5,Xmp(h3,h5)',h3,h5,Xmp(h3,h5)

      Xmp(h3,h5)=(
     &  A(2)*s12*za(l3,l5)*zb(l4,l6)*zab2(k1,k3,k4,k2)/zab2(k2,k3,k4,k1)
     & -(A(3)+A(4))*za(k1,l3)*za(k1,l5)*zb(k2,l4)*zb(k2,l6)
     & +A(5)*zab2(k1,k3,k4,k2)
     & *(za(l3,l5)*zb(k2,l4)*zb(l6,k2)/zb(k1,k2)
     &  +za(k1,l3)*za(k1,l5)*zb(l4,l6)/za(k1,k2)))/s34/s56/s12

c      Xpp(h3,h5)=iza(k1,k2)/s12/s34/s56*(
c     & +A(1)*s12*za(l3,l5)*zb(k1,k2)*zb(l4,l6)
c     & +A(3)*s12*za(l3,l5)*zb(k1,l6)*zb(l4,k2)
c     & +A(4)*s12*za(l3,l5)*zb(k1,l4)*zb(k2,l6)
c     & +A(5)*za(l3,l5)*zb(k2,l4)*zb(k2,l6)*zab2(k2,k3,k4,k1)
c     & -(A(5)+A(6))*za(l3,l5)*zb(k1,l4)*zb(k1,l6)*zab2(k1,k3,k4,k2))
c      write(6,*) 'old,h3,h5,Xpp(h3,h5)',h3,h5,Xpp(h3,h5)

      Xpp(h3,h5)=za(l3,l5)/(za(k1,k2)*s12*s34*s56)*(
     & +(A(3)+A(1))*s12*zb(k1,l6)*zb(l4,k2)
     & -(A(4)+A(1))*s12*zb(k2,l6)*zb(l4,k1)
     & +A(5)*(zb(k2,l4)*zb(k2,l6)*zab2(k2,k3,k4,k1)
     &       -zb(k1,l4)*zb(k1,l6)*zab2(k1,k3,k4,k2)))

c--- these are obtained by c.c. of above
c      Xpm(3-h3,3-h5)=(
c     & A(2)*s12*zb(l3,l5)*za(l4,l6)*zab2(k2,k3,k4,k1)*izab2(k1,k3,k4,k2)
c     & -(A(3)+A(4))*zb(k1,l3)*zb(k1,l5)*za(k2,l4)*za(k2,l6)
c     & -A(5)*zb(l3,l5)*za(k2,l4)*za(k2,l6)*zab2(k2,k3,k4,k1)*iza(k1,k2)
c     & +(A(5)+A(6))*zb(k1,l3)*zb(k1,l5)*za(l4,l6)*zab2(k2,k3,k4,k1)
c     & *izb(k1,k2))/s34/s56/s12
c      write(6,*) 'pm:old',Xpm(3-h3,3-h5)

      Xpm(3-h3,3-h5)=(
     &  A(2)*s12*zb(l3,l5)*za(l4,l6)*zba2(k1,k3,k4,k2)/zba2(k2,k3,k4,k1)
     & -(A(3)+A(4))*zb(k1,l3)*zb(k1,l5)*za(k2,l4)*za(k2,l6)
     & +A(5)*zba2(k1,k3,k4,k2)
     & *(zb(l3,l5)*za(k2,l4)*za(l6,k2)/za(k1,k2)
     &  +zb(k1,l3)*zb(k1,l5)*za(l4,l6)/zb(k1,k2)))/s34/s56/s12


c      Xmm(3-h3,3-h5)=izb(k1,k2)/s12/s34/s56*(
c     & +A(1)*s12*zb(l3,l5)*za(k1,k2)*za(l4,l6)
c     & +A(3)*s12*zb(l3,l5)*za(k1,l6)*za(l4,k2)
c     & +A(4)*s12*zb(l3,l5)*za(k1,l4)*za(k2,l6)
c     & +A(5)*zb(l3,l5)*za(k2,l4)*za(k2,l6)*zab2(k1,k3,k4,k2)
c     & -(A(5)+A(6))*zb(l3,l5)*za(k1,l4)*za(k1,l6)*zab2(k2,k3,k4,k1))
c      write(6,*) 'mm:old',Xmm(3-h3,3-h5)

      Xmm(3-h3,3-h5)=zb(l3,l5)/(zb(k1,k2)*s12*s34*s56)*(
     & +(A(3)+A(1))*s12*za(k1,l6)*za(l4,k2)
     & -(A(4)+A(1))*s12*za(k2,l6)*za(l4,k1)
     & +A(5)*(za(k2,l4)*za(k2,l6)*zba2(k2,k3,k4,k1)
     &       -za(k1,l4)*za(k1,l6)*zba2(k1,k3,k4,k2)))
      
      enddo
      enddo

      return
      end


