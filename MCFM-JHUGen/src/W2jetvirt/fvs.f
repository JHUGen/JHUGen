      double complex function Fvs(st,j1,j2,j3,j4,j5,j6,za,zb) 
C--- This function returns the result of BDK
C--- Published in Nucl.Phys.B513:3-86,1998.
C--- e-Print: hep-ph/9708239
C--- Eqs.(11.1) and Eqs.(11.6) as appropriate dependent on choice of "st"
      implicit none
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      character*9 st
      double complex L0,L1,Lsm1_2me,Lnrat
      double complex Brackpm,Brackpp,Brackppa
      double precision t
C---statement functions
      t(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)

C---This is the square bracket in Eq.(11.1) subject to exchange 16,25
      Brackppa(j1,j2,j3,j4,j5,j6)=
     & za(j1,j2)*(za(j5,j3)*zb(j3,j1))**2/za(j5,j6)/za(j3,j4)**2
     & *L1(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)**2
     & +za(j5,j2)*za(j5,j3)*zb(j3,j1)/za(j5,j6)/za(j3,j4)**2
     & *L0(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)

C---This is the whole expression in Eq.(11.1) subject to exchange 34
      Brackpp(j1,j2,j3,j4,j5,j6)=
     & -za(j2,j3)*za(j2,j4)*za(j3,j5)*za(j4,j5)
     & /za(j1,j2)/za(j5,j6)/za(j3,j4)**4
     & *Lsm1_2me(t(j1,j2,j3),t(j1,j2,j4),s(j1,j2),s(j5,j6))
     &  -(za(j2,j4)*za(j3,j5)+za(j2,j3)*za(j4,j5))
     & /za(j1,j2)/za(j5,j6)/za(j3,j4)**3
     & *(za(j5,j3)*zb(j3,j1)*za(j1,j2)
     & *L0(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)
     & +za(j5,j6)*zb(j6,j4)*za(j4,j2)
     & *L0(-t(j1,j2,j3),-s(j5,j6))/s(j5,j6)
     & +0.5d0*za(j5,j2)*(Lnrat(-t(j1,j2,j3),-t(j1,j2,j4))))
     & -Brackppa(j1,j2,j3,j4,j5,j6)-Brackppa(j6,j5,j3,j4,j2,j1)
     & +0.5d0/za(j3,j4)**2
     & *(za(j2,j5)**2/za(j1,j2)/za(j5,j6)
     & -zb(j1,j6)**2/zb(j1,j2)/zb(j5,j6))

C---end statement functions

      if(st.eq.'q+qb-g+g-') then
c---flip2:( 1<-->2, 3<-->4, 5<-->6, za<-->zb)
      Fvs=
     & +Brackpm(j1,j2,j3,j4,j5,j6,za,zb)
     & +Brackpm(j2,j1,j4,j3,j6,j5,zb,za)

      elseif(st.eq.'q+qb-g+g+') then

c---exch34:(3<-->4)
      Fvs=
     & +Brackpp(j1,j2,j3,j4,j5,j6)
     & +Brackpp(j1,j2,j4,j3,j5,j6) 

      endif
      return
      end
      
      double complex function Brackpm(j1,j2,j3,j4,j5,j6,za,zb) 
C---This is the whole expression in Eq.(11.6) subject to exchange flip_2
c---ie  flip2:( 1<-->2, 3<-->4, 5<-->6, za<-->zb)
      implicit none
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex Lsm1_2mh,I3m,Lnrat,zab2,Brackpma
      double precision t,delta,IDelta

C---statement functions
      t(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      delta(j1,j2,j3,j4,j5,j6)=s(j1,j2)-s(j3,j4)-s(j5,j6)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

      IDelta=1d0/(s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2d0*(+s(j1,j2)*s(j3,j4)+s(j1,j2)*s(j5,j6)+s(j5,j6)*s(j3,j4)))

      Brackpm=
     & -2d0*za(j2,j3)*zb(j4,j6)*zab2(j3,j1,j2,j6)
     & *zab2(j2,j5,j6,j4)*t(j1,j2,j3)
     & /za(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)**4
     & *Lsm1_2mh(s(j3,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))

     &+(2d0*za(j2,j3)*zb(j4,j6)*(zb(j1,j2)*za(j5,j6)*za(j2,j3)*zb(j4,j6)

     & -zab2(j5,j2,j4,j1)*zab2(j3,j1,j2,j4))*zab2(j4,j1,j2,j3)
     & *(t(j1,j2,j3)-t(j1,j2,j4))/zab2(j3,j1,j2,j4)**3*IDelta

     & -3d0*(s(j3,j4)*delta(j3,j4,j1,j2,j5,j6)*zab2(j2,j3,j4,j1)
     &  *zab2(j5,j3,j4,j6)*IDelta
     & -za(j2,j3)*zb(j3,j1)*za(j5,j4)*zb(j4,j6)
     & -za(j2,j4)*zb(j4,j1)*za(j5,j3)*zb(j3,j6))

     & *zab2(j4,j1,j2,j3)/zab2(j3,j1,j2,j4)*IDelta
     & -zb(j1,j4)*za(j3,j5)*za(j2,j3)*zb(j4,j6)
     & *zab2(j4,j1,j2,j3)**2/zab2(j3,j1,j2,j4)**2*IDelta
     & -zb(j1,j3)*za(j4,j5)*za(j2,j4)*zb(j3,j6)*IDelta)
     & *I3m(s(j1,j2),s(j3,j4),s(j5,j6))

     & +(2d0*za(j2,j3)*zb(j4,j6)*t(j1,j2,j3)*zab2(j2,j1,j4,j6)
     & /za(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)**3

     & -(zab2(j2,j1,j4,j6)**2
     &  +2d0*za(j2,j3)*zb(j3,j6)*za(j2,j4)*zb(j4,j6))
     & /za(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)**2)
     & *Lnrat(-t(j1,j2,j3),-s(j3,j4))

c---ie  flip3:( 1<-->5, 2<-->6, 3<-->4, za<-->zb)
     &  +Brackpma(j1,j2,j3,j4,j5,j6,za,zb)
     &  +Brackpma(j5,j6,j4,j3,j1,j2,zb,za)

     &  +zb(j1,j6)*zab2(j4,j1,j2,j3)
     &  *(zb(j1,j6)*delta(j3,j4,j5,j6,j1,j2)
     & -2d0*zb(j1,j2)*za(j2,j5)*zb(j5,j6))
     & /zb(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)*IDelta

     &  +(zab2(j2,j1,j4,j6)**2+za(j2,j1)*zb(j1,j6)*za(j2,j5)*zb(j5,j6))
     & /za(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)**2

       return
       end
       
       
      double complex function Brackpma(j1,j2,j3,j4,j5,j6,za,zb) 
C---This is the curly braces expression in Eq.(11.6) subject to flip_3
c---ie  flip3:( 1<-->5, 2<-->6, 3<-->4, za<-->zb)
      implicit none
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex L0,L1,Lnrat,zab2
      double precision t,delta,IDelta

C---statement functions
      t(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      delta(j1,j2,j3,j4,j5,j6)=s(j1,j2)-s(j3,j4)-s(j5,j6)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

      IDelta=1d0/(s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2d0*(s(j1,j2)*s(j3,j4)+s(j1,j2)*s(j5,j6)+s(j5,j6)*s(j3,j4)))

       Brackpma= 
     & 2d0*zb(j1,j2)*za(j2,j3)*zb(j3,j6)
     & /zb(j5,j6)/zab2(j3,j1,j2,j4)**2
     & *(zab2(j3,j1,j2,j6)*zab2(j2,j1,j3,j4)/zab2(j3,j1,j2,j4)
     & *L0(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)
     & -za(j2,j3)*zb(j3,j6)
     & *(L0(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)
     & -1/2d0*L1(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)))

     & -(3d0*delta(j5,j6,j3,j4,j1,j2)*zab2(j2,j3,j4,j1)
     & *zab2(j5,j3,j4,j6)*zab2(j4,j1,j2,j3)
     & /zab2(j3,j1,j2,j4)*IDelta**2

     & +2d0*za(j3,j2)*zb(j2,j1)*zb(j4,j6)*(t(j1,j2,j3)-t(j1,j2,j4))
     & *(za(j2,j5)*t(j1,j2,j3)+za(j2,j1)*zb(j1,j6)*za(j6,j5))
     & /zab2(j3,j1,j2,j4)**3*IDelta

     & -(2d0*za(j2,j3)*zb(j3,j6)*(t(j1,j2,j3)-t(j1,j2,j4))
     & +delta(j1,j2,j3,j4,j5,j6)
     & *(za(j2,j5)*t(j1,j2,j3)
     & +za(j2,j1)*zb(j1,j6)*za(j6,j5))/za(j5,j6))
     & *za(j5,j2)*zb(j2,j1)/zab2(j3,j1,j2,j4)**2*IDelta
 
     & +0.5d0*(za(j1,j2)*zb(j1,j6)**2/zb(j5,j6)
     &  +zb(j1,j2)*za(j2,j5)**2/za(j5,j6)
     & -2d0*za(j2,j5)*zb(j1,j6))*zab2(j4,j1,j2,j3)
     & /zab2(j3,j1,j2,j4)*IDelta)
     & *(Lnrat(-s(j1,j2),-s(j3,j4)))
     
      return
      end
      
      
