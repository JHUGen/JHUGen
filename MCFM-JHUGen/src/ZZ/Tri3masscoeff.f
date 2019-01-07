      double complex function Tri3masscoeff(j1,j2,j3,j4,j5,j6,za,zb) 
C-----Tri3massmass calculates the coefficient of the 3 mass triangle,
C-----in the massless limit using the formula of BDK (11.6) and (11.7)

      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6
      double complex Fvfbit1,Fvsbit1,Fvf,Fvs

c---flip2:( j4<-->j3, j2<-->j1, j5<-->j6, za<-->zb)
      Fvf=
     & +Fvfbit1(j1,j2,j3,j4,j5,j6,za,zb) 
     & +Fvfbit1(j2,j1,j4,j3,j6,j5,zb,za) 
      
      Fvs=
     & +Fvsbit1(j1,j2,j3,j4,j5,j6,za,zb) 
     & +Fvsbit1(j2,j1,j4,j3,j6,j5,zb,za) 


      Tri3masscoeff=(Fvf+Fvs)


      return
      end 

      double complex function Fvfbit1(j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4,j5,j6
      double complex Fvfbit0,zab2
      double precision t,ss,tt,m1sq,m2sq,Lsm1_2mhbit
C---This performs the flip3 on the explicit I3m piece.
C---and adds in the piece from the Lsm1_2mh
c---flip3:( j4<-->j5, j3<-->j6, j2<-->j1, za<-->zb)

C---statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      t(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      Lsm1_2mhbit(ss,tt,m1sq,m2sq)=0.5d0*(ss-m1sq-m2sq)+m1sq*m2sq/tt
C---end statement functions

      Fvfbit1=
     & +Fvfbit0(j1,j2,j3,j4,j5,j6,za,zb) 
     & +Fvfbit0(j2,j1,j6,j5,j4,j3,zb,za) 
     & +zab2(j3,j4,j2,j6)**2/(za(j4,j3)*zb(j5,j6)*zab2(j2,j3,j4,j1)**2)
     & *Lsm1_2mhbit(s(j1,j2),t(j2,j3,j4),s(j3,j4),s(j5,j6))
      return
      end
      


      double complex function Fvsbit1(j1,j2,j3,j4,j5,j6,za,zb) 
C---This is the whole expression in Eq.(11.6) subject to exchange flip_2
c---ie  flip2:( j1<-->j2, j3<-->j4, j5<-->j6, za<-->zb)
      implicit none
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex zab2
      double precision t,delta,IDelta,ss,tt,m1sq,m2sq,Lsm1_2mhbit

C---statement functions
      t(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      delta(j1,j2,j3,j4,j5,j6)=s(j1,j2)-s(j3,j4)-s(j5,j6)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      Lsm1_2mhbit(ss,tt,m1sq,m2sq)=0.5d0*(ss-m1sq-m2sq)+m1sq*m2sq/tt
C---end statement functions

      IDelta=1d0/(s(j3,j4)**2+s(j1,j2)**2+s(j5,j6)**2
     & -2d0*(+s(j3,j4)*s(j1,j2)+s(j3,j4)*s(j5,j6)+s(j5,j6)*s(j1,j2)))

      Fvsbit1=
     & +2d0*za(j3,j2)*zb(j1,j6)*zab2(j2,j4,j3,j6)
     & *zab2(j3,j5,j6,j1)*t(j2,j3,j4)
     & /(za(j4,j3)*zb(j5,j6)*zab2(j2,j3,j4,j1)**4)
     & *Lsm1_2mhbit(s(j1,j2),t(j2,j3,j4),s(j3,j4),s(j5,j6))

     &-(2d0*za(j3,j2)*zb(j1,j6)*(zb(j4,j3)*za(j5,j6)*za(j3,j2)*zb(j1,j6)
     & -zab2(j5,j3,j1,j4)*zab2(j2,j3,j4,j1))*zab2(j1,j3,j4,j2)
     & *(t(j2,j3,j4)-t(j1,j3,j4))/zab2(j2,j3,j4,j1)**3*IDelta

     & -3d0*(s(j1,j2)*delta(j2,j1,j4,j3,j5,j6)*zab2(j3,j2,j1,j4)
     &  *zab2(j5,j2,j1,j6)*IDelta
     & -za(j3,j2)*zb(j2,j4)*za(j5,j1)*zb(j1,j6)
     & -za(j3,j1)*zb(j1,j4)*za(j5,j2)*zb(j2,j6))
     & *zab2(j1,j3,j4,j2)/zab2(j2,j3,j4,j1)*IDelta

     & -zb(j4,j1)*za(j2,j5)*za(j3,j2)*zb(j1,j6)
     & *zab2(j1,j3,j4,j2)**2/zab2(j2,j3,j4,j1)**2*IDelta
     & -zb(j4,j2)*za(j1,j5)*za(j3,j1)*zb(j2,j6)*IDelta)

       return
       end
       
       
      double complex function Fvfbit0(j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4,j5,j6
      double complex zab2
C---statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
C---end statement functions
c
      Fvfbit0=0.5d0*zb(j4,j2)*za(j1,j5)
     & *zab2(j3,j2,j4,j6)/zab2(j2,j3,j4,j1)
     & /(s(j3,j4)+s(j3,j2)+s(j4,j2))
      return
      end
      


