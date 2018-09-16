      function BDK1211mp(k1,k2,k3,k4,k5,k6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: BDK1211mp

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: k1,k2,k3,k4,k5,k6,j2,j1
      complex(dp):: L0,L1,Lnrat,i3m,Lsm1_2mh,zab2,zba2
      complex(dp):: T1256,T6521,box1256,box6521,
     & squarebracketex,squarebracketfl3,
     & squarebracketfl3ex,squarebracket
      real(dp):: t,DELTA3,delta
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      zba2(k1,k2,k3,k4)=za(k4,k2)*zb(k2,k1)+za(k4,k3)*zb(k3,k1)
      DELTA3(k1,k2,k3,k4,k5,k6)=
     & +s(k1,k2)**2+s(k3,k4)**2+s(k5,k6)**2-2d0*s(k1,k2)*s(k3,k4)
     & -2d0*s(k3,k4)*s(k5,k6)-2d0*s(k5,k6)*s(k1,k2)
      delta(k1,k2,k3,k4,k5,k6)=s(k1,k2)-s(k3,k4)-s(k5,k6)

      j2=k2
      j1=k1

      T1256=
     & +0.5d0*(3d0*s(k3,k4)*delta(k3,k4,k1,k2,k5,k6)
     & -DELTA3(k1,k2,k3,k4,k5,k6))
     & *(t(k1,k2,k3)-t(k1,k2,k4))*zab2(j2,k3,k4,j1)*zab2(k5,k3,k4,k6)
     & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)**2

     & +0.5d0*s(k3,k4)*(t(k1,k2,k3)-t(k1,k2,k4))*za(j2,k5)*zb(j1,k6)
     & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

     & -zb(j1,k2)*za(k5,k6)*zab2(j2,k3,k4,k6)*zab2(k2,k3,k4,k6)
     & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

     & +za(j2,k3)*zb(k4,k6)*delta(k3,k4,k1,k2,k5,k6)
     & *(za(k5,k2)*zb(k2,j1)*delta(k1,k2,k3,k4,k5,k6)
     & -za(k5,k6)*zb(k6,j1)*delta(k5,k6,k3,k4,k1,k2))
     & /zab2(k3,k1,k2,k4)**2/DELTA3(k1,k2,k3,k4,k5,k6)

     & +za(j2,k3)*zb(k4,k6)*zab2(k5,k2,k3,j1)/zab2(k3,k1,k2,k4)**2

     & -2d0*zb(j1,k3)*za(k4,k5)*zab2(j2,k1,k3,k4)*zab2(k3,k1,k2,k6)
     & /t(k1,k2,k3)**2/zab2(k3,k1,k2,k4)

      T6521=
     & +0.5d0*(3d0*s(k3,k4)*delta(k3,k4,k1,k2,k5,k6)
     & -DELTA3(k1,k2,k3,k4,k5,k6))
     & *(t(k1,k2,k4)-t(k1,k2,k3))*zab2(k5,k3,k4,k6)*zab2(j2,k3,k4,j1)
     & /zab2(k3,k5,k6,k4)/DELTA3(k1,k2,k3,k4,k5,k6)**2

     & +0.5d0*s(k3,k4)*(t(k1,k2,k4)-t(k1,k2,k3))*za(k5,j2)*zb(k6,j1)
     & /zab2(k3,k5,k6,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

     & -zb(k6,k5)*za(j2,k1)*zab2(k5,k3,k4,k1)*zab2(k5,k3,k4,j1)
     & /zab2(k3,k5,k6,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

     & +za(k5,k3)*zb(k4,j1)*delta(k3,k4,k1,k2,k5,k6)
     & *(za(j2,k5)*zb(k5,k6)*delta(k5,k6,k3,k4,k1,k2)
     & -za(j2,k1)*zb(k1,k6)*delta(k1,k2,k3,k4,k5,k6))
     & /zab2(k3,k5,k6,k4)**2/DELTA3(k1,k2,k3,k4,k5,k6)

     & +za(k5,k3)*zb(k4,j1)*zab2(j2,k5,k3,k6)/zab2(k3,k5,k6,k4)**2

     & -2d0*zb(k6,k3)*za(k4,j2)*zab2(k5,k6,k3,k4)*zab2(k3,k6,k5,j1)
     & /t(k6,k5,k3)**2/zab2(k3,k5,k6,k4)

c      write(6,*) 'T1256',T1256
c      write(6,*) 'T6521',T6521
c      write(6,*) 'T1256+T6521',T1256+T6521
c      pause

C     T1256i=
C    &   +0.5d0*(3d0*s(k3,k4)*delta(k3,k4,k1,k2,k5,k6)
C    & -DELTA3(k1,k2,k3,k4,k5,k6))
C    & *(t(k1,k2,k3)-t(k1,k2,k4))*zab2(k2,k3,k4,k1)*zab2(k5,k3,k4,k6)
C    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)**2

C    & +0.5d0*s(k3,k4)*(t(k1,k2,k3)-t(k1,k2,k4))*za(k2,k5)*zb(k1,k6)
C    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

C    & -zb(k1,k2)*zab2(k2,k3,k4,k6)*(
C    & -zab2(k2,k3,k4,k1)*za(k5,k1)-zab2(k2,k3,k4,k2)*za(k5,k2)
C    & -zab2(k2,k3,k4,k3)*za(k5,k3)-zab2(k2,k3,k4,k4)*za(k5,k4))
C    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

C    & +za(k2,k3)*zb(k4,k6)*delta(k3,k4,k1,k2,k5,k6)
C    & *(za(k5,k2)*zb(k2,k1)*delta(k1,k2,k3,k4,k5,k6)
C    & +(+za(k5,k2)*zb(k2,k1)+za(k5,k3)*zb(k3,k1)+za(k5,k4)*zb(k4,k1))
C    & *delta(k5,k6,k3,k4,k1,k2))
C    & /zab2(k3,k1,k2,k4)**2/DELTA3(k1,k2,k3,k4,k5,k6)

C    & +za(k2,k3)*zb(k4,k6)*zab2(k5,k2,k3,k1)/zab2(k3,k1,k2,k4)**2

C    & -2d0*zb(k1,k3)*za(k4,k5)*zab2(k2,k1,k3,k4)*zab2(k3,k1,k2,k6)
C    & /t(k1,k2,k3)**2/zab2(k3,k1,k2,k4)

C     T6521i=
C    & -0.5d0*(3d0*s(k3,k4)*delta(k3,k4,k1,k2,k5,k6)
C    & -DELTA3(k1,k2,k3,k4,k5,k6))
C    & *(t(k1,k2,k4)-t(k1,k2,k3))*zab2(k5,k3,k4,k6)*zab2(k2,k3,k4,k1)
C    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)**2

C    & -0.5d0*s(k3,k4)*(t(k1,k2,k4)-t(k1,k2,k3))*za(k5,k2)*zb(k6,k1)
C    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

C    & -za(k2,k1)*(zb(k6,k1)*zab2(k1,k3,k4,k1)
C    & +zb(k6,k2)*zab2(k2,k3,k4,k1)
C    & +zb(k6,k3)*zab2(k3,k3,k4,k1)+zb(k6,k4)*zab2(k4,k3,k4,k1))
C    & *zab2(k5,k3,k4,k1)
C    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

C    & +za(k5,k3)*zb(k4,k1)*delta(k3,k4,k1,k2,k5,k6)
C    & *(-(za(k2,k1)*zb(k1,k6)+za(k2,k3)*zb(k3,k6)+za(k2,k4)*zb(k4,k6))
C    & *delta(k5,k6,k3,k4,k1,k2)
C    & -za(k2,k1)*zb(k1,k6)*delta(k1,k2,k3,k4,k5,k6))
C    & /zab2(k3,k1,k2,k4)**2/DELTA3(k1,k2,k3,k4,k5,k6)

C    & -za(k5,k3)*zb(k4,k1)*zab2(k2,k1,k4,k6)/zab2(k3,k1,k2,k4)**2

C    & -2d0*zb(k6,k3)*za(k4,k2)*(-zab2(k5,k1,k2,k4))*zab2(k3,k2,k4,k1)
C    & /t(k6,k5,k3)**2/zab2(k3,k1,k2,k4)

c      write(6,*) 'T1256',T1256
c      write(6,*) 'T6521',T6521
      squarebracket=
     & 0.5d0*zb(k6,k4)**2*za(k4,k2)**2
     & /za(k1,k2)/zb(k5,k6)/zab2(k3,k1,k2,k4)
     & *L1(-s(k5,k6),-t(k1,k2,k3))/t(k1,k2,k3)

     & +2d0*zb(k6,k4)*za(k4,k2)*zab2(k2,k1,k3,k6)
     & /za(k1,k2)/zb(k5,k6)/zab2(k3,k1,k2,k4)
     & *L0(-t(k1,k2,k3),-s(k5,k6))/s(k5,k6)

     & -za(k2,k3)*za(k2,k4)*zb(k6,k4)**2*t(k1,k2,k3)
     & /za(k1,k2)/zb(k5,k6)/zab2(k3,k1,k2,k4)**2
     & *L0(-t(k1,k2,k3),-s(k5,k6))/s(k5,k6)

     & -0.5d0*za(k2,k3)*zb(k6,k4)*zab2(k2,k1,k3,k6)
     & /za(k1,k2)/zb(k5,k6)/zab2(k3,k1,k2,k4)**2
     & *(Lnrat(-s(k3,k4),-s(k5,k6))+Lnrat(-t(k1,k2,k3),-s(k5,k6)))

     & -0.75d0*zab2(k2,k1,k3,k6)**2
     & /za(k1,k2)/zb(k5,k6)/t(k1,k2,k3)/zab2(k3,k1,k2,k4)
     & *(Lnrat(-s(k3,k4),-s(k5,k6))+Lnrat(-t(k1,k2,k3),-s(k5,k6)))

     & +(1.5d0*delta(k5,k6,k1,k2,k3,k4)*(t(k1,k2,k3)-t(k1,k2,k4))
     & *zab2(k2,k3,k4,k1)*zab2(k5,k3,k4,k6)
     & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)**2
     & -zb(k1,k2)*za(k2,k3)*zb(k4,k6)/zab2(k3,k1,k2,k4)**2
     & /DELTA3(k1,k2,k3,k4,k5,k6)
     & *(za(k2,k5)*(t(k1,k2,k3)-t(k1,k2,k4))
     & -2d0*za(k2,k1)*zb(k1,k6)*za(k6,k5))
     & +zb(k1,k2)*za(k2,k5)*(za(k2,k3)*zb(k3,k6)-za(k2,k4)*zb(k4,k6))
     & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)
     & +zb(k1,k6)/zb(k5,k6)/zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)
     & *(za(k2,k3)*zb(k3,k6)*t(k1,k2,k3)
     & -za(k2,k4)*zb(k4,k6)*t(k1,k2,k4)
     & +za(k2,k3)*zb(k4,k6)*delta(k3,k4,k1,k2,k5,k6)*t(k1,k2,k3)
     & /zab2(k3,k1,k2,k4))
     & )*Lnrat(-s(k1,k2),-s(k3,k4))

     & -0.25d0*zb(k1,k6)*(t(k1,k2,k3)-t(k1,k2,k4))
     & *(zb(k1,k6)*delta(k3,k4,k1,k2,k5,k6)
     & -2d0*zb(k1,k2)*za(k2,k5)*zb(k5,k6))
     & /zb(k1,k2)/zb(k5,k6)/zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)
     & ;

C     squarebracketi=
C    & 0.5d0*(-za(k1,k5)*zb(k1,k4)-za(k2,k5)*zb(k2,k4)
C    & -za(k3,k5)*zb(k3,k4))*zb(k6,k4)*za(k4,k2)**2
C    & /za(k1,k2)/s(k5,k6)/zab2(k3,k1,k2,k4)
C    & *L1(-s(k5,k6),-t(k1,k2,k3))/t(k1,k2,k3)

C    & +2d0*(-za(k1,k5)*zb(k1,k4)-za(k2,k5)*zb(k2,k4)
C    & -za(k3,k5)*zb(k3,k4))*za(k4,k2)*zab2(k2,k1,k3,k6)
C    & /za(k1,k2)/s(k5,k6)/zab2(k3,k1,k2,k4)
C    & *L0(-t(k1,k2,k3),-s(k5,k6))/s(k5,k6)

C    & -za(k2,k3)*za(k2,k4)*(-za(k1,k5)*zb(k1,k4)-za(k2,k5)*zb(k2,k4)
C    & -za(k3,k5)*zb(k3,k4))
C    & *zb(k6,k4)*t(k1,k2,k3)
C    & /za(k1,k2)/s(k5,k6)/zab2(k3,k1,k2,k4)**2
C    & *L0(-t(k1,k2,k3),-s(k5,k6))/s(k5,k6)

C    & -0.5d0*za(k2,k3)*(-zb(k1,k4)*za(k1,k5)-zb(k2,k4)*za(k2,k5)
C    * -zb(k3,k4)*za(k3,k5))*zab2(k2,k1,k3,k6)
C    & /za(k1,k2)/s(k5,k6)/zab2(k3,k1,k2,k4)**2
C    & *(Lnrat(-s(k3,k4),-s(k5,k6))+Lnrat(-t(k1,k2,k3),-s(k5,k6)))

C    & -0.75d0*zab2(k2,k1,k3,k6)
C    & *(-zab2(k2,k1,k3,k1)*za(k1,k5)-zab2(k2,k1,k3,k2)*za(k2,k5)
C    & -zab2(k2,k1,k3,k3)*za(k3,k5)-zab2(k2,k1,k3,k4)*za(k4,k5))
C    & /za(k1,k2)/s(k5,k6)/t(k1,k2,k3)/zab2(k3,k1,k2,k4)
C    & *(Lnrat(-s(k3,k4),-s(k5,k6))+Lnrat(-t(k1,k2,k3),-s(k5,k6)))

C    & +(1.5d0*delta(k5,k6,k1,k2,k3,k4)*(t(k1,k2,k3)-t(k1,k2,k4))
C    & *zab2(k2,k3,k4,k1)*zab2(k5,k3,k4,k6)
C    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)**2
C    & -zb(k1,k2)*za(k2,k3)*zb(k4,k6)/zab2(k3,k1,k2,k4)**2
C    & /DELTA3(k1,k2,k3,k4,k5,k6)
C    & *(za(k2,k5)*(t(k1,k2,k3)-t(k1,k2,k4))-2d0*za(k2,k1)
C    & *(-zb(k1,k2)*za(k2,k5)-zb(k1,k3)*za(k3,k5)-zb(k1,k4)*za(k4,k5)))
C    & +zb(k1,k2)*za(k2,k5)*(za(k2,k3)*zb(k3,k6)-za(k2,k4)*zb(k4,k6))
C    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)
C    & +(-zb(k1,k2)*za(k2,k5)-zb(k1,k3)*za(k3,k5)-zb(k1,k4)*za(k4,k5))
C    & /s(k5,k6)/zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)
C    & *(za(k2,k3)*zb(k3,k6)*t(k1,k2,k3)
C    & -za(k2,k4)*zb(k4,k6)*t(k1,k2,k4)
C    & +za(k2,k3)*zb(k4,k6)*delta(k3,k4,k1,k2,k5,k6)*t(k1,k2,k3)
C    & /zab2(k3,k1,k2,k4))
C    & )*Lnrat(-s(k1,k2),-s(k3,k4))

C    & -0.25d0*(-zb(k1,k2)*za(k2,k5)-zb(k1,k3)*za(k3,k5)
C    & -zb(k1,k4)*za(k4,k5))
C    & *(t(k1,k2,k3)-t(k1,k2,k4))
C    & *(zb(k1,k6)*delta(k3,k4,k1,k2,k5,k6)-2d0*zb(k1,k2)
C    & *(-za(k2,k1)*zb(k1,k6)-za(k2,k3)*zb(k3,k6)-za(k2,k4)*zb(k4,k6)))
C    & /zb(k1,k2)/s(k5,k6)/zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

      squarebracketex=
     & 0.5d0*zb(k1,k4)**2*za(k4,k5)**2
     & /za(k6,k5)/zb(k2,k1)/zab2(k3,k6,k5,k4)
     & *L1(-s(k2,k1),-t(k6,k5,k3))/t(k6,k5,k3)

     & +2d0*zb(k1,k4)*za(k4,k5)*zab2(k5,k6,k3,k1)
     & /za(k6,k5)/zb(k2,k1)/zab2(k3,k6,k5,k4)
     & *L0(-t(k6,k5,k3),-s(k2,k1))/s(k2,k1)

     & -za(k5,k3)*za(k5,k4)*zb(k1,k4)**2*t(k6,k5,k3)
     & /za(k6,k5)/zb(k2,k1)/zab2(k3,k6,k5,k4)**2
     & *L0(-t(k6,k5,k3),-s(k2,k1))/s(k2,k1)

     & -0.5d0*za(k5,k3)*zb(k1,k4)*zab2(k5,k6,k3,k1)
     & /za(k6,k5)/zb(k2,k1)/zab2(k3,k6,k5,k4)**2
     & *(Lnrat(-s(k3,k4),-s(k2,k1))+Lnrat(-t(k6,k5,k3),-s(k2,k1)))

     & -0.75d0*zab2(k5,k6,k3,k1)**2
     & /za(k6,k5)/zb(k2,k1)/t(k6,k5,k3)/zab2(k3,k6,k5,k4)
     & *(Lnrat(-s(k3,k4),-s(k2,k1))+Lnrat(-t(k6,k5,k3),-s(k2,k1)))

     & +(1.5d0*delta(k2,k1,k6,k5,k3,k4)*(t(k6,k5,k3)-t(k6,k5,k4))
     & *zab2(k5,k3,k4,k6)*zab2(k2,k3,k4,k1)
     & /zab2(k3,k6,k5,k4)/DELTA3(k6,k5,k3,k4,k2,k1)**2
     & -zb(k6,k5)*za(k5,k3)*zb(k4,k1)/zab2(k3,k6,k5,k4)**2
     * /DELTA3(k6,k5,k3,k4,k2,k1)
     & *(za(k5,k2)*(t(k6,k5,k3)-t(k6,k5,k4))
     & -2d0*za(k5,k6)*zb(k6,k1)*za(k1,k2))
     & +zb(k6,k5)*za(k5,k2)*(za(k5,k3)*zb(k3,k1)-za(k5,k4)*zb(k4,k1))
     & /zab2(k3,k6,k5,k4)/DELTA3(k6,k5,k3,k4,k2,k1)
     & +zb(k6,k1)/zb(k2,k1)/zab2(k3,k6,k5,k4)/DELTA3(k6,k5,k3,k4,k2,k1)
     & *(za(k5,k3)*zb(k3,k1)*t(k6,k5,k3)
     & -za(k5,k4)*zb(k4,k1)*t(k6,k5,k4)
     & +za(k5,k3)*zb(k4,k1)*delta(k3,k4,k6,k5,k2,k1)*t(k6,k5,k3)
     & /zab2(k3,k6,k5,k4))
     & )*Lnrat(-s(k6,k5),-s(k3,k4))

     & -0.25d0*zb(k6,k1)*(t(k6,k5,k3)-t(k6,k5,k4))
     & *(zb(k6,k1)*delta(k3,k4,k6,k5,k2,k1)
     & -2d0*zb(k6,k5)*za(k5,k2)*zb(k2,k1))
     & /zb(k6,k5)/zb(k2,k1)/zab2(k3,k6,k5,k4)/DELTA3(k6,k5,k3,k4,k2,k1)

C     squarebracketexi=
C    & 0.5d0*zb(k1,k4)**2*za(k4,k5)*(-za(k4,k1)*zb(k1,k6)
C    & -za(k4,k2)*zb(k2,k6)-za(k4,k3)*zb(k3,k6))
C    & /s(k5,k6)/zb(k2,k1)*(-1d0/zab2(k3,k1,k2,k4))
C    & *L1(-s(k2,k1),-t(k6,k5,k3))/t(k6,k5,k3)

C    & +2d0*zb(k1,k4)*(-za(k4,k1)*zb(k1,k6)-za(k4,k2)*zb(k2,k6)
C    & -za(k4,k3)*zb(k3,k6))*(-zab2(k5,k2,k4,k1))
C    & /s(k5,k6)/zb(k2,k1)*(-1d0/zab2(k3,k1,k2,k4))
C    & *L0(-t(k6,k5,k3),-s(k2,k1))/s(k2,k1)

C    & -(-za(k1,k3)*zb(k1,k6)-za(k2,k3)*zb(k2,k6)-za(k4,k3)*zb(k4,k6))
C    & *za(k5,k4)*zb(k1,k4)**2*t(k6,k5,k3)
C    & /s(k5,k6)/zb(k2,k1)/zab2(k3,k1,k2,k4)**2
C    & *L0(-t(k6,k5,k3),-s(k2,k1))/s(k2,k1)

C    & -0.5d0*(-za(k1,k3)*zb(k1,k6)-za(k2,k3)*zb(k2,k6)
C    & -za(k4,k3)*zb(k4,k6))*zb(k1,k4)
C    & *(-zab2(k5,k2,k4,k1))
C    & /s(k5,k6)/zb(k2,k1)/zab2(k3,k1,k2,k4)**2
C    & *(Lnrat(-s(k3,k4),-s(k2,k1))+Lnrat(-t(k6,k5,k3),-s(k2,k1)))

C    & -0.75d0*zab2(k5,k2,k4,k1)
C    & *(-zab2(k1,k2,k4,k1)*zb(k1,k6)-zab2(k2,k2,k4,k1)*zb(k2,k6)
C    & -zab2(k3,k2,k4,k1)*zb(k3,k6)-zab2(k4,k2,k4,k1)*zb(k4,k6))
C    & /s(k5,k6)/zb(k2,k1)/t(k6,k5,k3)*(-1d0/zab2(k3,k1,k2,k4))
C    & *(Lnrat(-s(k3,k4),-s(k2,k1))+Lnrat(-t(k6,k5,k3),-s(k2,k1)))

C    & +(1.5d0*delta(k2,k1,k6,k5,k3,k4)*(t(k6,k5,k3)-t(k6,k5,k4))
C    & *zab2(k5,k3,k4,k6)*zab2(k2,k3,k4,k1)
C    & *(-1d0/zab2(k3,k1,k2,k4))/DELTA3(k6,k5,k3,k4,k2,k1)**2
C    & -(-zb(k6,k1)*za(k1,k3)-zb(k6,k2)*za(k2,k3)-zb(k6,k4)*za(k4,k3))
C    & *zb(k4,k1)/zab2(k3,k1,k2,k4)**2/DELTA3(k6,k5,k3,k4,k2,k1)
C    & *(za(k5,k2)*(t(k6,k5,k3)-t(k6,k5,k4))
C    & -2d0*(-za(k5,k2)*zb(k2,k1)-za(k5,k3)*zb(k3,k1)
C    & -za(k5,k4)*zb(k4,k1))*za(k1,k2))
C    & +(-zb(k6,k1)*za(k1,k2)-zb(k6,k3)*za(k3,k2)-zb(k6,k4)*za(k4,k2))
C    & *(za(k5,k3)*zb(k3,k1)-za(k5,k4)*zb(k4,k1))
C    & *(-1d0/zab2(k3,k1,k2,k4))/DELTA3(k6,k5,k3,k4,k2,k1)
C    & +zb(k6,k1)/zb(k2,k1)*(-1d0/zab2(k3,k1,k2,k4))
C    & /DELTA3(k6,k5,k3,k4,k2,k1)
C    & *(za(k5,k3)*zb(k3,k1)*t(k6,k5,k3)
C    & -za(k5,k4)*zb(k4,k1)*t(k6,k5,k4)
C    & +za(k5,k3)*zb(k4,k1)*delta(k3,k4,k6,k5,k2,k1)*t(k6,k5,k3)
C    & *(-1d0/zab2(k3,k1,k2,k4)))
C    & )*Lnrat(-s(k6,k5),-s(k3,k4))

C    & -0.25d0*(-za(k5,k2)*zb(k2,k1)-za(k5,k3)*zb(k3,k1)
C    & -za(k5,k4)*zb(k4,k1))
C    & *(t(k6,k5,k3)-t(k6,k5,k4))
C    & *(zb(k6,k1)*delta(k3,k4,k6,k5,k2,k1)-2
C    & *(-zb(k6,k1)*za(k1,k2)-zb(k6,k3)*za(k3,k2)-zb(k6,k4)*za(k4,k2))
C    & *zb(k2,k1))
C    & /s(k6,k5)/zb(k2,k1)*(-1d0/zab2(k3,k1,k2,k4))
C    & /DELTA3(k6,k5,k3,k4,k2,k1)

C     squarebracketfl3i=
C    & 0.5d0*za(k2,k3)**2*zb(k3,k6)*(-zb(k3,k1)*za(k1,k5)
C    & -zb(k3,k2)*za(k2,k5)-zb(k3,k4)*za(k4,k5))
C    & /s(k5,k6)/za(k1,k2)*(-1d0/zba2(k4,k1,k2,k3))
C    & *L1(-s(k1,k2),-t(k5,k6,k4))/t(k5,k6,k4)

C    & +2d0*za(k2,k3)*(-zb(k3,k1)*za(k1,k5)-zb(k3,k2)*za(k2,k5)
C    & -zb(k3,k4)*za(k4,k5))*(-zba2(k6,k1,k3,k2))
C    & /s(k5,k6)/za(k1,k2)*(-1d0/zba2(k4,k1,k2,k3))
C    & *L0(-t(k5,k6,k4),-s(k1,k2))/s(k1,k2)

C    & -(-za(k1,k5)*zb(k1,k4)-za(k2,k5)*zb(k2,k4)-za(k3,k5)*zb(k3,k4))
C    & *zb(k6,k3)*za(k2,k3)**2*t(k5,k6,k4)
C    & /s(k5,k6)/za(k1,k2)*(-1d0/zba2(k4,k1,k2,k3))**2
C    & *L0(-t(k5,k6,k4),-s(k1,k2))/s(k1,k2)

C    & -0.5d0*(-za(k1,k5)*zb(k1,k4)-za(k2,k5)*zb(k2,k4)
C    & -za(k3,k5)*zb(k3,k4))*za(k2,k3)*(-zba2(k6,k1,k3,k2))
C    & /s(k5,k6)/za(k1,k2)*(-1d0/zba2(k4,k1,k2,k3))**2
C    & *(Lnrat(-s(k4,k3),-s(k1,k2))+Lnrat(-t(k5,k6,k4),-s(k1,k2)))

C    & -0.75d0*zba2(k6,k1,k3,k2)*(-zba2(k1,k1,k3,k2)*za(k1,k5)
C    & -zba2(k2,k1,k3,k2)*za(k2,k5)
C    & -zba2(k3,k1,k3,k2)*za(k3,k5)-zba2(k4,k1,k3,k2)*za(k4,k5))
C    & /s(k5,k6)/za(k1,k2)/t(k5,k6,k4)*(-1d0/zba2(k4,k1,k2,k3))
C    & *(Lnrat(-s(k4,k3),-s(k1,k2))+Lnrat(-t(k5,k6,k4),-s(k1,k2)))

C    & +(1.5d0*delta(k1,k2,k5,k6,k4,k3)*(t(k5,k6,k4)-t(k5,k6,k3))
C    & *zba2(k6,k4,k3,k5)*zba2(k1,k4,k3,k2)
C    & *(-1d0/zba2(k4,k1,k2,k3))/DELTA3(k5,k6,k4,k3,k1,k2)**2
C    & -(-za(k5,k1)*zb(k1,k4)-za(k5,k2)*zb(k2,k4)-za(k5,k3)*zb(k3,k4))
C    & *za(k3,k2)
C    & *(-1d0/zba2(k4,k1,k2,k3))**2/DELTA3(k5,k6,k4,k3,k1,k2)
C    & *(zb(k6,k1)*(t(k5,k6,k4)-t(k5,k6,k3))-2d0*(-zb(k6,k1)*za(k1,k2)
C    & -zb(k6,k3)*za(k3,k2)-zb(k6,k4)*za(k4,k2))*zb(k2,k1))
C    & +(-za(k5,k2)*zb(k2,k1)-za(k5,k3)*zb(k3,k1)-za(k5,k4)*zb(k4,k1))
C    & *(zb(k6,k4)*za(k4,k2)-zb(k6,k3)*za(k3,k2))
C    & *(-1d0/zba2(k4,k1,k2,k3))/DELTA3(k5,k6,k4,k3,k1,k2)
C    & +za(k5,k2)/za(k1,k2)*(-1d0/zba2(k4,k1,k2,k3))
C    & /DELTA3(k5,k6,k4,k3,k1,k2)
C    & *(zb(k6,k4)*za(k4,k2)*t(k5,k6,k4)
C    & -zb(k6,k3)*za(k3,k2)*t(k5,k6,k3)
C    & +zb(k6,k4)*za(k3,k2)*delta(k4,k3,k5,k6,k1,k2)*t(k5,k6,k4)
C    & *(-1d0/zba2(k4,k1,k2,k3)))
C    & )*Lnrat(-s(k5,k6),-s(k4,k3))

C    & -0.25d0*(-zb(k6,k1)*za(k1,k2)-zb(k6,k3)*za(k3,k2)
C    & -zb(k6,k4)*za(k4,k2))*(t(k5,k6,k4)-t(k5,k6,k3))
C    & *(za(k5,k2)*delta(k4,k3,k5,k6,k1,k2)-2d0*(-za(k5,k2)*zb(k2,k1)
C    & -za(k5,k3)*zb(k3,k1)-za(k5,k4)*zb(k4,k1))
C    & *za(k1,k2))
C    & /s(k5,k6)/za(k1,k2)
C    & *(-1d0/zba2(k4,k1,k2,k3))/DELTA3(k5,k6,k4,k3,k1,k2)

      squarebracketfl3=
     & 0.5d0*za(k2,k3)**2*zb(k3,k6)**2
     & /zb(k5,k6)/za(k1,k2)/zba2(k4,k5,k6,k3)
     & *L1(-s(k1,k2),-t(k5,k6,k4))/t(k5,k6,k4)

     & +2d0*za(k2,k3)*zb(k3,k6)*zba2(k6,k5,k4,k2)
     & /zb(k5,k6)/za(k1,k2)/zba2(k4,k5,k6,k3)
     & *L0(-t(k5,k6,k4),-s(k1,k2))/s(k1,k2)

     & -zb(k6,k4)*zb(k6,k3)*za(k2,k3)**2*t(k5,k6,k4)
     & /zb(k5,k6)/za(k1,k2)/zba2(k4,k5,k6,k3)**2
     & *L0(-t(k5,k6,k4),-s(k1,k2))/s(k1,k2)

     & -0.5d0*zb(k6,k4)*za(k2,k3)*zba2(k6,k5,k4,k2)
     & /zb(k5,k6)/za(k1,k2)/zba2(k4,k5,k6,k3)**2
     & *(Lnrat(-s(k4,k3),-s(k1,k2))+Lnrat(-t(k5,k6,k4),-s(k1,k2)))

     & -0.75d0*zba2(k6,k5,k4,k2)**2
     & /zb(k5,k6)/za(k1,k2)/t(k5,k6,k4)/zba2(k4,k5,k6,k3)
     & *(Lnrat(-s(k4,k3),-s(k1,k2))+Lnrat(-t(k5,k6,k4),-s(k1,k2)))

     & +(1.5d0*delta(k1,k2,k5,k6,k4,k3)*(t(k5,k6,k4)-t(k5,k6,k3))
     & *zba2(k6,k4,k3,k5)*zba2(k1,k4,k3,k2)
     & /zba2(k4,k5,k6,k3)/DELTA3(k5,k6,k4,k3,k1,k2)**2
     & -za(k5,k6)*zb(k6,k4)*za(k3,k2)/zba2(k4,k5,k6,k3)**2
     & /DELTA3(k5,k6,k4,k3,k1,k2)
     & *(zb(k6,k1)*(t(k5,k6,k4)-t(k5,k6,k3))
     & -2d0*zb(k6,k5)*za(k5,k2)*zb(k2,k1))
     & +za(k5,k6)*zb(k6,k1)*(zb(k6,k4)*za(k4,k2)-zb(k6,k3)*za(k3,k2))
     & /zba2(k4,k5,k6,k3)/DELTA3(k5,k6,k4,k3,k1,k2)
     & +za(k5,k2)/za(k1,k2)/zba2(k4,k5,k6,k3)/DELTA3(k5,k6,k4,k3,k1,k2)
     & *(zb(k6,k4)*za(k4,k2)*t(k5,k6,k4)
     & -zb(k6,k3)*za(k3,k2)*t(k5,k6,k3)
     & +zb(k6,k4)*za(k3,k2)*delta(k4,k3,k5,k6,k1,k2)*t(k5,k6,k4)
     & /zba2(k4,k5,k6,k3))
     & )*Lnrat(-s(k5,k6),-s(k4,k3))

     & -0.25d0*za(k5,k2)*(t(k5,k6,k4)-t(k5,k6,k3))
     & *(za(k5,k2)*delta(k4,k3,k5,k6,k1,k2)
     & -2d0*za(k5,k6)*zb(k6,k1)*za(k1,k2))
     & /za(k5,k6)/za(k1,k2)/zba2(k4,k5,k6,k3)/DELTA3(k5,k6,k4,k3,k1,k2)


      squarebracketfl3ex=
     & 0.5d0*za(k5,k3)**2*zb(k3,k1)**2/zb(k2,k1)/za(k6,k5)
     & /zba2(k4,k2,k1,k3)
     & *L1(-s(k6,k5),-t(k2,k1,k4))/t(k2,k1,k4)

     & +2d0*za(k5,k3)*zb(k3,k1)*zba2(k1,k2,k4,k5)
     & /zb(k2,k1)/za(k6,k5)/zba2(k4,k2,k1,k3)
     & *L0(-t(k2,k1,k4),-s(k6,k5))/s(k6,k5)

     & -zb(k1,k4)*zb(k1,k3)*za(k5,k3)**2*t(k2,k1,k4)
     & /zb(k2,k1)/za(k6,k5)/zba2(k4,k2,k1,k3)**2
     & *L0(-t(k2,k1,k4),-s(k6,k5))/s(k6,k5)

     & -0.5d0*zb(k1,k4)*za(k5,k3)*zba2(k1,k2,k4,k5)
     & /zb(k2,k1)/za(k6,k5)/zba2(k4,k2,k1,k3)**2
     & *(Lnrat(-s(k4,k3),-s(k6,k5))+Lnrat(-t(k2,k1,k4),-s(k6,k5)))

     & -0.75d0*zba2(k1,k2,k4,k5)**2
     & /zb(k2,k1)/za(k6,k5)/t(k2,k1,k4)/zba2(k4,k2,k1,k3)
     & *(Lnrat(-s(k4,k3),-s(k6,k5))+Lnrat(-t(k2,k1,k4),-s(k6,k5)))

     & +(1.5d0*delta(k6,k5,k2,k1,k4,k3)*(t(k2,k1,k4)-t(k2,k1,k3))
     & *zba2(k1,k4,k3,k2)*zba2(k6,k4,k3,k5)
     & /zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)**2
     & -za(k2,k1)*zb(k1,k4)*za(k3,k5)/zba2(k4,k2,k1,k3)**2
     & /DELTA3(k2,k1,k4,k3,k6,k5)
     & *(zb(k1,k6)*(t(k2,k1,k4)-t(k2,k1,k3))
     & -2d0*zb(k1,k2)*za(k2,k5)*zb(k5,k6))
     & +za(k2,k1)*zb(k1,k6)*(zb(k1,k4)*za(k4,k5)-zb(k1,k3)*za(k3,k5))
     & /zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)
     & +za(k2,k5)/za(k6,k5)/zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)
     & *(zb(k1,k4)*za(k4,k5)*t(k2,k1,k4)
     & -zb(k1,k3)*za(k3,k5)*t(k2,k1,k3)
     & +zb(k1,k4)*za(k3,k5)*delta(k4,k3,k2,k1,k6,k5)*t(k2,k1,k4)
     & /zba2(k4,k2,k1,k3))
     & )*Lnrat(-s(k2,k1),-s(k4,k3))

     & -0.25d0*za(k2,k5)*(t(k2,k1,k4)-t(k2,k1,k3))
     & *(za(k2,k5)*delta(k4,k3,k2,k1,k6,k5)
     * -2d0*za(k2,k1)*zb(k1,k6)*za(k6,k5))
     & /za(k2,k1)/za(k6,k5)/zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)


c      squarebracketfl3exi=
C    & 0.5d0*za(k5,k3)*(-za(k1,k3)*zb(k1,k6)-za(k2,k3)*zb(k2,k6)
C    & -za(k4,k3)*zb(k4,k6))*zb(k3,k1)**2
C    & /zb(k2,k1)/s(k5,k6)/zba2(k4,k2,k1,k3)
C    & *L1(-s(k6,k5),-t(k2,k1,k4))/t(k2,k1,k4)

C    & +2d0*(-za(k1,k3)*zb(k1,k6)-za(k2,k3)*zb(k2,k6)
C    & -za(k4,k3)*zb(k4,k6))*zb(k3,k1)*zba2(k1,k2,k4,k5)
C    & /zb(k2,k1)/s(k5,k6)/zba2(k4,k2,k1,k3)
C    & *L0(-t(k2,k1,k4),-s(k6,k5))/s(k6,k5)

C    & -zb(k1,k4)*zb(k1,k3)*za(k5,k3)*(-za(k1,k3)*zb(k1,k6)
C    & -za(k2,k3)*zb(k2,k6)-za(k4,k3)*zb(k4,k6))*t(k2,k1,k4)
C    & /zb(k2,k1)/s(k5,k6)/zba2(k4,k2,k1,k3)**2
C    & *L0(-t(k2,k1,k4),-s(k6,k5))/s(k6,k5)

C    & -0.5d0*zb(k1,k4)*(-za(k1,k3)*zb(k1,k6)-za(k2,k3)*zb(k2,k6)
C    & -za(k4,k3)*zb(k4,k6))*zba2(k1,k2,k4,k5)
C    & /zb(k2,k1)/s(k5,k6)/zba2(k4,k2,k1,k3)**2
C    & *(Lnrat(-s(k4,k3),-s(k6,k5))+Lnrat(-t(k2,k1,k4),-s(k6,k5)))

C    & -0.75d0*zba2(k1,k2,k4,k5)*(-zba2(k1,k2,k4,k1)*zb(k1,k6)
C    & -zba2(k1,k2,k4,k2)*zb(k2,k6)
C    & -zba2(k1,k2,k4,k3)*zb(k3,k6)-zba2(k1,k2,k4,k4)*zb(k4,k6))
C    & /zb(k2,k1)/s(k5,k6)/t(k2,k1,k4)/zba2(k4,k2,k1,k3)
C    & *(Lnrat(-s(k4,k3),-s(k6,k5))+Lnrat(-t(k2,k1,k4),-s(k6,k5)))

C    & +(1.5d0*delta(k6,k5,k2,k1,k4,k3)*(t(k2,k1,k4)-t(k2,k1,k3))
C    & *zba2(k1,k4,k3,k2)*zba2(k6,k4,k3,k5)
C    & /zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)**2
C    & -za(k2,k1)*zb(k1,k4)*za(k3,k5)/zba2(k4,k2,k1,k3)**2
C    & /DELTA3(k2,k1,k4,k3,k6,k5)
C    & *(zb(k1,k6)*(t(k2,k1,k4)-t(k2,k1,k3))-2d0*zb(k1,k2)
C    & *(-za(k2,k1)*zb(k1,k6)-za(k2,k3)*zb(k3,k6)-za(k2,k4)*zb(k4,k6)))
C    & +za(k2,k1)*zb(k1,k6)*(zb(k1,k4)*za(k4,k5)-zb(k1,k3)*za(k3,k5))
C    & /zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)
C    & +(-za(k2,k1)*zb(k1,k6)-za(k2,k3)*zb(k3,k6)-za(k2,k4)*zb(k4,k6))
C    * /s(k5,k6)/zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)
C    & *(zb(k1,k4)*za(k4,k5)*t(k2,k1,k4)
C    & -zb(k1,k3)*za(k3,k5)*t(k2,k1,k3)
C    & +zb(k1,k4)*za(k3,k5)*delta(k4,k3,k2,k1,k6,k5)*t(k2,k1,k4)
C    & /zba2(k4,k2,k1,k3))
C    & )*Lnrat(-s(k2,k1),-s(k4,k3))

C    & -0.25d0*(-za(k2,k1)*zb(k1,k6)-za(k2,k3)*zb(k3,k6)
C    & -za(k2,k4)*zb(k4,k6))*(t(k2,k1,k4)-t(k2,k1,k3))
C    & *(za(k2,k5)*delta(k4,k3,k2,k1,k6,k5)-2d0*za(k2,k1)
C    & *(-zb(k1,k2)*za(k2,k5)-zb(k1,k3)*za(k3,k5)-zb(k1,k4)*za(k4,k5)))
C    & /za(k2,k1)/s(k5,k6)/zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)


      box1256=(za(k4,k5)*zb(k1,k3))**2
     & /(zb(k1,k2)*za(k5,k6)*t(k1,k2,k3)*zab2(k4,k1,k2,k3))
     & -(zab2(k3,k1,k2,k6)*zab2(k2,k1,k3,k4))**2
     & /(za(k1,k2)*zb(k5,k6)*t(k1,k2,k3)*zab2(k3,k1,k2,k4)**3)

      box6521=(za(k4,k2)*zb(k6,k3))**2
     & /(zb(k6,k5)*za(k2,k1)*t(k6,k5,k3)*zab2(k4,k6,k5,k3))
     & -(zab2(k3,k6,k5,k1)*zab2(k5,k6,k3,k4))**2
     & /(za(k6,k5)*zb(k2,k1)*t(k6,k5,k3)*zab2(k3,k6,k5,k4)**3)

      BDK1211mp=
     & +Lsm1_2mh(s(k3,k4),t(k1,k2,k3),s(k1,k2),s(k5,k6))*box1256
     & +Lsm1_2mh(s(k3,k4),t(k6,k5,k3),s(k6,k5),s(k2,k1))*box6521
     & +i3m(s(k1,k2),s(k3,k4),s(k5,k6))*(T1256+T6521)
     & +squarebracket+squarebracketex
     & -squarebracketfl3-squarebracketfl3ex
      return
      end
