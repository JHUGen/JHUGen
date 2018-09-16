      function msqWaa(p1,p2,p3,p4,p5,p6,za,zb)
      implicit none
      include 'types.f'
      real(dp):: msqWaa
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      include 'ckm.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),
     & s34,s345,s346,s3456,s156,s256,s456,Qup,Qdo,s15,s16
      integer:: p1,p2,p3,p4,p5,p6,b5p,b5m,b6p,b6m
      complex(dp):: lo(2,2),iden,iza,izb,zab2,zab3,mwsq
      
c     statement function
      iden(s34)=1._dp/(s34-mwsq)
      iza(p1,p2)=1._dp/za(p1,p2)
      izb(p1,p2)=1._dp/zb(p1,p2)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zab3(p1,p2,p3,p4,p5)=
     & za(p1,p2)*zb(p2,p5)+za(p1,p3)*zb(p3,p5)+za(p1,p4)*zb(p4,p5)
c     end statement function

c--- implementation of complex mass scheme
      mwsq=wmass**2-im*wmass*wwidth

      s34=s(p3,p4)
      s16=s(p1,p6)
      s15=s(p1,p5)
      s345=s(p3,p4)+s(p3,p5)+s(p4,p5)
      s346=s(p3,p4)+s(p3,p6)+s(p4,p6)
      s3456=s(p3,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)+s(p5,p6)
      s156=s(p1,p5)+s(p1,p6)+s(p5,p6)
      s256=s(p2,p5)+s(p2,p6)+s(p5,p6)
      s456=s(p4,p5)+s(p4,p6)+s(p5,p6)

    
      lo(1,1)=zb(p1,p4)**2*(
     & +iden(s34)*iden(s345)*iden(s3456)*za(p3,p4)*za(p2,p6)
     & *zab2(p5,p3,p4,p1)
     & *izb(p1,p5)*izb(p1,p6)

     & +Qd*iden(s34)*iden(s345)*zab2(p5,p3,p4,p1)*za(p3,p4)
     & *izb(p1,p5)*izb(p1,p6)*izb(p2,p6)

     & +iden(s34)*iden(s346)*iden(s3456)*za(p3,p4)*za(p2,p5)
     & *zab2(p6,p3,p4,p1)
     & *izb(p1,p5)*izb(p1,p6)

     & +Qd*iden(s34)*iden(s346)*zab2(p6,p3,p4,p1)*za(p3,p4)
     & *izb(p1,p5)*izb(p1,p6)*izb(p2,p5)

     & +iden(s345)*iden(s3456)*za(p2,p6)*zab2(p3,p4,p5,p1)
     & *izb(p1,p5)*izb(p1,p6)*izb(p4,p5)

     & +iden(s346)*iden(s3456)*za(p2,p5)*zab2(p3,p4,p6,p1)
     & *izb(p1,p5)*izb(p1,p6)*izb(p4,p6)

     & +Qd*iden(s345)*zab2(p3,p4,p5,p1)
     & *izb(p1,p5)*izb(p1,p6)*izb(p2,p6)*izb(p4,p5)

     & +Qd*iden(s346)*zab2(p3,p4,p6,p1)
     & *izb(p1,p5)*izb(p1,p6)*izb(p2,p5)*izb(p4,p6)

     & +iden(s3456)*za(p2,p3)
     & *(zb(p4,p5)*zab2(p5,p4,p6,p1)+zb(p4,p6)*zab2(p6,p4,p5,p1))
     & /s456*izb(p1,p5)*izb(p1,p6)*izb(p4,p5)*izb(p4,p6)

     & -Qd**2*iden(s34)*za(p3,p4)*zb(p1,p2)
     & *izb(p1,p5)*izb(p1,p6)*izb(p2,p5)*izb(p2,p6))
     
      lo(2,2)=
     &  za(p2,p3)**2*(iden(s34)*zab2(p2,p3,p4,p6)*zb(p3,p4)
     & /(za(p1,p5)*za(p2,p5)*za(p2,p6))

     & *(iden(s346)*iden(s3456)*s15+Qu*iden(s346)+Qu**2/s156)
     & +iden(s34)*zab2(p2,p3,p4,p5)*zb(p3,p4)
     & /(za(p1,p6)*za(p2,p6)*za(p2,p5))

     & *(iden(s345)*iden(s3456)*s16+Qu**2/s156+Qu*iden(s345))

     & +(iden(s345)*iden(s3456)*s16+iden(s345)*Qu)
     & /(za(p1,p6)*za(p2,p6)*za(p2,p5)*za(p4,p5))*zab2(p2,p4,p5,p3)

     & +(iden(s346)*iden(s3456)*s15+iden(s346)*Qu)
     & /(za(p1,p5)*za(p2,p5)*za(p2,p6)*za(p4,p6))*zab2(p2,p4,p6,p3)

     & +iden(s3456)*zb(p1,p3)*za(p2,p4)
     & /(za(p4,p5)*za(p4,p6)*za(p2,p5)*za(p2,p6)))


      lo(1,2)=
     & -iden(s34)*iden(s345)*iden(s3456)*zab2(p2,p3,p5,p4)*zb(p1,p6)
     & *(za(p3,p5)*za(p2,p6)*zb(p1,p6)+za(p2,p5)*za(p3,p4)*zb(p1,p4)
     & )*iza(p2,p6)*izb(p1,p5)

     & + Qu*iden(s34)*iden(s345)*za(p4,p6)*(
     & - za(p1,p2)*za(p2,p3)*zab2(p5,p3,p4,p1)*zb(p1,p4)
     & + zab2(p2,p3,p4,p1)*za(p2,p6)*za(p3,p5)*zb(p4,p6)
     & + zab2(p2,p1,p6,p4)*za(p2,p5)*za(p3,p5)*zb(p1,p5)
     & + za(p2,p6)*za(p2,p5)*za(p3,p4)*zb(p1,p4)*zb(p4,p6)
     & )*iza(p1,p6)*iza(p2,p6)*izb(p1,p5)*iza(p4,p6)

     & + iden(s346)*iden(s3456)*za(p2,p5)*(
     & +za(p2,p5)*zab2(p3,p4,p6,p1)*zb(p1,p5)*zb(p4,p5)
     & -za(p2,p3)*zab2(p3,p4,p6,p5)*zb(p1,p3)*zb(p1,p4)
     & -za(p2,p3)*za(p3,p6)*zb(p1,p3)*zb(p1,p5)*zb(p4,p6)
     & )*iza(p2,p6)*izb(p1,p5)*iza(p4,p6)*izb(p4,p5)

     & +iden(s34)*iden(s346)*iden(s3456)*zab2(p3,p2,p5,p1)
     & *(zab2(p2,p4,p6,p1)*zb(p4,p6)-za(p2,p3)*zb(p1,p4)*zb(p3,p6))
     & *iza(p2,p6)*izb(p1,p5)*za(p2,p5)

     & + Qd*iden(s34)*iden(s346)*zab2(p3,p2,p5,p1)
     & *(za(p2,p5)*zb(p1,p5)*zb(p4,p6)
     &  -za(p2,p3)*zb(p1,p6)*zb(p3,p4))
     & *iza(p2,p6)*izb(p1,p5)*izb(p2,p5)

     & -iden(s345)*Qu*za(p2,p3)*zb(p1,p4)*zab2(p2,p1,p6,p4)
     & *iza(p1,p6)*iza(p2,p6)*izb(p1,p5)*izb(p4,p5)

     & -iden(s346)*Qd*zab2(p3,p2,p5,p1)*zab2(p2,p4,p6,p1)
     & *iza(p2,p6)*iza(p4,p6)*izb(p1,p5)*izb(p2,p5)

     & -iden(s34)*iden(s3456)*za(p2,p3)*za(p2,p5)*zb(p1,p4)*zb(p1,p6)
     & *iza(p2,p6)*izb(p1,p5)

     & +iden(s3456)*za(p2,p3)
     & *(za(p2,p3)*za(p4,p6)*zb(p4,p6)*zb(p1,p3)*zb(p1,p4)
     & +zab2(p2,p4,p6,p1)*zb(p4,p5)*zab2(p5,p4,p6,p1))
     & /s456*iza(p2,p6)*izb(p1,p5)*izb(p4,p5)*iza(p4,p6)

     & +iden(s345)*iden(s3456)*za(p2,p3)*zb(p1,p4)*zb(p1,p6)
     & *zab2(p2,p1,p6,p4)*iza(p2,p6)*izb(p1,p5)*izb(p4,p5)

     &+iden(s34)*(
     & -zab2(p2,p1,p6,p4)*zab2(p3,p2,p5,p1)*Qu*Qd
     & +za(p2,p6)*za(p2,p3)*zb(p1,p6)*zb(p2,p5)*zab2(p5,p1,p6,p4)*Qu**2
     & /s156
     & +za(p2,p5)*za(p1,p6)*zb(p1,p4)*zb(p1,p5)*zab2(p3,p2,p5,p6)*Qd**2
     & /s256
     & )*iza(p1,p6)*iza(p2,p6)*izb(p1,p5)*izb(p2,p5)

      lo(2,1)=
     & -iden(s34)*iden(s346)*iden(s3456)*zab2(p2,p3,p6,p4)*zb(p1,p5)
     & *(za(p3,p6)*za(p2,p5)*zb(p1,p5)+za(p2,p6)*za(p3,p4)*zb(p1,p4)
     & )*iza(p2,p5)*izb(p1,p6)

     & + Qu*iden(s34)*iden(s346)*za(p4,p5)*(
     & - za(p1,p2)*za(p2,p3)*zab2(p6,p3,p4,p1)*zb(p1,p4)
     & + zab2(p2,p3,p4,p1)*za(p2,p5)*za(p3,p6)*zb(p4,p5)
     & + zab2(p2,p1,p5,p4)*za(p2,p6)*za(p3,p6)*zb(p1,p6)
     & + za(p2,p5)*za(p2,p6)*za(p3,p4)*zb(p1,p4)*zb(p4,p5)
     & )*iza(p1,p5)*iza(p2,p5)*izb(p1,p6)*iza(p4,p5)

     & + iden(s345)*iden(s3456)*za(p2,p6)*(
     & +za(p2,p6)*zab2(p3,p4,p5,p1)*zb(p1,p6)*zb(p4,p6)
     & -za(p2,p3)*zab2(p3,p4,p5,p6)*zb(p1,p3)*zb(p1,p4)
     & -za(p2,p3)*za(p3,p5)*zb(p1,p3)*zb(p1,p6)*zb(p4,p5)
     & )*iza(p2,p5)*izb(p1,p6)*iza(p4,p5)*izb(p4,p6)

     & +iden(s34)*iden(s345)*iden(s3456)*zab2(p3,p2,p6,p1)
     & *(zab2(p2,p4,p5,p1)*zb(p4,p5)-za(p2,p3)*zb(p1,p4)*zb(p3,p5))
     & *iza(p2,p5)*izb(p1,p6)*za(p2,p6)

     & + Qd*iden(s34)*iden(s345)*zab2(p3,p2,p6,p1)
     & *(za(p2,p6)*zb(p1,p6)*zb(p4,p5)
     &  -za(p2,p3)*zb(p1,p5)*zb(p3,p4))
     & *iza(p2,p5)*izb(p1,p6)*izb(p2,p6)

     & -iden(s346)*Qu*za(p2,p3)*zb(p1,p4)*zab2(p2,p1,p5,p4)
     & *iza(p1,p5)*iza(p2,p5)*izb(p1,p6)*izb(p4,p6)

     & -iden(s345)*Qd*zab2(p3,p2,p6,p1)*zab2(p2,p4,p5,p1)
     & *iza(p2,p5)*iza(p4,p5)*izb(p1,p6)*izb(p2,p6)

     & -iden(s34)*iden(s3456)*za(p2,p3)*za(p2,p6)*zb(p1,p4)*zb(p1,p5)
     & *iza(p2,p5)*izb(p1,p6)

     & +iden(s3456)*za(p2,p3)
     & *(za(p2,p3)*za(p4,p5)*zb(p4,p5)*zb(p1,p3)*zb(p1,p4)
     & +zab2(p2,p4,p5,p1)*zb(p4,p6)*zab2(p6,p4,p5,p1))
     & /s456*iza(p2,p5)*izb(p1,p6)*izb(p4,p6)*iza(p4,p5)

     & +iden(s346)*iden(s3456)*za(p2,p3)*zb(p1,p4)*zb(p1,p5)
     & *zab2(p2,p1,p5,p4)*iza(p2,p5)*izb(p1,p6)*izb(p4,p6)

     &+iden(s34)*(
     & -zab2(p2,p1,p5,p4)*zab2(p3,p2,p6,p1)*Qu*Qd
     & +za(p2,p5)*za(p2,p3)*zb(p1,p5)*zb(p2,p6)*zab2(p6,p1,p5,p4)*Qu**2
     & /s156
     & +za(p2,p6)*za(p1,p5)*zb(p1,p4)*zb(p1,p6)*zab2(p3,p2,p6,p5)*Qd**2
     & /s256
     & )*iza(p1,p5)*iza(p2,p5)*izb(p1,p6)*izb(p2,p6)

      msqWaa=      
     &  abs(lo(1,1))**2+abs(lo(1,2))**2
     & +abs(lo(2,1))**2+abs(lo(2,2))**2
     
      return
      end
