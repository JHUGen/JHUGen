      double complex function Fvf(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      character*9 st
      double complex Lsm1_2mh,Lsm1_2me,I3m,zab2,I3m123456,I3m563412
      double precision t  
C---statement function
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      I3m123456=I3m(s(j1,j2),s(j3,j4),s(j5,j6))
      I3m563412=I3m123456

      if(st.eq.'q+qb-g+g-') then
      Fvf= 
     & -zab2(j5,j2,j4,j1)**2/(za(j5,j6)*zb(j1,j2)*zab2(j3,j1,j2,j4)**2)
     & *Lsm1_2mh(s(j3,j4),t(j1,j2,j4),s(j1,j2),s(j5,j6))
     & -zab2(j2,j1,j3,j6)**2/(za(j1,j2)*zb(j5,j6)*zab2(j3,j1,j2,j4)**2)
     & *Lsm1_2mh(s(j3,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))

     & -I3m123456*za(j4,j5)*zb(j1,j3)
     & *zab2(j2,j1,j3,j6)/(2d0*zab2(j3,j1,j2,j4)*t(j1,j2,j3))

     & -I3m123456*za(j2,j4)*zb(j3,j6)
     & *zab2(j5,j2,j4,j1)/(2d0*zab2(j3,j1,j2,j4)*t(j1,j2,j4))

     & +I3m563412*za(j2,j4)*zb(j3,j6)
     & *zab2(j5,j3,j6,j1)/(2d0*zab2(j3,j1,j2,j4)*t(j3,j5,j6))

     & +I3m563412*za(j4,j5)*zb(j1,j3)
     & *zab2(j2,j4,j5,j6)/(2d0*zab2(j3,j1,j2,j4)*t(j4,j5,j6))
      elseif(st.eq.'q+qb-g+g+') then
      Fvf=-za(j2,j5)**2/(za(j1,j2)*za(j5,j6)*za(j3,j4)**2)
     & *Lsm1_2me(t(j1,j2,j3),t(j1,j2,j4),s(j1,j2),s(j5,j6))
      endif
      return
      end 
