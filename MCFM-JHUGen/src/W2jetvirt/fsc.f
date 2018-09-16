      function Fsc(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: Fsc
      integer:: j1,j2,j3,j4,j5,j6
      character*9 st
      complex(dp):: Fsc1,Fsc2,Fsc3,Fsc4,Fsc5,Fsc6,Fsc7,Fsc8
      if(st=='q+g-g+qb-') then
      Fsc=Fsc1(j1,j2,j3,j4,j5,j6,za,zb)
      elseif(st=='q+g+g-qb-') then
      Fsc=Fsc2(j1,j2,j3,j4,j5,j6,za,zb)
      elseif(st=='q+g+g+qb-') then
      Fsc=Fsc3(j1,j2,j3,j4,j5,j6,za,zb)
      elseif(st=='q+g+qb-g-') then
      Fsc=Fsc4(j1,j2,j3,j4,j5,j6,za,zb)
      elseif(st=='q+g+qb-g+') then
      Fsc=Fsc5(j1,j2,j3,j4,j5,j6,za,zb)
      elseif(st=='q+qb-g-g+') then
      Fsc=Fsc6(j1,j2,j3,j4,j5,j6,za,zb)
      elseif(st=='q+qb-g+g-') then
      Fsc=Fsc7(j1,j2,j3,j4,j5,j6,za,zb)
      elseif(st=='q+qb-g+g+') then
      Fsc=Fsc8(j1,j2,j3,j4,j5,j6,za,zb)
      endif
      return
      end
