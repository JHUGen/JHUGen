      function Fcc(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: Fcc
      integer:: j1,j2,j3,j4,j5,j6
      character*9 st
      complex(dp):: Fcc_qpgmgpqm,Fcc_qpgpgmqm,Fcc_qpgpgpqm,
     & Fcc_qpgpqmgm,Fcc_qpgpqmgp,Fcc_qpqmgmgp,Fcc_qpqmgpgm,Fcc_qpqmgpgp

      if     (st=='q+g-g+qb-') then
        Fcc=Fcc_qpgmgpqm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st=='q+g+g-qb-') then
        Fcc=Fcc_qpgpgmqm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st=='q+g+g+qb-') then
        Fcc=Fcc_qpgpgpqm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st=='q+g+qb-g-') then
        Fcc=Fcc_qpgpqmgm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st=='q+g+qb-g+') then
        Fcc=Fcc_qpgpqmgp(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st=='q+qb-g-g+') then
        Fcc=Fcc_qpqmgmgp(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st=='q+qb-g+g-') then
        Fcc=Fcc_qpqmgpgm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st=='q+qb-g+g+') then
        Fcc=Fcc_qpqmgpgp(j1,j2,j3,j4,j5,j6,za,zb)
      else
        write(6,*) 'Error in Fcc: argument st is ',st
        stop
      endif

      return
      end

