      double complex function Fcc(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'zprods_decl.f'
      character*9 st
      double complex Fcc_qpgmgpqm,Fcc_qpgpgmqm,Fcc_qpgpgpqm,
     . Fcc_qpgpqmgm,Fcc_qpgpqmgp,Fcc_qpqmgmgp,Fcc_qpqmgpgm,Fcc_qpqmgpgp

      if     (st.eq.'q+g-g+qb-') then
        Fcc=Fcc_qpgmgpqm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st.eq.'q+g+g-qb-') then 
        Fcc=Fcc_qpgpgmqm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st.eq.'q+g+g+qb-') then
        Fcc=Fcc_qpgpgpqm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st.eq.'q+g+qb-g-') then
        Fcc=Fcc_qpgpqmgm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st.eq.'q+g+qb-g+') then
        Fcc=Fcc_qpgpqmgp(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st.eq.'q+qb-g-g+') then
        Fcc=Fcc_qpqmgmgp(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st.eq.'q+qb-g+g-') then
        Fcc=Fcc_qpqmgpgm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st.eq.'q+qb-g+g+') then
        Fcc=Fcc_qpqmgpgp(j1,j2,j3,j4,j5,j6,za,zb)
      else
        write(6,*) 'Error in Fcc: argument st is ',st
        stop
      endif

      return
      end
      
