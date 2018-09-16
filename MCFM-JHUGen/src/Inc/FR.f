!---------------------------------------------
! FR contains regular and plus pieces associated
! with radiation from of a photon by a quark qu,or antiquark
! qub
      integer:: quf,qubf
      parameter(quf=1,qubf=-1)

      real(dp):: FR(-1:1,-1:1,1:3)

      common/FR/FR
