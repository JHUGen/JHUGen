c--- ancillary array that is used for fragmentation
c--- contributions in gmgmjt process
      real(dp):: msqbits(12)
      common/msqbits/msqbits
      integer::
     & ddb_ddb,ddb_ssb,ddb_uub,
     & uub_uub,uub_ccb,uub_ddb,
     & dbd_ddb,dbd_ssb,dbd_uub,
     & ubu_uub,ubu_ccb,ubu_ddb
      parameter(
     & ddb_ddb=1,ddb_ssb=2,ddb_uub=3,
     & uub_uub=4,uub_ccb=5,uub_ddb=6,
     & dbd_ddb=7,dbd_ssb=8,dbd_uub=9,
     & ubu_uub=10,ubu_ccb=11,ubu_ddb=12)

!$omp threadprivate(/msqbits/)

