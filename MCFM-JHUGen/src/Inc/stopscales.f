      double precision initfacscale_H,initfacscale_L,
     .                 initrenscale_H,initrenscale_L
      double precision facscale_H,facscale_L,
     .                 renscale_H,renscale_L
      double precision msqLH(-nf:nf,-nf:nf),msqHL(-nf:nf,-nf:nf)
      double precision as_H,as_L
      common/stopscales/initfacscale_H,initfacscale_L,
     .                  initrenscale_H,initrenscale_L,
     .                  facscale_H,facscale_L,
     .                  renscale_H,renscale_L,
     .                  as_H,as_L,
     .                  msqLH,msqHL
!$omp threadprivate(/stopscales/)      
