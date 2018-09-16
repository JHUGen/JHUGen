c--- Grid of z-vaules for fragmentation functions

      integer:: num_z
      parameter(num_z = 45)

      real(dp):: z_grid(num_z)

      data z_grid  /
     & +0.4627e-01_dp,0.5131e-01_dp,0.5690e-01_dp,0.6311e-01_dp,0.6999e-01_dp,
     & +0.7762e-01_dp,
     & +0.8608e-01_dp,0.9547e-01_dp,0.1059e+00_dp,0.1174e+00_dp,0.1302e+00_dp,
     & +0.1444e+00_dp,
     & +0.1602e+00_dp,0.1776e+00_dp,0.1970e+00_dp,0.2185e+00_dp,0.2423e+00_dp,
     & +0.2687e+00_dp,
     & +0.2980e+00_dp,0.3305e+00_dp,0.3666e+00_dp,0.4065e+00_dp,0.4508e+00_dp,
     & +0.5000e+00_dp,
     & +0.5321e+00_dp,0.5643e+00_dp,0.5964e+00_dp,0.6286e+00_dp,0.6607e+00_dp,
     & +0.6929e+00_dp,
     & +0.7250e+00_dp,0.7571e+00_dp,0.7893e+00_dp,0.8214e+00_dp,0.8536e+00_dp,
     & +0.8857e+00_dp,
     & +0.9000e+00_dp,0.9179e+00_dp,0.9300e+00_dp,0.9400e+00_dp,0.9500e+00_dp,
     & +0.9600e+00_dp,
     & +0.9700e+00_dp,0.9800e+00_dp,0.99e+00_dp /

c--- Grid of M**2 Values for fragmentation functions

      integer:: num_M2_4Flav, num_M2_5Flav
      parameter (num_M2_4Flav=8,num_M2_5Flav=22)

      real(dp):: M2_4Flav(num_M2_4Flav), M2_5Flav(num_M2_5Flav)

      data M2_4Flav /
     & 2._dp,2.8_dp,3.9_dp,5.4_dp,7.5_dp,10.5_dp,14.5_dp,20.24_dp/

      data M2_5Flav /
     & +20.26_dp,33._dp,55._dp,93._dp,156._dp,261._dp,
     & +435._dp,732._dp,1224._dp,2048._dp,
     & +3427._dp,5733._dp,9592._dp,16047._dp,
     &  +26847._dp,44915._dp,75144._dp,1.257e+5_dp,
     & +2.1e+5_dp,3.52e+5_dp,5.89e+5_dp,9.848e+5_dp /






