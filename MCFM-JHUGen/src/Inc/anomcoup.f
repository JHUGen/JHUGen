      real(dp):: delg1_z,delg1_g,lambda_g,lambda_z,
     & h1Z,h2Z,h3Z,h4Z,h1gam,h2gam,h3gam,h4gam,
     & delk_g,delk_z,tevscale
      real(dp):: h1tZ,h2tZ,h3tZ,h4tZ,h1tgam,h2tgam,h3tgam,h4tgam
      logical:: anomtgc
      common/anomcoup/delg1_z,delg1_g,lambda_g,lambda_z,delk_g,delk_z,
     & h1Z,h2Z,h3Z,h4Z,h1gam,h2gam,h3gam,h4gam,
     & tevscale,anomtgc
!$omp threadprivate(/anomcoup/)

