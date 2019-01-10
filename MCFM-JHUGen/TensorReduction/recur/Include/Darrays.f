      integer di(3),dii(z2max),diii(z3max),diiii(z4max),
     . diiiii(z5max),diiiiii(z6max),diiiiiii(z7max),
     . dzzi(3),dzzii(z2max),dzziii(z3max),dzziiii(z4max),
     . dzzzzi(3),dzzzzii(z2max),
     . dzziiiii(z5max),dzzzziii(z3max),dzzzzzzi(3)
      common/darrays/di,dii,diii,diiii,diiiii,diiiiii,diiiiiii,
     . dzzi,dzzii,dzziii,dzziiii,dzzzzi,dzzzzii,
     . dzziiiii,dzzzziii,dzzzzzzi
!$omp threadprivate(/darrays/)
