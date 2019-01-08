      integer ci(2),cii(z2max),ciii(z3max),ciiii(z4max),
     . ciiiii(z5max),ciiiiii(z6max),ciiiiiii(z7max),
     . czzi(2),czzii(z2max),czziii(z3max),czziiii(z4max),
     . czzzzi(2),czzzzii(z2max),
     . czziiiii(z5max),czzzziii(z3max),czzzzzzi(3)
      common/carrays/ci,cii,ciii,ciiii,ciiiii,ciiiiii,ciiiiiii,
     . czzi,czzii,czziii,czziiii,czzzzi,czzzzii,
     . czziiiii,czzzziii,czzzzzzi
!$omp threadprivate(/carrays/)
