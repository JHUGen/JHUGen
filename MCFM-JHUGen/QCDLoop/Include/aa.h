        logical awrite,atest,aderiv
        common /aaflag/ awrite,atest,aderiv
!$omp threadprivate(/aaflag/)
