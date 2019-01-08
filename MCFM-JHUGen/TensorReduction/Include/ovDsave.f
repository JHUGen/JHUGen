      double complex
     & FD0save(Ndmax,-2:0),
     & FD1save(Ndmax,y1max,-2:0),
     & FD2save(Ndmax,y2max,-2:0),
     & FD3save(Ndmax,y3max,-2:0),
     & FD4save(Ndmax,y4max,-2:0)
      common/ovDsave/FD0save,FD1save,FD2save,FD3save,FD4save
!$omp threadprivate(/ovDsave/)      
     
      
