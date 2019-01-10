      double complex
     & FC0save(Ncmax,-2:0),
     & FC1save(Ncmax,y1max,-2:0),
     & FC2save(Ncmax,y2max,-2:0),
     & FC3save(Ncmax,y3max,-2:0),
     & C00save(Ncmax,-2:0),tau3save(Ncmax,4,-2:0)
      common/ovCsave/FC0save,FC1save,FC2save,FC3save,
     & C00save,tau3save
!$omp threadprivate(/ovCsave/)      
     
      
