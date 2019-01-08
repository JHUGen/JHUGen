      double complex
     & FB0save(Nbmax,-2:0),
     & FB1save(Nbmax,y1max,-2:0),
     & FB2save(Nbmax,y2max,-2:0),
     & B00save(Nbmax,-2:0)
      common/ovBsave/FB0save,FB1save,FB2save,B00save
!$omp threadprivate(/ovBsave/)      
