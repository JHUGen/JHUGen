      integer:: jets
      character*2 jetlabel(mxpart)
      common/parts_int/jets
      common/parts_char/jetlabel
!$omp threadprivate(/parts_int/,/parts_char/)      
