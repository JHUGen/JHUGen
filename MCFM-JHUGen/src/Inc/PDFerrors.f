      logical:: PDFerrors
      integer:: maxPDFsets,currentPDF
c--- 40 is my choice for the maximum number of dipoles
c--- 50 is the old choice for the maximum number of PDF error sets,
c--- that has been increased to 1000 in order to include NNPDF sets
      real(dp):: PDFxsec(0:1000),PDFxsec_nd(0:1000,0:40),
     & PDFwgt(0:1000)
      common/PDFerrors/PDFerrors,maxPDFsets,PDFxsec
      common/PDFweight/PDFwgt  
!$omp threadprivate(/PDFweight/)
