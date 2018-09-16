c--- Common block that identifies the heavy flavour being used
c--- in the current process
      integer:: flav
      common/heavyflav/flav   
!$omp threadprivate(/heavyflav/)
