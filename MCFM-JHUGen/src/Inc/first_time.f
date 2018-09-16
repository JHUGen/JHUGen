      logical:: first_time, new_event
      common/first_time/first_time,new_event
!$omp threadprivate(/first_time/)
