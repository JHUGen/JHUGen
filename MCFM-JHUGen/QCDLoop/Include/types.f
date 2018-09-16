c      integer, parameter  :: qp = selected_real_kind(18)
c      integer, parameter  :: dp = selected_real_kind(8)
c      integer, parameter  :: sp = selected_real_kind(4)

      integer, parameter  :: sp = kind(1.0)
c      integer, parameter  :: 
c     & dp = selected_real_kind(2*precision(1.0_sp))
c      integer, parameter  :: 
c     & qp = selected_real_kind(2*precision(1.0_sp))

      integer, parameter  :: 
     & dp = kind(1d0)
      integer, parameter  :: 
c     & qp = selected_real_kind(2*precision(1._dp))
     & qp = kind(1d0)


