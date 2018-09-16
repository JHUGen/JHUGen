      logical:: clustering,inclusive
      integer:: jetalgorithm
      character*4 algorithm
      common/clustering/clustering,inclusive
      common/algorithm/algorithm
      common/jetalgorithm/jetalgorithm
      integer, parameter ::
     & kt=1, antikt=2, Rsepcone=3, hqrk=4, noclustering=5
