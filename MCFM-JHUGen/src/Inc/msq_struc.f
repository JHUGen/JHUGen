c--- contains the common block and definitions necessary to
c--- separate matrix elements into their constituent parts - for
c--- example, by species of parton in the final state or by
c--- colour orderings (or both)
      real(dp):: msq_struc(8,-nf:nf,-nf:nf),
     &                msq_strucv(8,-nf:nf,-nf:nf)
      integer:: igg_ab,igg_ba,igg_sym,iqq_a,iqq_b,iqq_i,
     & igggg_a,igggg_b,igggg_c,iqr
      common/msq_struc/msq_struc,msq_strucv
!$omp threadprivate(/msq_struc/)
      parameter(igg_ab=4,igg_ba=5,igg_sym=6)
c--- Note that the 4-quark and 4-gluon pieces are never needed simultaneously,
c--- so we can reuse the same indices to save memory
      parameter(iqq_a=1,iqq_b=2,iqq_i=3)
      parameter(igggg_a=1,igggg_b=2,igggg_c=3)
c--- One extra parameter needed for non-identical quark pieces
      parameter(iqr=7)

