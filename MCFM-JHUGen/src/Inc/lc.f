      integer:: colourchoice
c--- 'colourchoice' allows calculation by colour structure
c--- For Gflag=.true. [QQGG, QQGGG processes]
c--- 1) Only leading colour ( NCF . N )
c--- 2) Only sub-leading ( NCF . 1/N )
c--- 3) Only sub-sub-leading ( NCF . 1/N**3 )  [QQGGG only]
c--- 0) The total
c--- For Qflag=.true. [QQBQQB process]
c--- 1) Only leading colour ( NCF . 1 )
c--- 2) Only sub-leading ( NCF . 1/N ) [Identical quarks only]
c--- 0) The total
      common/ColC/colourchoice
!$omp threadprivate(/ColC/)
