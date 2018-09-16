      subroutine stop_def(xsqV,xsqR,q,mc,ms)
      implicit none
      include 'types.f'
c--- Subroutine to be run for each event to calculate log's, bubbles,
c--- triangles and boxes (A0,B0,C0,C00,C001,C002,._dp0) and fill the
c--- common block in stopf1inc.f
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: xsqV,xsqR,q(mxpart,4),ms,mc
      real(dp):: muR,ms2,mc2,xsn,xsd,s,t,u,qsq,dot,cDs
      include 'stopf1inc.f'
      complex(dp):: lnrat,LLs1,LLs2,tr1f,tr2f,tr3,tr3c00f,tr3c001f
     .,tr3c002f,tr3s00f,tr3s001f,tr3s002f,tr4,tr5,B0xf,Bfun
      external lnrat,LLs1,LLs2,tr1f,tr2f,tr3,tr3c00f,tr3c001f
     .,tr3c002f,tr3s00f,tr3s001f,tr3s002f,tr4,tr5,B0xf,Bfun

      ms2 = ms**2
      mc2 = mc**2
      cDs = dot(q,3,4)+mc2*dot(q,4,2)/2._dp/dot(q,3,2)
     &  +ms2*dot(q,3,2)/2._dp/dot(q,4,2)
      qsq = ms2+mc2+2._dp*cDs+2._dp*dot(q,3,2)+2._dp*dot(q,4,2)
      s   = ms2+2._dp*dot(q,4,2)
      t   = mc2+2._dp*dot(q,3,2)
      u   = mc2+ms2+2._dp*cDs
      xsn = (1._dp-sqrt(1._dp-4._dp*ms*mc/(u-(ms-mc)**2)))
      xsd = (1._dp+sqrt(1._dp-4._dp*ms*mc/(u-(ms-mc)**2)))
      muR = sqrt(xsqR)

      lRc1      = lnrat(mc*muR,mc2-t)
      lRs1      = lnrat(ms*muR,ms2-s)
      lRc2      = lRc1**2
      lRs2      = lRs1**2
      lRcs      = lnrat(xsn,-xsd)*lnrat(xsqR,mc*ms)
      lc        = lnrat(mc2-t,mc2)
      ls        = lnrat(ms2-s,ms2)
      lp        = lnrat(xsn,-xsd)
      lVc       = lnrat(xsqV,mc2)
      lVs       = lnrat(xsqV,ms2)
      B0csf     = B0xf(xsqV,u,ms2,mc2)
      B0cgsf    = B0xf(xsqV,qsq,ms2,mc2)
      BfunX     = Bfun(xsqV,qsq,ms2,mc2)
      tr1Xfc    = tr1f(xsqR,t,mc2)
      tr1Xfs    = tr1f(xsqR,s,ms2)
      tr2fu     = tr2f(xsqR,u,ms2,mc2)
      tr3Xc     = tr3(s,qsq,mc2,ms2)
      tr3Xs     = tr3(t,qsq,ms2,mc2)
      tr3s00ft  = tr3s00f(xsqV,t,qsq,ms2,mc2)
      tr3s001ft = tr3s001f(xsqV,t,qsq,ms2,mc2)
      tr3s002ft = tr3s002f(xsqV,t,qsq,ms2,mc2)
      tr3c00fs  = tr3c00f(xsqV,s,qsq,mc2,ms2)
      tr3c001fs = tr3c001f(xsqV,s,qsq,mc2,ms2)
      tr3c002fs = tr3c002f(xsqV,s,qsq,mc2,ms2)
      tr4Xc     = tr4(t,mc2)
      tr4Xs     = tr4(s,ms2)
      tr5Xc     = tr5(u,qsq,mc2,ms2)
      tr5Xs     = tr5(u,qsq,ms2,mc2)
      LsA       = LLs1(xsqR,s,t,qsq,ms2,mc2)
      LsB1      = LLs2(xsqR,s,u,qsq,mc2,ms2)
      LsB2      = LLs2(xsqR,t,u,qsq,ms2,mc2)

      return
      end



