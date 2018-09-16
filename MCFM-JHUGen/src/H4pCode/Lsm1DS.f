      function Lsm1DS(s,t,msq)
      implicit none
      include 'types.f'
      complex(dp):: Lsm1DS
c--- This is an implementation of Eq. (B.2) from DS, 0906.0008

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: s,t,msq,r1,r2,omr1,omr2,ddilog
      complex(dp):: dilog1,dilog2,Lnrat
      r1=s/msq
      r2=t/msq
      omr1=one-r1
      omr2=one-r2
      if (omr1 > one) then
       dilog1=cplx1(pisqo6-ddilog(r1))
     & -Lnrat(-s,-msq)*cplx1(log(omr1))
      else
       dilog1=cplx1(ddilog(omr1))
      endif
      if (omr2 > one) then
       dilog2=cplx1(pisqo6-ddilog(r2))
     & -Lnrat(-t,-msq)*cplx1(log(omr2))
      else
       dilog2=cplx1(ddilog(omr2))
      endif
      Lsm1DS=dilog1+dilog2+Lnrat(-s,-msq)*Lnrat(-t,-msq)-cplx1(pisqo6)
      return
      end

