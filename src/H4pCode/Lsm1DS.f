      double complex function Lsm1DS(s,t,msq)
c--- This is an implementation of Eq. (B.2) from DS, 0906.0008
      implicit none
      include 'constants.f'
      double precision s,t,msq,r1,r2,omr1,omr2,ddilog
      double complex dilog1,dilog2,Lnrat
      r1=s/msq
      r2=t/msq
      omr1=one-r1
      omr2=one-r2
      if (omr1 .gt. one) then 
       dilog1=dcmplx(pisqo6-ddilog(r1))
     & -Lnrat(-s,-msq)*dcmplx(log(omr1))
      else
       dilog1=dcmplx(ddilog(omr1))
      endif
      if (omr2 .gt. one) then 
       dilog2=dcmplx(pisqo6-ddilog(r2))
     & -Lnrat(-t,-msq)*dcmplx(log(omr2))
      else
       dilog2=dcmplx(ddilog(omr2))
      endif
      Lsm1DS=dilog1+dilog2+Lnrat(-s,-msq)*Lnrat(-t,-msq)-dcmplx(pisqo6)
      return
      end

