      double complex function C0fa2m(t,qsq,msq)
C     C0(Pc,Pg,0,msq,msq)=
C     C0(tsq,0,qsq,0,msq,msq) (LT notation) 
C     result for qsq<0,t<0 is 
C     C0fa2m(t,qsq,msq)=(Li2(qsq/msq)-Li2(t/msq))/(t-qsq)

      implicit none 
      include 'constants.f'
      double precision t,qsq,msq,r,omr,ddilog
      double complex lnrat,wlog,dilogt,dilogq

      r=1d0-qsq/msq
      omr=qsq/msq
      wlog=lnrat(msq-qsq,msq)
      if (omr .gt. one) then 
         dilogq=dcmplx(pisqo6-ddilog(r))-wlog*dcmplx(log(omr))
      else
         dilogq=dcmplx(ddilog(omr))
      endif

      r=1d0-t/msq
      omr=t/msq
      wlog=lnrat(msq-t,msq)
      if (omr .gt. one) then 
         dilogt=dcmplx(pisqo6-ddilog(r))-wlog*dcmplx(log(omr))
      else
         dilogt=dcmplx(ddilog(omr))
      endif
      C0fa2m=(dilogq-dilogt)/(t-qsq)
      return
      end
