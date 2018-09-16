      function C0fa2m(t,qsq,msq)
      implicit none
      include 'types.f'
      complex(dp):: C0fa2m
C     C0(Pc,Pg,0,msq,msq)=
C     C0(tsq,0,qsq,0,msq,msq) (LT notation) 
C     result for qsq<0,t<0 is 
C     C0fa2m(t,qsq,msq)=(Li2(qsq/msq)-Li2(t/msq))/(t-qsq)

       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: t,qsq,msq,r,omr,ddilog
      complex(dp):: lnrat,wlog,dilogt,dilogq

      r=one-qsq/msq
      omr=qsq/msq
      wlog=lnrat(msq-qsq,msq)
      if (omr > one) then 
         dilogq=cplx1(pisqo6-ddilog(r))-wlog*cplx1(log(omr))
      else
         dilogq=cplx1(ddilog(omr))
      endif

      r=one-t/msq
      omr=t/msq
      wlog=lnrat(msq-t,msq)
      if (omr > one) then 
         dilogt=cplx1(pisqo6-ddilog(r))-wlog*cplx1(log(omr))
      else
         dilogt=cplx1(ddilog(omr))
      endif
      C0fa2m=(dilogq-dilogt)/(t-qsq)
      return
      end
