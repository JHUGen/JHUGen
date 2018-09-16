      double complex function C0fb2m(u,msq)
C     C0(Pc,Pg,0,msq,msq)=
C     C0(msq,0,u,0,msq,msq) (LT notation) 
C  C0fb2m(u,msq)=(Pi^2/6-Li2[u/m2])/(u-msq) with u=(Pc+Pg)^2;
      implicit none 
      include 'constants.f'
      double precision u,msq,ubar,r,omr,ddilog
      double complex lnrat,wlogu,dilogu
      ubar=u-msq
      r=-ubar/msq
      omr=1d0-r
      if (omr .gt. one) then 
         wlogu=lnrat(-ubar,msq)
         dilogu=dcmplx(pisqo6-ddilog(r))-wlogu*dcmplx(log(omr))
      else
         dilogu=dcmplx(ddilog(omr))
      endif
      C0fb2m=(dcmplx(pisqo6)-dilogu)/ubar
      return
      end
