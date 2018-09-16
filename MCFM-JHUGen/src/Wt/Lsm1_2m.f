      function Lsm1_2m(s,tbar,qsq,msq)
      implicit none
      include 'types.f'
      complex(dp):: Lsm1_2m
C this corresponds to (in looptools notation)
c      Lsm1_2m=
c     &  s*(t-msq)*D0(0,0,msq,qsq,s,t,0,0,0,msq)
c     & -s*C0(0,s,0,0,0,0)
c     & -(t-qsq)*C0(0,t,qsq,0,0,msq)
c     & -(t-msq)*C0(0,msq,t,0,0,msq)
C this corresponds to (in normal notation) Pg^2=0,Ps^2=0,Pc^2=m^2
c      Lsm1_2m=
c     &  2*Pg.Ps*(2*Ps.Pc)*D0(Pg,Ps,Pc,0,0,0,msq)
c     & -2*Pg.Ps*C0(Pg,Ps,0,0,0)
c     & +2*Pg.(Ps+Pc)*C0(Pg,Ps+Pc,0,0,msq)
c     & -(2*Ps.Pc)*C0(Ps,Pc,0,0,msq)

C      D0(Pg,Ps,Pc,0,0,0,msq)=
c      Lsm1_2m(2*Pg.Ps,2*Pc.Ps,qsq,msq)/4/Pg.Ps/Ps.Pc
c     & +C0(Pg,Ps,0,0,0)/2/Ps.Pc
c     & -C0(Pg,Ps+Pc,0,0,msq)*Pg.(Ps+Pc)/2/Pg.Ps/Ps.Pc
c     & +C0(Ps,Pc,0,0,msq)/2/Pg.Ps


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: s,tbar,msq,qsq,ddilog,r,omr
      complex(dp):: lnrat,dilogqt,dilogt,dilogq,wlogt,wlogs
C      t=tbar+msq
      r=-(msq-qsq)/tbar
      omr=one-r
      wlogt=lnrat(msq-qsq,-tbar)
      if (omr > one) then
         dilogqt=cplx1(pisqo6-ddilog(r))-wlogt*cplx1(log(omr))
      else
         dilogqt=cplx1(ddilog(omr))
      endif

      r=(msq-qsq)/msq
      omr=one-r
      wlogt=lnrat(msq-qsq,msq)
      if (omr > one) then
         dilogq=cplx1(pisqo6-ddilog(r))-wlogt*cplx1(log(omr))
      else
         dilogq=cplx1(ddilog(omr))
      endif

      wlogt=lnrat(-tbar,msq)
      wlogs=lnrat(-s,msq)

      r=-tbar/msq
      omr=one-r
      if (omr > one) then
         dilogt=cplx1(pisqo6-ddilog(r))-wlogt*cplx1(log(omr))
      else
         dilogt=cplx1(ddilog(omr))
      endif


      Lsm1_2m=-cplx1(three*pisqo6)+dilogq-two*(dilogqt+dilogt)
     & +two*wlogt*(wlogs-wlogt)-half*wlogs**2

      return
      end
