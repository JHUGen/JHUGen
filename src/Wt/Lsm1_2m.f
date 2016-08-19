      double complex function Lsm1_2m(s,tbar,qsq,msq)
C this corresponds to (in looptools notation)
c      Lsm1_2m=
c     .  s*(t-msq)*D0(0,0,msq,qsq,s,t,0,0,0,msq)
c     . -s*C0(0,s,0,0,0,0)
c     . -(t-qsq)*C0(0,t,qsq,0,0,msq)
c     . -(t-msq)*C0(0,msq,t,0,0,msq)
C this corresponds to (in normal notation) Pg^2=0,Ps^2=0,Pc^2=m^2
c      Lsm1_2m=
c     .  2*Pg.Ps*(2*Ps.Pc)*D0(Pg,Ps,Pc,0,0,0,msq)
c     . -2*Pg.Ps*C0(Pg,Ps,0,0,0)
c     . +2*Pg.(Ps+Pc)*C0(Pg,Ps+Pc,0,0,msq)
c     . -(2*Ps.Pc)*C0(Ps,Pc,0,0,msq)

C      D0(Pg,Ps,Pc,0,0,0,msq)=
c      Lsm1_2m(2*Pg.Ps,2*Pc.Ps,qsq,msq)/4/Pg.Ps/Ps.Pc
c     . +C0(Pg,Ps,0,0,0)/2/Ps.Pc
c     . -C0(Pg,Ps+Pc,0,0,msq)*Pg.(Ps+Pc)/2/Pg.Ps/Ps.Pc
c     . +C0(Ps,Pc,0,0,msq)/2/Pg.Ps

      implicit none
      include 'constants.f'
      double precision s,tbar,msq,qsq,ddilog,r,omr
      double complex lnrat,dilogqt,dilogt,dilogq,wlogt,wlogs 
C      t=tbar+msq
      r=-(msq-qsq)/tbar
      omr=1d0-r
      wlogt=lnrat(msq-qsq,-tbar)
      if (omr .gt. one) then 
         dilogqt=dcmplx(pisqo6-ddilog(r))-wlogt*dcmplx(log(omr))
      else
         dilogqt=dcmplx(ddilog(omr))
      endif

      r=(msq-qsq)/msq
      omr=1d0-r
      wlogt=lnrat(msq-qsq,msq)
      if (omr .gt. one) then 
         dilogq=dcmplx(pisqo6-ddilog(r))-wlogt*dcmplx(log(omr))
      else
         dilogq=dcmplx(ddilog(omr))
      endif

      wlogt=lnrat(-tbar,msq)
      wlogs=lnrat(-s,msq)

      r=-tbar/msq
      omr=1d0-r
      if (omr .gt. one) then 
         dilogt=dcmplx(pisqo6-ddilog(r))-wlogt*dcmplx(log(omr))
      else
         dilogt=dcmplx(ddilog(omr))
      endif


      Lsm1_2m=-dcmplx(3d0*pisqo6)+dilogq-2d0*(dilogqt+dilogt)
     . +2d0*wlogt*(wlogs-wlogt)-0.5d0*wlogs**2

      return
      end
