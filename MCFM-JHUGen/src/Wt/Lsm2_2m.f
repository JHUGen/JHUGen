      double complex function Lsm2_2m(t13bar,t23bar,qsq,msq)
      implicit none
      include 'constants.f'
      double precision qsq,msq,ddilog,r,omr,t13bar,t23bar
      double complex lnrat,dilogt,dilogu,wlog,dilog1,dilog2,dilogq 

C     This function is equal to (in looptools notation)
c      Lsm2_2m_ff(t13-msq,t23-msq,qsq,msq)=
c     . (t13-msq)*(t23-msq)*D0(0,msq,0,qsq;t13,t23;0,0,msq,msq)
c     . -(t13-msq)*C0(0,msq,t13;0,0,msq)
c     . -(t23-qsq)*C0(0,t23,qsq;0,0,msq)

C     This function is equal to (in normal notation) Pg^2=Ps^2=0,Pc^2=m^2
c      Lsm2_2m_ff(2*Pg.Pc,2*Ps.Pc,qsq,msq)=
c     . (2*Pg.Pc)*(2*Ps.Pc)*D0(Pg,Pc,Ps,0,0,msq,msq)
c     . -(2*Pg.Pc)*C0(Pg,Pc,0,0,msq)
c     . +(2*Pg.(Pc+Ps))*C0(Pg,Pc+Ps;0,0,msq)

C     This function is equal to (in normal notation) Pg^2=Ps^2=0,Pc^2=m^2
c      D0(Pg,Pc,Ps,0,0,msq,msq)=
C      +C0(Pg,Pc,0,0,msq)/2/Ps.Pc
C      -C0(Pg,Pc+Ps;0,0,msq)*Pg.(Pc+Ps)/2/Pg.Pc/Ps.Pc
C      +Lsm2_2m_ff(2*Pg.Pc,2*Ps.Pc,qsq,msq)/4/Pg.Pc/Ps.Pc
      
      

c      dilogt ~ ddilog(t13/msq))
c      dilogu ~ ddilog(t23/msq))
c      dilogq ~ ddilog(qsq/msq)
c      dilog1 ~ ddilog((t13-qsq)/(t13-msq))
c      dilog2 ~ ddilog((t23-qsq)/(t23-msq))
 
C      t13=t13bar+msq 
C      t23=t23bar+msq 

      r=-t13bar/msq
      omr=1d0-r
      wlog=lnrat(-t13bar,msq)
      if (omr .gt. one) then 
         dilogt=dcmplx(pisqo6-ddilog(r))-wlog*dcmplx(log(omr))
      else
         dilogt=dcmplx(ddilog(omr))
      endif

      r=-t23bar/msq
      omr=1d0-r
      wlog=lnrat(-t23bar,msq)
      if (omr .gt. one) then 
         dilogu=dcmplx(pisqo6-ddilog(r))-wlog*dcmplx(log(omr))
      else
         dilogu=dcmplx(ddilog(omr))
      endif

      r=-(msq-qsq)/t13bar
      omr=1d0-r
      wlog=lnrat(msq-qsq,-t13bar)
      if (omr .gt. one) then 
         dilog1=dcmplx(pisqo6-ddilog(r))-wlog*dcmplx(log(omr))
      else
         dilog1=dcmplx(ddilog(omr))
      endif

      r=-(msq-qsq)/t23bar
      omr=1d0-r
      wlog=lnrat(msq-qsq,-t23bar)
      if (omr .gt. one) then 
         dilog2=dcmplx(pisqo6-ddilog(r))-wlog*dcmplx(log(omr))
      else
         dilog2=dcmplx(ddilog(omr))
      endif

      r=(msq-qsq)/msq
      omr=1d0-r
      wlog=lnrat(msq-qsq,msq)
      if (omr .gt. one) then 
         dilogq=dcmplx(pisqo6-ddilog(r))-wlog*dcmplx(log(omr))
      else
         dilogq=dcmplx(ddilog(omr))
      endif


      Lsm2_2m=-lnrat(-t13bar,-t23bar)**2-dcmplx(pisqo6)
     . +dilogq-dilogt-dilogu-2d0*(dilog1+dilog2)

      return
      end
