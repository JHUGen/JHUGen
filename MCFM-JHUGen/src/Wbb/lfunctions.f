************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998/July 2005/June 2007                                   *
************************************************************************

      double complex function L0(x,y)
      implicit none
      include 'constants.f'
      double complex Lnrat
      double precision x,y,denom
      denom=one-x/y
      if (abs(denom) .lt. 1d-7) then
      L0=-cone-dcmplx(denom*(half+denom/3d0))
      else
      L0=Lnrat(x,y)/dcmplx(denom)
      endif
      return
      end

      double complex function L1(x,y)
      implicit none
      include 'constants.f'
      double precision x,y,denom
      double complex L0
      denom=one-x/y
      if (abs(denom) .lt. 1d-7) then
      L1=-half*cone-dcmplx(denom/3d0*(one+0.75d0*denom))
      else
      L1=(L0(x,y)+cone)/dcmplx(denom)
      endif
      return
      end

      double complex function L2(x,y)
      implicit none
      include 'constants.f'
      double complex Lnrat
      double precision x,y,r,denom
      r=x/y
      denom=one-r
      if (abs(denom) .lt. 1d-7) then
      L2=(dcmplx(10d0)+denom*(dcmplx(15d0)+dcmplx(18d0)*denom))
     . /dcmplx(60d0)
      else
      L2=(Lnrat(x,y)-dcmplx(0.5d0*(r-1d0/r)))/dcmplx(denom)**3
      endif  
      return
      end

      double complex function L0old(x,y)
      implicit none
      include 'constants.f'
      double complex Lnrat
      double precision x,y,denom
      denom=one-x/y
      L0old=Lnrat(x,y)/dcmplx(denom)
      return
      end


      double complex function L1old(x,y)
      implicit none
      include 'constants.f'
      double precision x,y,denom
      double complex L0old
      denom=one-x/y
      L1old=(L0old(x,y)+cone)/dcmplx(denom)
      return
      end

      double complex function Ls0(x1,y1,x2,y2)
      implicit none
      include 'constants.f'
      double precision x1,x2,y1,y2,r1,r2
      double complex Lsm1
      r1=x1/y1
      r2=x2/y2
      Ls0=Lsm1(x1,y1,x2,y2)/dcmplx(one-r1-r2)
      return
      end

      double complex function Ls1(x1,y1,x2,y2)
      implicit none
      include 'constants.f'
      double precision x1,x2,y1,y2,r1,r2
      double complex Ls0,L0
      r1=x1/y1
      r2=x2/y2
      Ls1=(Ls0(x1,y1,x2,y2)+L0(x1,y1)+L0(x2,y2))/dcmplx(one-r1-r2)
      return
      end



      double complex function Lsm1(x1,y1,x2,y2)
      implicit none
      include 'constants.f'
      double precision x1,x2,y1,y2,r1,r2,omr1,omr2,ddilog
      double complex dilog1,dilog2,Lnrat
      r1=x1/y1
      r2=x2/y2
      omr1=one-r1
      omr2=one-r2
      if (omr1 .gt. one) then 
         dilog1=dcmplx(pisqo6-ddilog(r1))-Lnrat(x1,y1)*dcmplx(log(omr1))
      else
         dilog1=dcmplx(ddilog(omr1))
      endif
      if (omr2 .gt. one) then 
         dilog2=dcmplx(pisqo6-ddilog(r2))-Lnrat(x2,y2)*dcmplx(log(omr2))
      else
          dilog2=dcmplx(ddilog(omr2))
      endif
      lsm1=dilog1+dilog2+Lnrat(x1,y1)*Lnrat(x2,y2)-dcmplx(pisqo6)
      return
      end

      double complex function Lsm1_2mh(s,t,m1sq,m2sq)
      implicit none
      include 'constants.f'
      double precision s,t,m1sq,m2sq
      double complex lsm1_2mht,I3m
      Lsm1_2mh=Lsm1_2mht(s,t,m1sq,m2sq)
     & +(half*(s-m1sq-m2sq)+m1sq*m2sq/t)*I3m(s,m1sq,m2sq)
      return
      end

      double complex function Lsm1_2mht(s,t,m1sq,m2sq)
      implicit none
      include 'constants.f'
      double precision s,t,m1sq,m2sq,ddilog,r1,r2,omr1,omr2
      double complex Lnrat,dilog1,dilog2
      r1=m1sq/t
      r2=m2sq/t
      omr1=one-r1
      omr2=one-r2

      if (omr1 .gt. one) then 
      dilog1=dcmplx(pisqo6-ddilog(r1))-Lnrat(-m1sq,-t)*dcmplx(log(omr1))
      else
      dilog1=dcmplx(ddilog(omr1))
      endif
      if (omr2 .gt. one) then 
      dilog2=dcmplx(pisqo6-ddilog(r2))-Lnrat(-m2sq,-t)*dcmplx(log(omr2))
      else
      dilog2=dcmplx(ddilog(omr2))
      endif
      lsm1_2mht=-dilog1-dilog2
     & +half*(Lnrat(-s,-m1sq)*Lnrat(-s,-m2sq)-Lnrat(-s,-t)**2)
      return
      end


      double complex function Lsm1_2me(s,t,m1sq,m3sq)
      implicit none
c---- formula taken from G.~Duplancic and B~Nizic [arXiv:hep-ph/0006249 v2]
c---- Eq 71
c---- Lsm1_2me notation follows from 
c----  Z.~Bern, L.~J.~Dixon and D.~A.~Kosower,
c----  %``Dimensionally regulated pentagon integrals,''
c----  Nucl.\ Phys.\ B {\bf 412}, 751 (1994)
c----  [arXiv:hep-ph/9306240].
c----  %%CITATION = HEP-PH 9306240;%%
c----  Eqs. (I.13)
C---- analytic continuation has been checked by calculating numerically.
      include 'constants.f'
      integer j
      double precision s,t,m1sq,m3sq,ddilog,arg(4),omarg(4),f2me,htheta
      double complex Li2(4),wlog(4)
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(s)=half+half*sign(one,s)

      f2me=(s+t-m1sq-m3sq)/(s*t-m1sq*m3sq)
     
      arg(1)=f2me*s
      arg(2)=f2me*t
      arg(3)=f2me*m1sq
      arg(4)=f2me*m3sq

      do j=1,4
         omarg(j)=one-arg(j)
         wlog(j)=log(abs(arg(j)))    
     .     +impi*dcmplx(htheta(-arg(j))*sign(one,f2me))
         if (omarg(j) .gt. one) then 
             Li2(j)=dcmplx(pisqo6-ddilog(arg(j)))
     .       -wlog(j)*dcmplx(log(omarg(j)))
          else
             Li2(j)=dcmplx(ddilog(omarg(j)))
         endif
      enddo
      Lsm1_2me=Li2(1)+Li2(2)-Li2(3)-Li2(4)

      return
      end

c      double complex function Lsm1_2me(s,t,m1sq,m3sq)
c      implicit none
c      include 'constants.f'
c      integer j
c      double precision s,t,m1sq,m3sq,ddilog,arg(5),omarg(5)
c      double complex Lnrat,Li2(5),wlog(5)
c      arg(1)=m1sq/s
c      wlog(1)=Lnrat(-m1sq,-s)

c      arg(2)=m1sq/t
c      wlog(2)=Lnrat(-m1sq,-t)

c      arg(3)=m3sq/s
c      wlog(3)=Lnrat(-m3sq,-s)

c      arg(4)=m3sq/t
c      wlog(4)=Lnrat(-m3sq,-t)

c      arg(5)=arg(1)*arg(4)
c      wlog(5)=Lnrat(-m1sq,-s)+Lnrat(-m3sq,-t)

c      do j=1,5
c         omarg(j)=one-arg(j)
c         if (omarg(j) .gt. one) then 
c             Li2(j)=dcmplx(pisqo6-ddilog(arg(j)))-wlog(j)*log(omarg(j))
c          else
c             Li2(j)=dcmplx(ddilog(omarg(j)))
c         endif
c      enddo
c      Lsm1_2me=Li2(5)-Li2(1)-Li2(2)-Li2(3)-Li2(4)-half*Lnrat(-s,-t)**2
c      return
c      end

