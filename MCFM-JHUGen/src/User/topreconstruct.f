      subroutine topreconstruct(q,
     &  failed,ptt,yt,pttb,ytb,ptttb,yttb,mttb,mwp,mwm)
      implicit none
      include 'types.f'
c--- Given the usual momentum array, try to reconstruct which
c--- set of momenta reconstruct the top and antitop quark;
c--- given those assignments, the routine returns the following quantities:
c---
c---    ptt, yt            pt and rapidity of top quark
c---    pttb, ytb      pt and rapidity of anti-top quark
c---    ptttb, yttb      pt and rapidity of (top,anti-top) system
c---    mttb            invariant mass of (top,anti-top) system
c---    mwp,wwm            invariant masses of W+ and W- candidates
c---
c--- For real radiation events, the jet algorithm may result in 
c--- events for which the invariant mass is not exactly equal to mt,
c--- despite the phase space producing tops exactly on-shell. In that
c--- case, allow |s-mt| <= 'toler' [GeV]  
c---
c--- If no consistent reconstruction is found (i.e. no combinations
c--- satisfying above condition), then failed is equal to .true.
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'plabel.f'
      logical:: failed
      integer:: iorder(mxpart),i,j
      real(dp):: tiny,check678,check6789,check345,check3459,
     & ptthree,yrapthree,ptfour,yrapfour,ptsix,yrapsix,yrapseven,
     & ptt,yt,pttb,ytb,ptttb,yttb,mttb,mwp,mwm,p(mxpart,4),toler,
     & mwtest(mxpart),dot,q(mxpart,4)
      parameter (tiny=1.e-8_dp,toler=5d9)

c--- copy array from p to q
c--- q may now be reshuffled to find best W candidate, if necessary      
      do i=1,9
      do j=1,4
      p(i,j)=q(i,j)
      enddo
      enddo

c--- default value
      iorder(4)=-1
c--- hadronic decay of top      
      if ((plabel(3) == 'pp') .and. (p(9,4) > tiny)) then
c--- W- candidate could be 34,39,49 or 349
        mwtest(1)=2._dp*dot(p,3,4) 
        mwtest(2)=2._dp*dot(p,3,9) 
        mwtest(3)=2._dp*dot(p,4,9) 
        mwtest(4)=mwtest(1)+mwtest(2)+mwtest(3)
      do j=1,4
        mwtest(j)=abs(mwtest(j)-wmass**2)
        iorder(j)=j
      enddo
      call arraysort(4,mwtest,iorder)
      if (iorder(4) == 2) then
c----- candidate is (3,9) -> swap 4,9
           do j=1,4
             p(4,j)=q(9,j)
             p(9,j)=q(4,j)
           enddo
      endif
      if (iorder(4) == 3) then
c----- candidate is (4,9) -> swap 3,9
           do j=1,4
             p(3,j)=q(9,j)
             p(9,j)=q(3,j)
           enddo
      endif
      endif
      if ((plabel(3) == 'pp') .and. (iorder(4) == 4)) then
      mwp=sqrt((p(3,4)+p(4,4)+p(9,4))**2-(p(3,1)+p(4,1)+p(9,1))**2         
     &           -(p(3,2)+p(4,2)+p(9,2))**2-(p(3,3)+p(4,3)+p(9,3))**2)
      else
      mwp=sqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2         
     &           -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)
      endif      

c--- hadronic decay of antitop      
      if ((plabel(7) == 'pp') .and. (p(9,4) > tiny)) then
c--- W- candidate could be 78,79,89 or 789
        mwtest(1)=2._dp*dot(p,7,8) 
        mwtest(2)=2._dp*dot(p,7,9) 
        mwtest(3)=2._dp*dot(p,8,9) 
        mwtest(4)=mwtest(1)+mwtest(2)+mwtest(3)
      do j=1,4
        mwtest(j)=abs(mwtest(j)-wmass**2)
        iorder(j)=j
      enddo
      call arraysort(4,mwtest,iorder)
      if (iorder(4) == 2) then
c----- candidate is (7,9) -> swap 8,9
           do j=1,4
             p(8,j)=q(9,j)
             p(9,j)=q(8,j)
           enddo
      endif
      if (iorder(4) == 3) then
c----- candidate is (8,9) -> swap 7,9
           do j=1,4
             p(7,j)=q(9,j)
             p(9,j)=q(7,j)
           enddo
      endif
      endif
      if ((plabel(7) == 'pp') .and. (iorder(4) == 4)) then
      mwm=sqrt((p(7,4)+p(8,4)+p(9,4))**2-(p(7,1)+p(8,1)+p(9,1))**2         
     &           -(p(7,2)+p(8,2)+p(9,2))**2-(p(7,3)+p(8,3)+p(9,3))**2)
      else
      mwm=sqrt((p(7,4)+p(8,4))**2-(p(7,1)+p(8,1))**2         
     &           -(p(7,2)+p(8,2))**2-(p(7,3)+p(8,3))**2)
      endif      
      
c--- default: reconstruction is okay
      failed=.false.      
      
      check345=abs(sqrt(
     & (p(3,4)+p(4,4)+p(5,4))**2      
     &-(p(3,1)+p(4,1)+p(5,1))**2      
     &-(p(3,2)+p(4,2)+p(5,2))**2      
     &-(p(3,3)+p(4,3)+p(5,3))**2) - mt ) 
      check678=abs(sqrt(
     & (p(6,4)+p(7,4)+p(8,4))**2      
     &-(p(6,1)+p(7,1)+p(8,1))**2      
     &-(p(6,2)+p(7,2)+p(8,2))**2      
     &-(p(6,3)+p(7,3)+p(8,3))**2)  -mt )
     
      if ((plabel(3) == 'pp') .and. (iorder(4) == 4)) then
c--- W candidate is 349, so top=345 is not acceptable
        check345=1d9
      endif       
     
      if ((plabel(7) == 'pp') .and. (iorder(4) == 4)) then
c--- W candidate is 789, so top=678 is not acceptable
        check678=1d9
      endif       
     
      if (p(9,4) > tiny) then
c--- these are also possible candidates
        check3459=abs(sqrt(
     &   (p(3,4)+p(4,4)+p(5,4)+p(9,4))**2      
     &  -(p(3,1)+p(4,1)+p(5,1)+p(9,1))**2      
     &  -(p(3,2)+p(4,2)+p(5,2)+p(9,2))**2      
     &  -(p(3,3)+p(4,3)+p(5,3)+p(9,3))**2) - mt)
        check6789=abs(sqrt(
     &   (p(6,4)+p(7,4)+p(8,4)+p(9,4))**2      
     &  -(p(6,1)+p(7,1)+p(8,1)+p(9,1))**2      
     &  -(p(6,2)+p(7,2)+p(8,2)+p(9,2))**2      
     &  -(p(6,3)+p(7,3)+p(8,3)+p(9,3))**2) - mt )
      else
c--- these values will automatically fail the next checks
        check3459=1d9
        check6789=1d9
      endif
      
      if     ((check345 < toler) .and. (check678 < toler)
     &  .and. (check345+check678 < check345+check6789)
     &  .and. (check345+check678 < check3459+check678)) then
c--- top = 345, anti-top = 678
        ptt=ptthree(3,4,5,p)
        yt=yrapthree(3,4,5,p)
        pttb=ptthree(6,7,8,p)
        ytb=yrapthree(6,7,8,p)
        ptttb=ptsix(3,4,5,6,7,8,p)
       yttb=yrapsix(3,4,5,6,7,8,p)
      mttb=(p(3,4)+p(4,4)+p(5,4)+p(6,4)+p(7,4)+p(8,4))**2       
     &      -(p(3,1)+p(4,1)+p(5,1)+p(6,1)+p(7,1)+p(8,1))**2       
     &      -(p(3,2)+p(4,2)+p(5,2)+p(6,2)+p(7,2)+p(8,2))**2       
     &      -(p(3,3)+p(4,3)+p(5,3)+p(6,3)+p(7,3)+p(8,3))**2       
        mttb=sqrt(max(mttb,zip))
      mwm=sqrt((p(7,4)+p(8,4))**2-(p(7,1)+p(8,1))**2         
     &           -(p(7,2)+p(8,2))**2-(p(7,3)+p(8,3))**2)
      elseif ((check3459 < toler) .and. (check678 < toler)
     &  .and. (check3459+check678 < check345+check678)
     &  .and. (check3459+check678 < check345+check6789)) then
c--- top = 3459, anti-top = 678
        ptt=ptfour(3,4,5,9,p)
        yt=yrapfour(3,4,5,9,p)
        pttb=ptthree(6,7,8,p)
        ytb=yrapthree(6,7,8,p)
      ptttb=zip
      yttb=yrapseven(3,4,5,6,7,8,9,p)
      mttb=(p(3,4)+p(4,4)+p(5,4)+p(6,4)+p(7,4)+p(8,4)+p(9,4))**2      
     &      -(p(3,1)+p(4,1)+p(5,1)+p(6,1)+p(7,1)+p(8,1)+p(9,1))**2      
     &      -(p(3,2)+p(4,2)+p(5,2)+p(6,2)+p(7,2)+p(8,2)+p(9,2))**2      
     &      -(p(3,3)+p(4,3)+p(5,3)+p(6,3)+p(7,3)+p(8,3)+p(9,3))**2      
      mttb=sqrt(max(mttb,zip))
      mwm=sqrt((p(7,4)+p(8,4))**2-(p(7,1)+p(8,1))**2         
     &           -(p(7,2)+p(8,2))**2-(p(7,3)+p(8,3))**2)
      elseif ((check345 < toler) .and. (check6789 < toler)
     &  .and. (check345+check6789 < check345+check678)
     &  .and. (check345+check6789 < check3459+check678)) then
c--- top = 345, anti-top = 6789
        ptt=ptthree(3,4,5,p)
        yt=yrapthree(3,4,5,p)
        pttb=ptfour(6,7,8,9,p)
        ytb=yrapfour(6,7,8,9,p)
        ptttb=zip
        yttb=yrapseven(3,4,5,6,7,8,9,p)
        mttb=(p(3,4)+p(4,4)+p(5,4)+p(6,4)+p(7,4)+p(8,4)+p(9,4))**2      
     &      -(p(3,1)+p(4,1)+p(5,1)+p(6,1)+p(7,1)+p(8,1)+p(9,1))**2      
     &      -(p(3,2)+p(4,2)+p(5,2)+p(6,2)+p(7,2)+p(8,2)+p(9,2))**2      
     &      -(p(3,3)+p(4,3)+p(5,3)+p(6,3)+p(7,3)+p(8,3)+p(9,3))**2      
        mttb=sqrt(max(mttb,zip))
        mwm=sqrt((p(7,4)+p(8,4))**2-(p(7,1)+p(8,1))**2         
     &           -(p(7,2)+p(8,2))**2-(p(7,3)+p(8,3))**2)
      else
c--- reconstruction failed: set flag and use out-of-range returns
        failed=.true.
        ptt=-one
        yt=100._dp
        pttb=-one
        ytb=100._dp
        ptttb=-one
        yttb=100._dp
        mttb=-one
        mwp=-one
        mwm=-one
      endif      
    
      return
      end
      
