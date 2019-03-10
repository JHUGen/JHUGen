      subroutine ovCtensor(p1,p2,xm1s,xm2s,xm3s,
     & FC0,FC1,FC2,FC3,C00,tau3)
      implicit none
C     p1,p2 are the external momenta,
C     m1s,m2s,m3s are the squares of the internal masses
C     FC0...FC3 are the rank 0,...3 triangle functions
C     Lorentz indices are stored as linear array, 
C     thus FC2(y2(n1,n2),ep), etc
C     Author: R.K.Ellis (January 2013)
      include 'TRconstants.f'
      include 'TRydef.f'
      include 'TRscale.f'
      include 'TRbadpoint.f'
      include 'TRmaxindex.f'
      include 'TRonshellcutoff.f'
      include 'TRclear.f'
      include 'ovCnames.f'
      include 'ovCsave.f'
      double precision p1(4),p2(4),p3(4),p12(4),xm1s,xm2s,xm3s,
     & m1s,m2s,m3s,Gram2,s1Dp1,s1Dp2,s2Dp2,p1Dp1,p2Dp2,p3Dp3,p1Dp2,
     & ovw2,d,Gram(2,2),inGram(2,2),vmat(2,2),wvec(2),wmax
      double complex FC0(-2:0),FC1(y1max,-2:0),FC2(y2max,-2:0),
     & FC3(y3max,-2:0),
     & FB0_1(-2:0),FB1_1(y1max,-2:0),FB2_1(y2max,-2:0),
     & FB0_2(-2:0),FB1_2(y1max,-2:0),FB2_2(y2max,-2:0),
     & FB0_3(-2:0),FB1_3(y1max,-2:0),FB2_3(y2max,-2:0),
     & C00(-2:0),B00_1(-2:0),B00_2(-2:0),B00_3(-2:0),trI3,
     & tau3(4,-2:0),RHS(2),inRHS(2)
      integer n1,n2,n3,ep,indx(2)
      logical failed,iterate,dosvd
      double precision para(Pcc)
      integer jtable,j,Ntrue
      logical,save:: first=.true.
      double precision,save:: tableC(Pcc,Ncmax)      
      integer,save :: Nstore=0
!$omp threadprivate(first,tableC,Nstore)

c--- controls whether or not to iterate the solution by LU decomposition
      iterate=.false.
c--- controls whether to use SVD instead of LU decomposition
c--- (default behaviour is to use LU, then resort to SVD if the
c---  resulting tensor fails the consistency check)
      dosvd=.false.

      if (clear(3)) then
      clear(3)=.false.
      Nstore=0
      endif

      if (Nstore .gt. Ncmax) then
      print * 
      print *, 'ovCtensor: Nstore .gt. Ncmax'
      print *, 'Nstore,Ncmax',Nstore,Ncmax
      print *, 'Either adjust Ncmax in Cnames.f and recompile'
      print *, 'or call clearcache to clear the cache.'
      stop
      endif
      do j=1,4
      para(j)=p1(j)
      para(4+j)=p2(j)
      enddo
      para(9)=xm1s
      para(10)=xm2s
      para(11)=xm3s
C if parameter set is found set pvBcache equal to the starting
C value
      if (Nstore .eq. 0) go to 20
      do jtable=1,Nstore
      Ntrue=0
        do j=1,Pcc
        if (abs(para(j)-tableC(j,jtable)) .lt. 1d-8) then
          Ntrue=Ntrue+1
        else
          exit
        endif 
        enddo
        if (Ntrue .eq. Pcc) then
c--- retrieve from cache
c          write(6,*) 'Retrieving from cache: ',jtable
          do ep=-2,0
            FC0(ep)=FC0save(jtable,ep)
            C00(ep)=C00save(jtable,ep)
            do j=1,y1max
              FC1(j,ep)=FC1save(jtable,j,ep)
            enddo
            do j=1,y2max
              FC2(j,ep)=FC2save(jtable,j,ep)
            enddo
            do j=1,y3max
              FC3(j,ep)=FC3save(jtable,j,ep)
            enddo
            do j=1,4
              tau3(j,ep)=tau3save(jtable,j,ep)
            enddo
          enddo
          return
        endif
      enddo

C    if parameter set is not found we have to calculate
 20   continue
      Nstore=Nstore+1
      do j=1,Pcc
      tableC(j,Nstore)=para(j)
      enddo
c      write(6,*) 'Computing new Nstore: ',Nstore

      if (first) then
      first=.false.
      call ovarraysetup
      endif
      
      p12(:)=p1(:)+p2(:)
      p3(:)=-p12(:)
      p1Dp1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      p2Dp2=p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2
      p3Dp3=p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2
      p1Dp2=p1(4)*p2(4)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      
      if (abs(p1Dp1) .lt. onshellcutoff) p1Dp1=0d0     
      if (abs(p2Dp2) .lt. onshellcutoff) p2Dp2=0d0     
      if (abs(p3Dp3) .lt. onshellcutoff) p3Dp3=0d0  
      if (abs(xm1s) .lt. onshellcutoff) then
        m1s=0d0
      else
        m1s=xm1s
      endif
      if (abs(xm2s) .lt. onshellcutoff) then
        m2s=0d0
      else
        m2s=xm2s
      endif
      if (abs(xm3s) .lt. onshellcutoff) then
        m3s=0d0
      else
        m3s=xm3s
      endif
      
      inGram(1,1)=p1Dp1
      inGram(2,1)=p1Dp2
      inGram(1,2)=p1Dp2
      inGram(2,2)=p2Dp2
     
C----calculate bubble tensor integrals 
      call ovBtensor(p2,m2s,m3s,FB0_1,FB1_1,FB2_1,B00_1)
      call ovBtensor(p12,m1s,m3s,FB0_2,FB1_2,FB2_2,B00_2)
      call ovBtensor(p1,m1s,m2s,FB0_3,FB1_3,FB2_3,B00_3)

      s1Dp1=0.5d0*(m2s-m1s-p1Dp1)
      s2Dp2=0.5d0*(m3s-m2s-p2Dp2)
      s1Dp2=s2Dp2-p1Dp2

c--- point to restart from, if necessary
   66 continue

      if (dosvd) then
        call ovdsvdcmp(inGram,Gram,2,2,wvec,vmat)
        Gram2=-wvec(1)*wvec(2)
        wmax=max(wvec(1),wvec(2))
        do n1=1,2
        if (wvec(n1) .lt. wmax*1d-6) wvec(n1)=0d0
        enddo
      else
        call ludcmp(inGram,2,indx,d,Gram)   
        Gram2=d*Gram(1,1)*Gram(2,2)
      endif
      
      do ep=-2,0
      FC0(ep)=trI3(p1Dp1,p2Dp2,p3Dp3,m1s,m2s,m3s,musq,ep)     
      inRHS(1)=s1Dp1*FC0(ep)+half*FB0_2(ep)-half*FB0_1(ep)
      inRHS(2)=s1Dp2*FC0(ep)+half*FB0_3(ep)-half*FB0_2(ep)
      if (dosvd) then
        call zsvbksb(Gram,wvec,vmat,2,2,inRHS,RHS)
      else
        call zlubksb(Gram,2,indx,inRHS,RHS)
        if (iterate) call mprove(inGram,Gram,2,indx,inRHS,RHS)
      endif
      
      do n1=1,4
      FC1(n1,ep)=p1(n1)*RHS(1)+p2(n1)*RHS(2)
      enddo 
      enddo

      if (maxcindex .eq. 1) goto 99
      
      do ep=-2,0
      C00(ep)=m1s*FC0(ep)+FB0_1(ep)
      do n2=1,4
      inRHS(1)=s1Dp1*FC1(n2,ep)
     & +half*FB1_2(n2,ep)-half*(FB1_1(n2,ep)-p1(n2)*FB0_1(ep))
      inRHS(2)=s1Dp2*FC1(n2,ep)+half*FB1_3(n2,ep)-half*FB1_2(n2,ep)
      if (dosvd) then
        call zsvbksb(Gram,wvec,vmat,2,2,inRHS,RHS)
      else
        call zlubksb(Gram,2,indx,inRHS,RHS)
        if (iterate) call mprove(inGram,Gram,2,indx,inRHS,RHS)
      endif
      do n1=n2,4
!     Setup temporary value
      FC2(y2(n1,n2),ep)=p1(n1)*RHS(1)+p2(n1)*RHS(2)
      if(n1.eq.n2) then
        if (n1 .lt. 4) then
          C00(ep)=C00(ep)+FC2(y2(n2,n2),ep)
        else
          C00(ep)=C00(ep)-FC2(y2(n2,n2),ep)
        endif
      endif
      enddo 
      enddo 
      enddo 

c--- to account for factor of 1/(n-2)
      C00(0)=0.5d0*(C00(0)+C00(-1))
      C00(-1)=0.5d0*C00(-1)

      do ep=-2,0
      do n2=1,4
      do n1=n2,4
      FC2(y2(n1,n2),ep)=FC2(y2(n1,n2),ep)
     & +ovw2(n1,n2,p1,p2,Gram2)*C00(ep)
      enddo 
      enddo 
      enddo 

      if (maxcindex .eq. 2) goto 99
      
      do ep=-2,0
      inRHS(1)=s1Dp1*C00(ep)+half*(B00_2(ep)-B00_1(ep))
      inRHS(2)=s1Dp2*C00(ep)+half*(B00_3(ep)-B00_2(ep))
      if (dosvd) then
        call zsvbksb(Gram,wvec,vmat,2,2,inRHS,RHS)
      else
        call zlubksb(Gram,2,indx,inRHS,RHS)
        if (iterate) call mprove(inGram,Gram,2,indx,inRHS,RHS)
      endif
      do n1=1,4
      tau3(n1,ep)=RHS(1)*p1(n1)+RHS(2)*p2(n1)
      enddo
      enddo 

      do ep=-2,0
      do n3=1,4
      do n2=n3,4
      inRHS(1)=s1Dp1*FC2(y2(n2,n3),ep)
     & +half*(FB2_2(y2(n2,n3),ep)
     & -(FB2_1(y2(n2,n3),ep)-p1(n2)*FB1_1(n3,ep)
     &  -p1(n3)*FB1_1(n2,ep)+p1(n2)*p1(n3)*FB0_1(ep)))
      inRHS(2)=s1Dp2*FC2(y2(n2,n3),ep)
     & +half*(FB2_3(y2(n2,n3),ep)-FB2_2(y2(n2,n3),ep))
      if (dosvd) then
        call zsvbksb(Gram,wvec,vmat,2,2,inRHS,RHS)
      else
        call zlubksb(Gram,2,indx,inRHS,RHS)
        if (iterate) call mprove(inGram,Gram,2,indx,inRHS,RHS)
      endif
      do n1=n2,4
      FC3(y3(n1,n2,n3),ep)=
     & p1(n1)*RHS(1)+p2(n1)*RHS(2)
     & +ovw2(n1,n3,p1,p2,Gram2)*tau3(n2,ep)
     & +ovw2(n1,n2,p1,p2,Gram2)*tau3(n3,ep)
      enddo 
      enddo 
      enddo 
      enddo 

c--- before returning, check tensor computed correctly        
   99 continue
      
       call ovCcheck(maxcindex,p1,p12,m1s,m2s,m3s,
     &  FC0,FC1,FC2,FC3,failed)
      if (failed) then
c        write(6,*) 'badpoint set in ovCtensor'
        pvbadpoint=.true.
      endif

c--- if we found a bad point and haven't used svd yet, try it      
      if ((pvbadpoint) .and. (dosvd .eqv. .false.)) then
        pvbadpoint=.false.
        dosvd=.true.
        goto 66
      endif
      
c--- store in cache
      do ep=-2,0
        FC0save(Nstore,ep)=FC0(ep)
        C00save(Nstore,ep)=C00(ep)
        do j=1,y1max
          FC1save(Nstore,j,ep)=FC1(j,ep)
        enddo
        do j=1,y2max
          FC2save(Nstore,j,ep)=FC2(j,ep)
        enddo
        do j=1,y3max
          FC3save(Nstore,j,ep)=FC3(j,ep)
        enddo
        do j=1,4
          tau3save(Nstore,j,ep)=tau3(j,ep)
        enddo
      enddo
      
      return
      end

