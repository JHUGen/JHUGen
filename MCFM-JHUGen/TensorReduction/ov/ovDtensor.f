      subroutine ovDtensor(p1,p2,p3,m1s,m2s,m3s,m4s,
     & FD0,FD1,FD2,FD3,FD4)
      implicit none
C     p1,p2,p3 are the external momenta,
C     m1s,m2s,m3s,m4s are the squares of the internal masses
C     FD0...FD4 are the rank 0,...4 box functions
C     Lorentz indices are stored as linear array, thus FD2(y2(n1,n2),ep)
C     Author: R.K.Ellis (January 2013)
      include 'TRconstants.f'
      include 'TRydef.f'
      include 'TRscale.f'
      include 'TRbadpoint.f'
      include 'TRmaxindex.f'
      include 'TRclear.f'
      include 'ovDnames.f'
      include 'ovDsave.f'
      double precision p1(4),p2(4),p3(4),p4(4),m1s,m2s,m3s,m4s,
     & p12(4),p23(4),p123(4),
     & p1Dp1,p2Dp2,p3Dp3,p4Dp4,p1Dp2,p1Dp3,p2Dp3,
     & s1Ds1,s1Dp1,s2Dp2,s3Dp3,s2Dp3,s1Dp2,s1Dp3,
     & s12,s23,ovw3,w(4,4),inGram(3,3),Gram(3,3),d,Gram3,
     & vmat(3,3),wvec(3),wmax
      double complex FD0(-2:0),FD1(y1max,-2:0),FD2(y2max,-2:0),
     &FD3(y3max,-2:0),FD4(y4max,-2:0),
     &FC0_1(-2:0),FC1_1(y1max,-2:0),FC2_1(y2max,-2:0),FC3_1(y3max,-2:0),
     &FC0_2(-2:0),FC1_2(y1max,-2:0),FC2_2(y2max,-2:0),FC3_2(y3max,-2:0),
     &FC0_3(-2:0),FC1_3(y1max,-2:0),FC2_3(y2max,-2:0),FC3_3(y3max,-2:0),
     &FC0_4(-2:0),FC1_4(y1max,-2:0),FC2_4(y2max,-2:0),FC3_4(y3max,-2:0),
     &D00(-2:0),C00_1(-2:0),C00_2(-2:0),C00_3(-2:0),C00_4(-2:0),
     &D0000(-2:0),trI4,tau4(4,-2:0),
     &tau3_1(4,-2:0),tau3_2(4,-2:0),tau3_3(4,-2:0),tau3_4(4,-2:0),
     &tmp(4,-2:0),RHS(3),inRHS(3),tmpRHS(3,y3max)
      integer n1,n2,n3,n4,ep,indx(3)
      logical failed,iterate,dosvd
      double precision para(Pdd)
      logical,save:: first=.true.
      double precision,save::tableD(Pdd,Ndmax)      
      integer jtable,j,Ntrue
      integer, save:: Nstore=0
!$omp threadprivate(first,Nstore,tableD)

c--- controls whether or not to iterate the solution by LU decomposition
      iterate=.false.
c--- controls whether to use SVD instead of LU decomposition
c--- (default behaviour is to use LU, then resort to SVD if the
c---  resulting tensor fails the consistency check)
      dosvd=.false.

      if (clear(4)) then
      clear(4)=.false.
      Nstore=0
      endif

      if (Nstore .gt. Ndmax) then
      print * 
      print *, 'ovDtensor: Nstore .gt. Ndmax'
      print *, 'Nstore,Ndmax',Nstore,Ndmax
      print *, 'Either adjust Ndmax in Dnames.f and recompile'
      print *, 'or call clearcache to clear the cache.'
      stop
      endif
      do j=1,4
      para(j)=p1(j)
      para(4+j)=p2(j)
      para(8+j)=p3(j)
      enddo
      para(13)=m1s
      para(14)=m2s
      para(15)=m3s
      para(16)=m4s
C if parameter set is found set pvBcache equal to the starting
C value
      if (Nstore .eq. 0) go to 20
      do jtable=1,Nstore
      Ntrue=0
        do j=1,Pdd
        if (abs(para(j)-tableD(j,jtable)) .lt. 1d-8) then
          Ntrue=Ntrue+1
        else
          exit
        endif 
        enddo
        if (Ntrue .eq. Pdd) then
c--- retrieve from cache
c          write(6,*) 'Retrieving from cache: ',jtable
          do ep=-2,0
            FD0(ep)=FD0save(jtable,ep)
            do j=1,y1max
              FD1(j,ep)=FD1save(jtable,j,ep)
            enddo
            do j=1,y2max
              FD2(j,ep)=FD2save(jtable,j,ep)
            enddo
            do j=1,y3max
              FD3(j,ep)=FD3save(jtable,j,ep)
            enddo
            do j=1,y4max
              FD4(j,ep)=FD4save(jtable,j,ep)
            enddo
          enddo
          return
        endif
      enddo

C    if parameter set is not found we have to calculate
 20   continue
      Nstore=Nstore+1
      do j=1,Pdd
      tableD(j,Nstore)=para(j)
      enddo
c      write(6,*) 'Computing new Nstore: ',Nstore

      if (first) then
      first=.false.
      call ovarraysetup
      endif
      
      p12(:)=p1(:)+p2(:)
      p23(:)=p2(:)+p3(:)
      p123(:)=p12(:)+p3(:)
      p4(:)=-p123(:)
      p1Dp1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      p2Dp2=p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2
      p3Dp3=p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2
      p4Dp4=p4(4)**2-p4(1)**2-p4(2)**2-p4(3)**2
      s12=p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2
      s23=p23(4)**2-p23(1)**2-p23(2)**2-p23(3)**2
      p1Dp2=0.5d0*(s12-p1Dp1-p2Dp2)
      p2Dp3=0.5d0*(s23-p2Dp2-p3Dp3)
      p1Dp3=p1(4)*p3(4)-p1(1)*p3(1)-p1(2)*p3(2)-p1(3)*p3(3)

      inGram(1,1)=p1Dp1
      inGram(2,2)=p2Dp2
      inGram(3,3)=p3Dp3
      inGram(1,2)=p1Dp2
      inGram(2,3)=p2Dp3
      inGram(1,3)=p1Dp3
      inGram(3,1)=inGram(1,3)
      inGram(2,1)=inGram(1,2)
      inGram(3,2)=inGram(2,3)
     
c--- point to restart from, if necessary
   66 continue

      if (dosvd) then
        call ovdsvdcmp(inGram,Gram,3,3,wvec,vmat)
        Gram3=wvec(1)*wvec(2)*wvec(3)
        wmax=max(wvec(1),wvec(2),wvec(3))
        do n1=1,3
        if (wvec(n1) .lt. wmax*1d-6) wvec(n1)=0d0
        enddo
      else
        call ludcmp(inGram,3,indx,d,Gram)   
        Gram3=d*Gram(1,1)*Gram(2,2)*Gram(3,3)
      endif
      
C----calculate triangle tensor integrals
      call ovCtensor(p2,p3,m2s,m3s,m4s,
     & FC0_1,FC1_1,FC2_1,FC3_1,C00_1,tau3_1)
      call ovCtensor(p12,p3,m1s,m3s,m4s,
     & FC0_2,FC1_2,FC2_2,FC3_2,C00_2,tau3_2)
      call ovCtensor(p1,p23,m1s,m2s,m4s,
     & FC0_3,FC1_3,FC2_3,FC3_3,C00_3,tau3_3)
      call ovCtensor(p1,p2,m1s,m2s,m3s,
     & FC0_4,FC1_4,FC2_4,FC3_4,C00_4,tau3_4)
     
      s1Ds1=m1s
      s1Dp1=0.5d0*(m2s-m1s-p1Dp1)
      s2Dp2=0.5d0*(m3s-m2s-p2Dp2)
      s3Dp3=0.5d0*(m4s-m3s-p3Dp3)
      s2Dp3=s3Dp3-p2Dp3
      s1Dp2=s2Dp2-p1Dp2
      s1Dp3=s2Dp3-p1Dp3

      do ep=-2,0
      FD0(ep)=trI4(p1Dp1,p2Dp2,p3Dp3,p4Dp4,s12,s23,
     & m1s,m2s,m3s,m4s,musq,ep)     
      enddo

      do ep=-2,0
      inRHS(1)=s1Dp1*FD0(ep)+half*(FC0_2(ep)-FC0_1(ep))
      inRHS(2)=s1Dp2*FD0(ep)+half*(FC0_3(ep)-FC0_2(ep))
      inRHS(3)=s1Dp3*FD0(ep)+half*(FC0_4(ep)-FC0_3(ep))
      if (dosvd) then
        call zsvbksb(Gram,wvec,vmat,3,3,inRHS,RHS)
      else
        call zlubksb(Gram,3,indx,inRHS,RHS)
        if (iterate) call mprove(inGram,Gram,3,indx,inRHS,RHS)
      endif
      do n1=1,4
      FD1(n1,ep)=p1(n1)*RHS(1)+p2(n1)*RHS(2)+p3(n1)*RHS(3)
      enddo 
      enddo 

      if (maxdindex .eq. 1) goto 99
      
c--- Start of rank 2      
      do n1=1,4
      do n2=1,4
      w(n1,n2)=ovw3(n1,n2,p1,p2,p3,Gram3)
      enddo
      enddo
      
      do ep=-2,0
      do n2=1,4
      inRHS(1)=s1Dp1*FD1(n2,ep)
     & +half*(FC1_2(n2,ep)-(FC1_1(n2,ep)-p1(n2)*FC0_1(ep)))
      inRHS(2)=s1Dp2*FD1(n2,ep)+half*(FC1_3(n2,ep)-FC1_2(n2,ep))
      inRHS(3)=s1Dp3*FD1(n2,ep)+half*(FC1_4(n2,ep)-FC1_3(n2,ep))
      if (dosvd) then
        call zsvbksb(Gram,wvec,vmat,3,3,inRHS,RHS)
      else
        call zlubksb(Gram,3,indx,inRHS,RHS)
        if (iterate) call mprove(inGram,Gram,3,indx,inRHS,RHS)
      endif
      do n1=n2,4
      FD2(y2(n1,n2),ep)=p1(n1)*RHS(1)+p2(n1)*RHS(2)+p3(n1)*RHS(3)
      enddo
      enddo
      enddo
      
      D00(-2:-1)=czip
      ep=0
      D00(ep)=s1Ds1*FD0(ep)+FC0_1(ep)
      do n1=1,4
      if (n1 .lt. 4) then
        D00(ep)=D00(ep)+FD2(y2(n1,n1),ep)
      else
        D00(ep)=D00(ep)-FD2(y2(n1,n1),ep)
      endif
      enddo 

c--- Only need ep=0 because D00 finite
      ep=0
      do n2=1,4
      do n1=n2,4
      FD2(y2(n1,n2),ep)=FD2(y2(n1,n2),ep)+w(n1,n2)*D00(ep)
      enddo 
      enddo 

      if (maxdindex .eq. 2) goto 99
      
c--- Start of rank 3      
      do ep=-2,0
      inRHS(1)=s1Dp1*D00(ep)+half*(C00_2(ep)-C00_1(ep))
      inRHS(2)=s1Dp2*D00(ep)+half*(C00_3(ep)-C00_2(ep))
      inRHS(3)=s1Dp3*D00(ep)+half*(C00_4(ep)-C00_3(ep))
      if (dosvd) then
        call zsvbksb(Gram,wvec,vmat,3,3,inRHS,RHS)
      else
        call zlubksb(Gram,3,indx,inRHS,RHS)
        if (iterate) call mprove(inGram,Gram,3,indx,inRHS,RHS)
      endif
      do n1=1,4
      tau4(n1,ep)=RHS(1)*p1(n1)+RHS(2)*p2(n1)+RHS(3)*p3(n1)
      enddo
      enddo 

      do ep=-2,0
      do n3=1,4
      do n2=n3,4
      inRHS(1)=s1Dp1*FD2(y2(n2,n3),ep)
     &+half*(FC2_2(y2(n2,n3),ep)-(FC2_1(y2(n2,n3),ep)
     &-p1(n2)*FC1_1(n3,ep)-p1(n3)*FC1_1(n2,ep)+p1(n2)*p1(n3)*FC0_1(ep)))
      inRHS(2)=s1Dp2*FD2(y2(n2,n3),ep)
     &+half*(FC2_3(y2(n2,n3),ep)-FC2_2(y2(n2,n3),ep))
      inRHS(3)=s1Dp3*FD2(y2(n2,n3),ep)
     &+half*(FC2_4(y2(n2,n3),ep)-FC2_3(y2(n2,n3),ep))
      if (dosvd) then
        call zsvbksb(Gram,wvec,vmat,3,3,inRHS,RHS)
      else
        call zlubksb(Gram,3,indx,inRHS,RHS)
        if (iterate) call mprove(inGram,Gram,3,indx,inRHS,RHS)
      endif
      do n1=n2,4
      FD3(y3(n1,n2,n3),ep)=
     & p1(n1)*RHS(1)+p2(n1)*RHS(2)+p3(n1)*RHS(3)
     & +w(n1,n2)*tau4(n3,ep)
     & +w(n1,n3)*tau4(n2,ep)
      enddo 
      enddo 
      enddo 
      enddo 

      if (maxdindex .eq. 3) goto 99

c--- Start of rank 4     
      D0000=czip
      do ep=-2,0
      do n2=1,4
      inRHS(1)=s1Dp1*tau4(n2,ep)+half*(tau3_2(n2,ep)-tau3_1(n2,ep))
      inRHS(2)=s1Dp2*tau4(n2,ep)+half*(tau3_3(n2,ep)-tau3_2(n2,ep))
      inRHS(3)=s1Dp3*tau4(n2,ep)+half*(tau3_4(n2,ep)-tau3_3(n2,ep))
      if (dosvd) then
        call zsvbksb(Gram,wvec,vmat,3,3,inRHS,RHS)
      else
        call zlubksb(Gram,3,indx,inRHS,RHS)
        if (iterate) call mprove(inGram,Gram,3,indx,inRHS,RHS)
      endif
!     Setup temporary value
      tmp(n2,ep)=p1(n2)*RHS(1)+p2(n2)*RHS(2)+p3(n2)*RHS(3)
      if (n2 .lt. 4) then
        D0000(ep)=D0000(ep)-tmp(n2,ep)
      else
        D0000(ep)=D0000(ep)+tmp(n2,ep)
      endif
      enddo 
      D0000(ep)=-D0000(ep)+m1s*D00(ep)+half*C00_1(ep)
      enddo            
C   factor of 1/[n-1]
      D0000(0)=1d0/3d0*(D0000(0)+2d0/3d0*D0000(-1))
      D0000(-1)=1d0/3d0*D0000(-1)+2d0/3d0*D0000(-2)

      do ep=-2,0
      do n4=1,4
      do n3=n4,4
      do n2=n3,4
      inRHS(1)=s1Dp1*FD3(y3(n2,n3,n4),ep)
     &+half*(FC3_2(y3(n2,n3,n4),ep)-(FC3_1(y3(n2,n3,n4),ep)
     & -p1(n2)*FC2_1(y2(n3,n4),ep)
     & -p1(n3)*FC2_1(y2(n4,n2),ep)
     & -p1(n4)*FC2_1(y2(n2,n3),ep)
     & +p1(n2)*p1(n3)*FC1_1(n4,ep)
     & +p1(n3)*p1(n4)*FC1_1(n2,ep)
     & +p1(n4)*p1(n2)*FC1_1(n3,ep)
     & -p1(n2)*p1(n3)*p1(n4)*FC0_1(ep)))
      inRHS(2)=s1Dp2*FD3(y3(n2,n3,n4),ep)
     &+half*(FC3_3(y3(n2,n3,n4),ep)-FC3_2(y3(n2,n3,n4),ep))
      inRHS(3)=s1Dp3*FD3(y3(n2,n3,n4),ep)
     &+half*(FC3_4(y3(n2,n3,n4),ep)-FC3_3(y3(n2,n3,n4),ep))
      if (dosvd) then
        call zsvbksb(Gram,wvec,vmat,3,3,inRHS,RHS)
      else
        call zlubksb(Gram,3,indx,inRHS,RHS)
        if (iterate) call mprove(inGram,Gram,3,indx,inRHS,RHS)
      endif
      do n1=n2,4
      FD4(y4(n1,n2,n3,n4),ep)=
     & p1(n1)*RHS(1)+p2(n1)*RHS(2)+p3(n1)*RHS(3)
      enddo 
      enddo 
      enddo 
      enddo 
      enddo 

      do ep=-2,0
      do n4=1,4
      do n3=n4,4
      do n2=n3,4
      do n1=n2,4
      if (n2 .eq. n3) then   ! Setup on first pass only
      inRHS(1)=
     & +s1Dp1*(w(n1,n3)*tau4(n4,ep)+w(n1,n4)*tau4(n3,ep))
     & -half*(w(n1,n3)*(tau3_1(n4,ep)-p1(n4)*C00_1(ep))
     &       +w(n1,n4)*(tau3_1(n3,ep)-p1(n3)*C00_1(ep)))
     & +half*(w(n1,n3)*tau3_2(n4,ep)+w(n1,n4)*tau3_2(n3,ep))
      inRHS(2)=
     & +s1Dp2*(w(n1,n3)*tau4(n4,ep)+w(n1,n4)*tau4(n3,ep))
     & -half*(w(n1,n3)*tau3_2(n4,ep)+w(n1,n4)*tau3_2(n3,ep))
     & +half*(w(n1,n3)*tau3_3(n4,ep)+w(n1,n4)*tau3_3(n3,ep))
      inRHS(3)=
     & +s1Dp3*(w(n1,n3)*tau4(n4,ep)+w(n1,n4)*tau4(n3,ep))
     & -half*(w(n1,n3)*tau3_3(n4,ep)+w(n1,n4)*tau3_3(n3,ep))
     & +half*(w(n1,n3)*tau3_4(n4,ep)+w(n1,n4)*tau3_4(n3,ep))
      if (dosvd) then
        call zsvbksb(Gram,wvec,vmat,3,3,inRHS,RHS)
      else
        call zlubksb(Gram,3,indx,inRHS,RHS)
        if (iterate) call mprove(inGram,Gram,3,indx,inRHS,RHS)
      endif
      tmpRHS(:,y3(n1,n3,n4))=RHS(:)
      endif
      FD4(y4(n1,n2,n3,n4),ep)=FD4(y4(n1,n2,n3,n4),ep)
     & +p1(n2)*tmpRHS(1,y3(n1,n3,n4))
     & +p2(n2)*tmpRHS(2,y3(n1,n3,n4))
     & +p3(n2)*tmpRHS(3,y3(n1,n3,n4))
      enddo 
      enddo 
      enddo 
      enddo 
      enddo 

      do ep=-2,0
      do n4=1,4
      do n3=n4,4
      do n2=n3,4
      do n1=n2,4
      if (n3 .eq. n4) then   ! Setup on first pass only
      inRHS(1)=
     & +s1Dp1*w(n1,n2)*tau4(n4,ep)
     & -half*w(n1,n2)*(tau3_1(n4,ep)-p1(n4)*C00_1(ep))
     & +half*w(n1,n2)*tau3_2(n4,ep)
      inRHS(2)=
     & +s1Dp2*w(n1,n2)*tau4(n4,ep)
     & -half*w(n1,n2)*tau3_2(n4,ep)
     & +half*w(n1,n2)*tau3_3(n4,ep)
      inRHS(3)=
     & +s1Dp3*w(n1,n2)*tau4(n4,ep)
     & -half*w(n1,n2)*tau3_3(n4,ep)
     & +half*w(n1,n2)*tau3_4(n4,ep)
      if (dosvd) then
        call zsvbksb(Gram,wvec,vmat,3,3,inRHS,RHS)
      else
        call zlubksb(Gram,3,indx,inRHS,RHS)
        if (iterate) call mprove(inGram,Gram,3,indx,inRHS,RHS)
      endif
      tmpRHS(:,y3(n1,n2,n4))=RHS(:)
      endif
      FD4(y4(n1,n2,n3,n4),ep)=FD4(y4(n1,n2,n3,n4),ep)
     & +p1(n3)*tmpRHS(1,y3(n1,n2,n4))
     & +p2(n3)*tmpRHS(2,y3(n1,n2,n4))
     & +p3(n3)*tmpRHS(3,y3(n1,n2,n4))
      enddo 
      enddo 
      enddo 
      enddo 
      enddo 

      do ep=-2,0
      do n4=1,4
      do n3=n4,4
      do n2=n3,4
      do n1=n2,4
      FD4(y4(n1,n2,n3,n4),ep)=FD4(y4(n1,n2,n3,n4),ep)
     &+(w(n1,n2)*w(n3,n4)+w(n1,n3)*w(n2,n4)+w(n1,n4)*w(n2,n3))*D0000(ep)
      enddo 
      enddo 
      enddo 
      enddo 
      enddo 
            
      
c--- before returning, check tensor computed correctly        
   99 continue
         
      call ovDcheck(maxdindex,p1,p12,p123,m1s,m2s,m3s,m4s,
     &  FD0,FD1,FD2,FD3,FD4,failed)
      if (failed) then
c        write(6,*) 'badpoint set in ovDtensor'
        pvbadpoint=.true.
      endif

c--- if we found a bad point and haven't used svd yet, try it      
      if ((pvbadpoint) .and. (dosvd .eqv. .false.)) then
        pvbadpoint=.false.
        dosvd=.true.
        goto 66
      endif
      
c      if (pvbadpoint) pause 'end of ovDtensor'

c--- store in cache
      do ep=-2,0
        FD0save(Nstore,ep)=FD0(ep)
        do j=1,y1max
          FD1save(Nstore,j,ep)=FD1(j,ep)
        enddo
        do j=1,y2max
          FD2save(Nstore,j,ep)=FD2(j,ep)
        enddo
        do j=1,y3max
          FD3save(Nstore,j,ep)=FD3(j,ep)
        enddo
        do j=1,y4max
          FD4save(Nstore,j,ep)=FD4(j,ep)
        enddo
      enddo
      
      return
      end

