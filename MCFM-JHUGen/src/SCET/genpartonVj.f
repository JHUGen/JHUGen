      subroutine genpartonVj(p,r1,r2,r3,r4,tcut,ph,wt,pass)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      integer branch,clust
      logical pass
      real(dp):: r1,r2,r3,r4,wt,tcut
      real(dp):: Q(4),pah(4),pbh(4),Qh(4),p1(4),p2(4),pj(4)
      real(dp):: p(mxpart,4),ph(mxpart,4),p1h(4),p2h(4)
!
!     Radiate an additional massless particle (from p_branch).
!     using p1+p2 --> Q+pj to generate ph1+ph2 --> Qh+ph1+ph2 

      p1(:)=p(1,:)
      p2(:)=p(2,:)
      Q(:) =p(3,:)+p(4,:)
      pj(:)=p(5,:)
      branch=1+int(r4*4._dp)
      if (branch.le.2) then
         call gen_initVj(branch,p1,p2,Q,pj,r1,r2,r3,tcut,
     &                   pah,pbh,Qh,p1h,p2h,wt,pass)
         call cluster(p1h,p2h,clust)
c         if (clust.ne.1) pass=.false.
         if (clust.eq.2) pass=.false.
      else
         branch=branch-2
         call gen_final(branch,p1,p2,Q,pj,r1,r2,r3,tcut,
     &                  pah,pbh,Qh,p1h,p2h,wt,pass) 
         call cluster(p1h,p2h,clust)
         if (clust.ne.2) pass=.false.
c         if (clust.eq.1) pass=.false.
      endif
      wt=wt*2._dp

      ph(1,:)=pah(:)
      ph(2,:)=pbh(:)
      ph(3,:)=Qh(:)
      ph(4,:)=p1h(:)
      ph(5,:)=p2h(:)
      
      return
      end
      

! ------------------------------------------------------------------------------
! ---------------------   INITIAL STATE BRANCHING ------------------------------
! ------------------------------------------------------------------------------
 
      subroutine gen_initVj(b,p1,p2,Q,pj,r1,r2,r3,t0,
     & pah,pbh,Qh,p1h,p2h,wt,pass)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'energy.f'
      include 'mxpart.f'
      logical pass
      integer b
      real(dp):: wt,t0,xa,xb,ta,tb,Jac,r1,r2,r3,tmax,p2tDQ
      real(dp):: phi,pt,al,xah,xbh,p2DQ,pl2
      real(dp):: Q(4),pab(4),Qh(4),pah(4),pbh(4),p1h(4),p2h(4)
      real(dp):: p1(4),p2(4),pj(4)

      pass=.false.
      pab=p1+p2
      xa=(pab(4)+pab(3))/sqrts
      xb=(pab(4)-pab(3))/sqrts

!     pick tau_a and tau_b of emitted particle
      tmax=sqrts
      if (tmax.lt.t0) return
      if (b.eq.1) then
         call pick(2,tb,t0,tmax,r1,wt)
         call pick(2,ta,t0,tmax,r2,wt)
!         call pick_tau(t0,t0,xb,xa,r1,tb,wt)
!         if (tb.le.sqrts*(1-xb)*(1-xa*xb)/(2-xa-xb)) then
!            call pick_tau(t0,tb,xa,xb,r2,ta,wt)
!         else
!            call pick_tau1(t0,tb,xa,xb,r2,ta,wt)
!         endif
      endif
      if (b.eq.2) then
         call pick(2,ta,t0,tmax,r1,wt)
         call pick(2,tb,t0,tmax,r2,wt)
!         call pick_tau(t0,t0,xa,xb,r1,ta,wt)
!         if (ta.le.sqrts*(1-xa)*(1-xa*xb)/(2-xa-xb)) then
!            call pick_tau(t0,ta,xb,xa,r2,tb,wt)
!         else
!            call pick_tau1(t0,ta,xb,xa,r2,tb,wt)
!         endif
      endif
      wt=0.25_dp*wt   
      call pick(1,phi,0._dp,2._dp*pi,r3,wt)
!     construct p2h
      pt=sqrt(abs(ta*tb))
      p2h(4)=0.5_dp*(ta+tb)
      p2h(3)=0.5_dp*(ta-tb)
      p2h(2)=pt*sin(phi)
      p2h(1)=pt*cos(phi)
      p1h=pj
!     jet momentum is left unaltered
!     construct intial state particles
      p2DQ=p2h(4)*Q(4)-p2h(3)*Q(3)
      p2tDQ=p2h(1)*Q(1)+p2h(2)*Q(2)
      pl2=pt**2
      Jac=pl2**2-two*pl2*p2tDQ+p2DQ**2
      al=(p2DQ-pl2-sqrt(abs(Jac)))/pl2
      Pab(4)=Pab(4)-al*p2h(4)
      Pab(3)=Pab(3)-al*p2h(3)
      Pab(2)=zip
      Pab(1)=zip
      xah=(Pab(4)+Pab(3))/sqrts
      xbh=(Pab(4)-Pab(3))/sqrts
      if ((xah.gt.1d0).or.(xbh.gt.1d0)) return
      pah(4)=half*sqrts*xah
      pah(3)=half*sqrts*xah
      pah(2)=zip
      pah(1)=zip
      pbh(4)=half*sqrts*xbh
      pbh(3)=-half*sqrts*xbh
      pbh(2)=zip
      pbh(1)=zip
!     construct Qhat
      Qh(4)=Q(4)-(al+1d0)*p2h(4)
      Qh(3)=Q(3)-(al+1d0)*p2h(3)
      Qh(2)=Q(2)-p2h(2)
      Qh(1)=Q(1)-p2h(1)
!     calculate Jacobian
      p2DQ=p2h(4)*Qh(4)-p2h(3)*Qh(3)
      p2tDQ=p2h(1)*Qh(1)+p2h(2)*Qh(2)
      Jac=sqrt(abs((pl2**2+2d0*pl2*p2tDQ+p2DQ**2)/Jac))

      wt=wt*Jac/twopi**3
!     Change flux factor
!      wt=wt*(xa/xah)*(xb/xbh)
      pass=.true.

!      write(*,*) '----------------------------------------------------'
!      write(*,*) 'p1',p1,dotpr(p1,p1)
!      write(*,*) 'p2',p2,dotpr(p2,p2)
!      write(*,*) 'Q ',Q,sqrt(dotpr(Q,Q))
!      write(*,*) 'pj',pj,dotpr(pj,pj)
!      write(*,*) '++',p1+p2-Q-pj
!      write(*,*)
!      write(*,*)
!      write(*,*) 'pah',pah,dotpr(pah,pah)
!      write(*,*) 'pbh',pbh,dotpr(pbh,pbh)
!      write(*,*) 'Qh ',Qh,sqrt(dotpr(Qh,Qh))
!      write(*,*) 'p1h',p1h,dotpr(p1h,p1h)
!      write(*,*) 'p2h',p2h,dotpr(p2h,p2h)
!      write(*,*) '++',pah+pbh-Qh-p1h-p2h
!      write(*,*)
!      write(*,*)
      return
      end


! ------------------------------------------------------------------------------
! ----------------------   FINAL STATE BRANCHING -------------------------------
! ------------------------------------------------------------------------------

      subroutine gen_final(b,p1,p2,Q,pj,r1,r2,r3,t0,
     & pah,pbh,Qh,p1h,p2h,wt,pass)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'energy.f'
      logical pass
      integer b
      real(dp):: t0,t1,t2,tmax,wt,dotpr,phi,pt2,Jac
      real(dp):: paDpb,pjDpa,pjDpb,pjbDn1,pjbDpa,pjbDpj,pjDpt
      real(dp):: pjDpl,pl2,exp1,sqarg,a,bm,bp,s12,xa,xb,xah,xbh
      real(dp):: r1,r2,r3,pa(4),pb(4),Q(4),p1(4),p2(4),pjb(4)
      real(dp):: pah(4),pbh(4),Qh(4),p1h(4),p2h(4),n1(4),n2(4)
      real(dp):: pt(4),pj(4),pl(4),p12(4),pab(4)

      pass=.false.
      Qh(:)=Q(:)
      xa=two*p1(4)/sqrts
      xb=two*p2(4)/sqrts
!      write(*,*) xa,xb
      if (b.eq.1) then
         pa(:)=p1(:)
         pb(:)=p2(:)
      else
         pa(:)=p2(:)
         pb(:)=p1(:)
      endif
      
* caculate kinematics and construct momenta
      tmax=sqrts
      if (tmax.lt.t0) return
      call pick(2,t2,t0,tmax,r1,wt)
      call pick(2,t1,t0,tmax,r2,wt)
      call pick(1,phi,0._dp,2._dp*pi,r3,wt)
      wt=0.5_dp*wt
c      write(6,*) 't1,t2',t1,t2

      pjDpa=dotpr(pj,pa)
      pjDpb=dotpr(pj,pb)
      paDpb=dotpr(pa,pb)
      n1(:)=pjDpb*pa(:)-pjDpa*pb(:)+paDpb*pj(:)
      n1(:)=n1(:)/sqrt(abs(dotpr(n1,n1)))

      pjb(1)=pj(2)
      pjb(2)=pj(3)
      pjb(3)=pj(1)
      pjb(4)=pj(4)
      pjbDn1=dotpr(pjb,n1)
      pjbDpj=dotpr(pjb,pj)
      pjbDpa=dotpr(pjb,pa)
      n2(:)=-pjDpa*pjbDn1*n1(:)+pjbDpj*pa(:)+pjbDpa*pj(:)-pjDpa*pjb(:)
      n2(:)=n2(:)/sqrt(abs(dotpr(n2,n2)))
      pt2=2._dp*t1*t2*pa(4)*pj(4)/pjDpa
      pt(:)=sqrt(pt2)*(cos(phi)*n1(:)+sin(phi)*n2(:))
      p1h(:)=(t1*pa(4)*pj(:)+t2*pj(4)*pa(:))/pjDpa+pt(:)

      pt(1)=p1h(1)
      pt(2)=p1h(2)
      pt(3)=0._dp
      pt(4)=0._dp
      pl(1)=0._dp
      pl(2)=0._dp
      pl(3)=pj(3)
      pl(4)=pj(4)
      pl2=dotpr(pl,pl)
      pjDpl=dotpr(pl,p1h)
      pjDpt=-dotpr(pt,pj)
      exp1=pl2-2._dp*pjDpt
      sqarg=pjDpl**2+exp1*pl2
         
      bm=(pjDpl-sqrt(sqarg))/pl2
      bp=(pjDpl+sqrt(sqarg))/pl2
      p2h(1)=pj(1)-p1h(1)
      p2h(2)=pj(2)-p1h(2)
      p2h(3)=bp*pj(3)-p1h(3)
      p2h(4)=bp*pj(4)-p1h(4)
      pab(:)=pa(:)+pb(:)-(1._dp-bp)*pl(:)
      xah=(Pab(4)+Pab(3))/sqrts
      xbh=(Pab(4)-Pab(3))/sqrts
      if ((xah.ge.1._dp).or.(xbh.ge.1._dp)) return
      pah(4)=0.5_dp*sqrts*xah
      pah(3)=0.5_dp*sqrts*xah
      pah(2)=0._dp
      pah(1)=0._dp
      pbh(4)=0.5_dp*sqrts*xbh
      pbh(3)=-0.5_dp*sqrts*xbh
      pbh(2)=0._dp
      pbh(1)=0._dp


* calculate the event weight
      Jac=1._dp/(pl2*abs(bp-bm))
      pl(1)=0._dp
      pl(2)=0._dp
      pl(3)=p1h(3)+p2h(3)
      pl(4)=p1h(4)+p2h(4)
      pl2=dotpr(pl,pl)
      p12(:)=p1h(:)+p2h(:)
      s12=dotpr(p12,p12)
      a=2._dp*sqrt(pl2*(pl2-s12))
      wt=wt*a*Jac/twopi**3
* construct initial state momenta
!     Change flux factor
!      wt=wt*(xa/xah)*(xb/xbh)
      pass=.true.
      return
      end



      subroutine cluster(p1,p2,clust)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'jetcuts.f'
      integer clust,accjets
      real(dp):: Deta,Dphi,R
      real(dp):: pt1,pt2,eta1,eta2,Rcut
      real(dp):: p1(4),p2(4)
      common/Rcut/Rcut

      pt1=sqrt(p1(1)**2+p1(2)**2)
      pt2=sqrt(p2(1)**2+p2(2)**2)
      eta1=half*log((p1(4)-p1(3))/(p1(4)+p1(3)))
      eta2=half*log((p2(4)-p2(3))/(p2(4)+p2(3)))
      Deta=eta1-eta2
      Dphi=acos((p1(1)*p2(1)+p1(2)*p2(2))/pt1/pt2)
      R=sqrt(Deta**2+Dphi**2)
      clust=1

!      write(6,*) 'ptjetmin,Rcut',ptjetmin,Rcut
      
      if (R.lt.Rcut) then
         clust=2
      else
         accjets=0
         if ((pt1>ptjetmin).and.(abs(eta1)<etajetmax))
     &     accjets=accjets+1
         if ((pt2>ptjetmin).and.(abs(eta2)<etajetmax))
     &     accjets=accjets+1
         if (accjets == 2) clust=3
         if (accjets == 0) clust=0
c         if ((pt1.gt.ptjetmin).and.(pt2.gt.ptjetmin)) clust=3
c         if ((pt1.lt.ptjetmin).and.(pt2.lt.ptjetmin)) clust=0
      endif
      
      return
      end
