      subroutine genQ(mQ,r1,tcut,p,wt)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'energy.f'
      include 'mxpart.f'
      real(dp):: r1,wt,S,Q2,tau,y,ymin,ymax,xa,xb,mQ,tcut
      real(dp):: p(mxpart,4)

!     generate phase space: pa+pb -> Q
      p=0._dp
      S=sqrts**2
      Q2=mQ**2
!     Generate pa+pb -> Q
      tau=Q2/S
!      ymin=-abs(half*log(tau))
!      ymin=-abs(log((0.999999_dp-tcut/sqrts)/sqrt(tau)))
      ymin=-abs(log((one-tcut/sqrts)/sqrt(tau)))
      ymax=-ymin
      wt=one
      call pick(1,y,ymin,ymax,r1,wt)
      wt=wt/S
      xa=min(sqrt(tau)*exp(y),one)
      xb=min(sqrt(tau)*exp(-y),one)

      p(1,4)=0.5_dp*sqrts*xa
      p(1,3)=p(1,4)
      p(2,4)=0.5_dp*sqrts*xb
      p(2,3)=-p(2,4)
      p(3,:)=p(1,:)+p(2,:)
      wt=wt*two*pi
!     Add flux factor
      wt=wt/xa/xb/sqrts/sqrts
      return
      end

      subroutine genparton(n,p,r1,r2,r3,r4,tcut,ph,wt)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      integer branch,n
      real(dp):: r1,r2,r3,r4,wt,tcut
      real(dp):: Q(4),pah(4),pbh(4),Qh(4),pr(4)
      real(dp):: p(mxpart,4),ph(mxpart,4)
!
!     Radiate an additional massless partice (from p_branch).
!     using p1+p2 --> Q --> p3+...+p_n to generate ph1+ph2 --> Qh --> ph3+...+ph{n}+ph{n+1}
      Q(:)=p(1,:)+p(2,:)
      branch=1+int(r4*n)
      if (branch.le.2) then
         call gen_init(branch,Q,r1,r2,r3,tcut,pah,pbh,Qh,pr,wt)
         ph(1,:)=pah(:)
         ph(2,:)=pbh(:)
         ph(n+2,:)=pr(:)
         if (n.eq.2) then
            ph(3,:)=Qh(:)
         else
!     lorentz boost Q=p3+...+pn to Qh=ph3+...+ph{n}
         endif
      else
!     call gen_final()
      endif
      wt=wt*n! We should sum over the n events. But we pick one randomly, so multiply weight with n
      return
      end
      
      subroutine genparton2(n,p,r1,r2,r3,r4,r5,r6,r7,
     & t0,ph,wt,pass)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      logical pass
      integer b1,b2,n
      real(dp):: r(6),wt,t0,r1,r2,r3,r4,r5,r6,r7
      real(dp):: Q(4),pah(4),pbh(4),Qh(4),pr1(4),pr2(4)
      real(dp):: p(mxpart,4),ph(mxpart,4)
      
      r(1)=r1
      r(2)=r2
      r(3)=r3
      r(4)=r4
      r(5)=r5
      r(6)=r6

!
!     Radiate 2 additional massless partices (from p_branch).
!     using p1+p2 --> Q --> p3+...+p_n to generate ph1+ph2 --> Qh --> ph3+...+ph{n}+ph{n+1}+ph{n+2}
      Q(:)=p(1,:)+p(2,:)
!      b1=1+int(rand()*n)
!      b2=1+int(rand()*n)
      if (n /= 2) then
        write(6,*) 'genparton2 only conceived for n=2; n=',n
        stop
      endif
      b1=1+int(r7*4)
      if (b1 < 3) then
        b2=1
      else
        b1=b1-2
        b2=2
      endif
      if ((b1.le.2).and.(b2.le.2)) then
         call gen_init2(b1,b2,Q,r,t0,pah,pbh,Qh,pr1,pr2,wt,pass)
         if (pass) then
            ph(1,:)=pah(:)
            ph(2,:)=pbh(:)
            ph(n+2,:)=pr1(:)
            ph(n+3,:)=pr2(:)
            if (n.eq.2) then
               ph(3,:)=Qh(:)
            else
!     lorentz boost Q=p3+...+pn to Qh=ph3+...+ph{n}
            endif
         endif
      else
        write(6,*) 'No final state branching implemented yet'
        stop
!     call gen_final()
      endif
!      wt=wt*n*n! We should sum over the n events. But we pick one randomly, so multiply weight with n
      return
      end

      subroutine gen_init(branch,Q,r1,r2,r3,tcut,pah,pbh,Qh,pr,wt)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'energy.f'
      integer branch
      real(dp):: wt,tcut,xa,xb,ta,tb,r1,r2,r3,Jac
      real(dp):: phi,pt,al,xah,xbh,prDQ
      real(dp):: Q(4),pab(4),Qh(4),pah(4),pbh(4),pr(4)

      xa=(Q(4)+Q(3))/sqrts
      xb=(Q(4)-Q(3))/sqrts
!     pick tau_a and tau_b of emitted particle
      if (branch.eq.1) then
         call pick_tau(tcut,tcut,xb,xa,r1,tb,wt)
         if (tb.le.sqrts*(1-xb)*(1-xa*xb)/(2-xa-xb)) then
            call pick_tau(tcut,tb,xa,xb,r2,ta,wt)
         else
            call pick_tau1(tcut,tb,xa,xb,r2,ta,wt)
         endif
      endif
      if (branch.eq.2) then
         call pick_tau(tcut,tcut,xa,xb,r1,ta,wt)
         if (ta.le.sqrts*(1-xa)*(1-xa*xb)/(2-xa-xb)) then
            call pick_tau(tcut,ta,xb,xa,r2,tb,wt)
         else
            call pick_tau1(tcut,ta,xb,xa,r2,tb,wt)
         endif
      endif
      wt=0.25_dp*wt   
      call pick(1,phi,0._dp,2._dp*pi,r3,wt)
!     construct pr
      pt=sqrt(abs(ta*tb))
      pr(4)=0.5_dp*(ta+tb)
      pr(3)=0.5_dp*(ta-tb)
      pr(2)=pt*sin(phi)
      pr(1)=pt*cos(phi)
!     construct intial state particles
      prDQ=(pr(4)*Q(4)-pr(3)*Q(3))/pt**2
      al=prDQ-1._dp-sqrt(abs(prDQ**2+1._dp))
      Pab(4)=Q(4)-al*pr(4)
      Pab(3)=Q(3)-al*pr(3)
      Pab(2)=0._dp
      Pab(1)=0._dp
      xah=min((Pab(4)+Pab(3))/sqrts,1._dp)
      xbh=min((Pab(4)-Pab(3))/sqrts,1._dp)
!      write(*,*) xah,xbh,branch,ta,tb
      pah(4)=0.5_dp*sqrts*xah
      pah(3)=0.5_dp*sqrts*xah
      pah(2)=0._dp
      pah(1)=0._dp
      pbh(4)=0.5_dp*sqrts*xbh
      pbh(3)=-0.5_dp*sqrts*xbh
      pbh(2)=0._dp
      pbh(1)=0._dp
!     construct Qhat
      Qh(4)=Q(4)-(al+1._dp)*pr(4)
      Qh(3)=Q(3)-(al+1._dp)*pr(3)
      Qh(2)=-pr(2)
      Qh(1)=-pr(1)
!     calculate Jacobian
      Jac=(pr(4)*Qh(4)-pr(3)*Qh(3))**2-pt**4
      Jac=Jac/((pr(4)*Q(4)-pr(3)*Q(3))**2+pt**4)
      Jac=sqrt(abs(Jac))/(2._dp*pi)**3
      wt=wt*Jac
!     Change flux factor
      wt=wt*(xa/xah)*(xb/xbh)
      return
      end

      subroutine gen_init2(b1,b2,Q,r,t0,pah,pbh,Qh,pr1,pr2,wt,pass)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'energy.f'
      logical pass
      integer b1,b2
      real(dp):: wt,xa,xb,ta1,tb1,ta2,tb2,Jac
      real(dp):: phi1,phi2,pt1,pt2,pl2,al,xah,xbh,prDQ
      real(dp):: r1,r2,r3,r4,t0,soma,somb
      real(dp):: Q(4),pab(4),Qh(4),pah(4),pbh(4),r(6)
      real(dp):: pr(4),pr1(4),pr2(4)
!
! a requirement is t0<tcut/2
!
      pass=.false.
      xa=(Q(4)+Q(3))/sqrts
      xb=(Q(4)-Q(3))/sqrts
      r1=r(1)
      r2=r(2)
      r3=r(3)
      r4=r(4)
      
      soma=sqrts*(1._dp-xa)
      somb=sqrts*(1._dp-xb)
      if (soma.lt.2._dp*t0) return
      if (somb.lt.2._dp*t0) return
      call pick(2,ta1,t0,soma-t0 ,r1,wt)
      call pick(2,tb1,t0,somb-t0 ,r2,wt)
      call pick(2,ta2,t0,soma-ta1,r3,wt)
      call pick(2,tb2,t0,somb-tb1,r4,wt)

!     construct pr1
      wt=0.25_dp*wt
      r1=r(5)
      call pick(1,phi1,0._dp,2._dp*pi,r1,wt)
      pt1=sqrt(abs(ta1*tb1))
      pr1(4)=0.5_dp*(ta1+tb1)
      pr1(3)=0.5_dp*(ta1-tb1)
      pr1(2)=pt1*sin(phi1)
      pr1(1)=pt1*cos(phi1)

!     construct pr2
      wt=0.25_dp*wt
      r1=r(6)
      call pick(1,phi2,0._dp,2._dp*pi,r1,wt)
      pt2=sqrt(abs(ta2*tb2))
      pr2(4)=0.5_dp*(ta2+tb2)
      pr2(3)=0.5_dp*(ta2-tb2)
      pr2(2)=pt2*sin(phi2)
      pr2(1)=pt2*cos(phi2)
!     construct intial state particles
      pr=pr1+pr2
      pt2=pr(1)**2+pr(2)**2
      pl2=pr(4)**2-pr(3)**2
      prDQ=(pr(4)*Q(4)-pr(3)*Q(3))/pl2
      al=prDQ-1._dp-sqrt(abs(prDQ**2+pt2/pl2))
      Pab(4)=Q(4)-al*pr(4)
      Pab(3)=Q(3)-al*pr(3)
      Pab(2)=0._dp
      Pab(1)=0._dp
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
!     construct Qhat
      Qh(4)=Q(4)-(al+1._dp)*pr(4)
      Qh(3)=Q(3)-(al+1._dp)*pr(3)
      Qh(2)=-pr(2)
      Qh(1)=-pr(1)
!     calculate Jacobian
      prDQ=(pr(4)*Q(4)-pr(3)*Q(3))**2
      Jac=sqrt(abs(prDQ/(prDQ+pt2*pl2)))
      wt=wt*Jac/(2._dp*pi)**6
!     Change flux factor
      wt=wt*(xa/xah)*(xb/xbh)
      pass=.true.
      return
      end


      subroutine pick_tau(tcut,tb,xa,xb,r,ta,wt)
      implicit none
      include 'types.f'
      include 'energy.f'
      real(dp):: tcut,ta,tb,xa,xb,r,wt,tamin,tamax,S

      S=sqrts*sqrts
      tamin=tcut
      tamax=tb*(xa-2._dp)+sqrts*(1._dp-xa)*xb
      tamax=(tamax+sqrt(4._dp*sqrts*tb*(1._dp-xa)*xb+tamax**2))/2._dp/xb
      if (tb.gt.tcut) tamax=min(tamax,tb)
      call pick(2,ta,tamin,tamax,r,wt)

      return
      end
      

      subroutine pick_tau1(tcut,ta,xb,xa,r,tb,wt)
      implicit none
      include 'types.f'
      include 'energy.f'
      real(dp):: tcut,ta,tb,xa,xb,r,wt,tbmin,tbmax

      tbmin=tcut
      tbmax=ta*xb*(ta-sqrts*(1._dp-xa))/(ta*(xa-2._dp)+sqrts*(1._dp-xa))
      tbmax=min(tbmax,ta)
      call pick(2,tb,tbmin,tbmax,r,wt)
      
      return
      end

      subroutine pick(itype,s,smin,smax,r1,wt)
      implicit none
      include 'types.f'
      integer itype
      real(dp):: s,smin,smax,r1,wt
*
* itype=1  pick linearly
* itype=2  pick logarithmically
*
      if(itype.eq.1)then
        s =smin+r1*(smax-smin)
        wt=wt*(smax-smin)
      endif
      if(itype.eq.2)then
        s =smin*(smax/smin)**r1
        wt=wt*log(smax/smin)*s
      endif
      return
      end

