      subroutine qqb_zbb_gvec(p,n,in,msq)
      implicit none
      include 'types.f'

C-----Nov 10 99 --- checked that gives the right matrix element
C     when summed over polarizations.

c----Matrix element for Z+bbar production
C----averaged over initial colours and spins
c    line "in" contracted with the vector n(mu)
c     g(-p1)+g(-p2)--> bb(p5)+ b(p6)+Z(f(p3)+af(p4))

c---but this routine contains no factor of 1/2 for identical gluons
c   in the final state.

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'msqv_cs.f'
      include 'mmsqv_cs.f'
      include 'heavyflav.f'
      include 'nflav.f'

C ip is the label of the emitter
C in is the label of the contracted line
      integer:: j,k,pq,pl,in,ics
      real(dp):: fac,n(4)
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart),prop
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),ggqqb(2,2),tmp


      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call spinoru(6,p,za,zb)
C   zab=<i-|k|j-> zba=<i+|k|j+> where k is an arbitrary 4-vector
c---Conventions of Bern, Dixon, Kosower, Weinzierl,
c---ie, za(i,j)*zb(j,i)=s(i,j)
      call spinork(6,p,zab,zba,n)


C---exclude the photon pole, 4*mbsq choosen as a scale approx above upsilon
c      if (s(3,4) < 4._dp*mbsq) return

      do pq=1,2
      do pl=1,2
      ggqqb(pq,pl) =0._dp
      enddo
      enddo

c      write(6,*) 'in in qqb_zbb_gvec',in

C arguments 1-4 represent (1) incoming quark line
C                         (2) incoming quark line
C                         (3) outgoing gluon line
C                         (4) outgoing gluon line contracted with n
      if     (in == 1) then
        call zbbsqn(5,6,2,1,p,n,za,zb,zab,zba,ggqqb)
      elseif (in == 2)  then
        call zbbsqn(5,6,1,2,p,n,za,zb,zab,zba,ggqqb)
c--- since we have interchanged 1 and 2 to get the gg matrix element,
c---  the colour structures should be interchanged too
        do pq=1,2
        do pl=1,2
        tmp=mmsqv_cs(1,pq,pl)
        mmsqv_cs(1,pq,pl)=mmsqv_cs(2,pq,pl)
        mmsqv_cs(2,pq,pl)=tmp
        enddo
        enddo
      endif

      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)
      fac=v*xn/four*(esq*gsq)**2*two

      do pq=1,2
      do pl=1,2
      ggqqb(pq,pl) =avegg*fac*ggqqb(pq,pl)
      enddo
      enddo
      do j=-nflav,nflav
      do k=-nflav,nflav
      if     ((j == 0) .and. (k == 0)) then
       msq(j,k)=
     & +abs(Q(flav)*q1+L(flav)*l1*prop)**2*ggqqb(1,1)
     & +abs(Q(flav)*q1+R(flav)*r1*prop)**2*ggqqb(2,2)
     & +abs(Q(flav)*q1+L(flav)*r1*prop)**2*ggqqb(1,2)
     & +abs(Q(flav)*q1+R(flav)*l1*prop)**2*ggqqb(2,1)
       do ics=0,2
       msqv_cs(ics,j,k)=
     & +abs(Q(flav)*q1+L(flav)*l1*prop)**2*avegg*fac*mmsqv_cs(ics,1,1)
     & +abs(Q(flav)*q1+R(flav)*l1*prop)**2*avegg*fac*mmsqv_cs(ics,2,1)
     & +abs(Q(flav)*q1+L(flav)*r1*prop)**2*avegg*fac*mmsqv_cs(ics,1,2)
     & +abs(Q(flav)*q1+R(flav)*r1*prop)**2*avegg*fac*mmsqv_cs(ics,2,2)
       enddo
      endif

      enddo
      enddo

      return
      end

      subroutine zbbsqn(i1,i2,i5,i6,p,n,za,zb,zab,zba,msq)
      implicit none
      include 'types.f'

C-----Apart from overall factors returns the matrix element squared
C-----msq dependent on the helicities pq and pl of the quark and
C-----lepton lines for
C-----q(-p1)+qbar(-p2)-->l(p3)+al(p4)+g(p5)+g(p6) where
C-----where gluon 6 has been contracted with the vector n
Cargument 1-4 represent (i1) incoming quark line
C                       (i2) incoming quark line
C                       (i5) outgoing gluon line
C                       (i6) outgoing gluon line contracted with n
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'mmsqv_cs.f'
      complex(dp):: qcdabn(2,2,2),qcdban(2,2,2),qedn(2,2,2)
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      real(dp):: msq(2,2),n(4),p(mxpart,4),nDp5
      integer:: i1,i2,i3,i4,i5,i6,pg,pq,pl

      i3=3
      i4=4

      nDp5=n(4)*p(i5,4)-n(3)*p(i5,3)-n(2)*p(i5,2)-n(1)*p(i5,1)

      call checkndotp(p,n,i6)

      call subqcdn(i1,i2,i3,i4,i5,i6,nDp5,za,zb,zab,zba,qcdabn,qcdban)

C--first argument is gluon line
C--second argument is polarization of i5 line pq
C--third argument is polarization of lepton line pl
C  1=L,2=R

      do pq=1,2
      do pl=1,2
      mmsqv_cs(0,pq,pl)=0._dp
      mmsqv_cs(1,pq,pl)=0._dp
      mmsqv_cs(2,pq,pl)=0._dp

      do pg=1,2
      qedn(pg,pq,pl)=qcdabn(pg,pq,pl)+qcdban(pg,pq,pl)
      mmsqv_cs(0,pq,pl)=mmsqv_cs(0,pq,pl)-ninth*abs(qedn(pg,pq,pl))**2
      mmsqv_cs(1,pq,pl)=mmsqv_cs(1,pq,pl)+abs(qcdabn(pg,pq,pl))**2
      mmsqv_cs(2,pq,pl)=mmsqv_cs(2,pq,pl)+abs(qcdban(pg,pq,pl))**2
      enddo
      msq(pq,pl)=mmsqv_cs(0,pq,pl)+mmsqv_cs(1,pq,pl)+mmsqv_cs(2,pq,pl)
      enddo
      enddo

      return
      end



