      subroutine qqb_tbb_g(p,msq)
      implicit none
c     Matrix element for real corrections to single top production
C     (nwz=+1)
c      u(-p1)+dbar(-p2)-->t(=> n(p3)+e^+(p4)+b(p5))+F(p6)+G(p7)
C     or for
C     (nwz=-1)
c      ubar(-p1)+d(-p2)-->t~(=> e^-(p3)+n(p4)+bbar(p5))+F(p6)+G(p7)
C     averaged(summed) over initial(final) colours and spins

c--- G(p7) represents either a gluon or a light quark
c---  For the case isub=1, F(p6) represents a light quark
c---   and thus this corresponds to t-channel W-exchange
c---  For the case isub=2, F(p6) represents b-bar and thus
c---   this corresponds to s-channel production

c--- added by JC, 4/9/08: common block "stopbmass" holds a logical
c---  flag "masslessb" that determines whether the extra b~ that
c---  appears in the real matrix elements is massless or not:
c---
c---       u -------------- d
c---               $
c---               $
c---               $
c---    g ~~~~~------------ t
c---           \
c---            \ 
c---              b~


      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'zprods_com.f'
      include 'stopbmass.f'
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),
     . ubtdg_l,ubtdg_h
      double precision fac,qqb,qbq,qg,gq,qbg,gqb,ub,bu,ubarb,bubar,gb,bg
      integer j,k,isub

      common/isub/isub

      call spinoru(7,p,za,zb)
      fac=2d0*gsq*cf*gw**8*xn**2
      
      if (nwz .eq. +1) then
c--- t production
      
      if     (isub .eq. 1) then
        qqb=0d0
        qbq=0d0
        qg=0d0
        gq=0d0
        qbg=0d0
        gqb=0d0
        ub=aveqq*fac
     .   *(ubtdg_l(1,2,3,4,5,6,7,p)+ubtdg_h(1,2,3,4,5,6,7,p))
        bu=aveqq*fac
     .   *(ubtdg_l(2,1,3,4,5,6,7,p)+ubtdg_h(2,1,3,4,5,6,7,p))
        ubarb=aveqq*fac
     .   *(ubtdg_l(6,2,3,4,5,1,7,p)+ubtdg_h(6,2,3,4,5,1,7,p))
        bubar=aveqq*fac
     .   *(ubtdg_l(6,1,3,4,5,2,7,p)+ubtdg_h(6,1,3,4,5,2,7,p))
c--- Note that the two '_h' radiation terms below correspond to
c---  diagrams of the form g+b-->W(s+c)+t(->W+b) and are a large
c---  contribution at the LHC that is not included in B. Harris et al.  
        bg=aveqg*fac
     .   *(ubtdg_l(7,1,3,4,5,6,2,p)+zip*ubtdg_h(7,1,3,4,5,6,2,p))
        gb=aveqg*fac
     .   *(ubtdg_l(7,2,3,4,5,6,1,p)+zip*ubtdg_h(7,2,3,4,5,6,1,p))
c In this version of the matrix elements, the b-quark is in
c position 6, for potential merging with the s-channel process
c        qg= aveqg*fac
c     .   *(zip*ubtdg_l(1,6,3,4,5,7,2,p)+ubtdg_h(1,6,3,4,5,7,2,p))
c        gq= aveqg*fac
c     .   *(zip*ubtdg_l(2,6,3,4,5,7,1,p)+ubtdg_h(2,6,3,4,5,7,1,p))
c        qbg=aveqg*fac
c     .   *(zip*ubtdg_l(7,6,3,4,5,1,2,p)+ubtdg_h(7,6,3,4,5,1,2,p))
c        gqb=aveqg*fac
c     .   *(zip*ubtdg_l(7,6,3,4,5,2,1,p)+ubtdg_h(7,6,3,4,5,2,1,p))
c In this version of the matrix elements, the b-quark is in
c position 7 with label "pp", for counting the extra b-quark
c contribution - Z. Sullivan 1/25/05
        if (masslessb) then
        qg= aveqg*fac
     .   *(zip*ubtdg_l(1,7,3,4,5,6,2,p)+ubtdg_h(1,7,3,4,5,6,2,p))
        gq= aveqg*fac
     .   *(zip*ubtdg_l(2,7,3,4,5,6,1,p)+ubtdg_h(2,7,3,4,5,6,1,p))
        qbg=aveqg*fac
     .   *(zip*ubtdg_l(6,7,3,4,5,1,2,p)+ubtdg_h(6,7,3,4,5,1,2,p))
        gqb=aveqg*fac
     .   *(zip*ubtdg_l(6,7,3,4,5,2,1,p)+ubtdg_h(6,7,3,4,5,2,1,p))
        endif
      elseif (isub .eq. 2) then
        qqb=aveqq*fac
     .   *(ubtdg_l(1,6,3,4,5,2,7,p)+ubtdg_h(1,6,3,4,5,2,7,p))
        qbq=aveqq*fac
     .   *(ubtdg_l(2,6,3,4,5,1,7,p)+ubtdg_h(2,6,3,4,5,1,7,p))
        qg= aveqg*fac
     .   *(ubtdg_l(1,6,3,4,5,7,2,p)+zip*ubtdg_h(1,6,3,4,5,7,2,p))
        gq= aveqg*fac
     .   *(ubtdg_l(2,6,3,4,5,7,1,p)+zip*ubtdg_h(2,6,3,4,5,7,1,p))
        qbg=aveqg*fac
     .   *(ubtdg_l(7,6,3,4,5,1,2,p)+zip*ubtdg_h(7,6,3,4,5,1,2,p))
        gqb=aveqg*fac
     .   *(ubtdg_l(7,6,3,4,5,2,1,p)+zip*ubtdg_h(7,6,3,4,5,2,1,p))
        ub=0d0
        bu=0d0
        ubarb=0d0
        bubar=0d0
        bg=0d0
        gb=0d0
      else
        write(6,*) 'Value of isub is wrong in qqb_tbb_g.f: isub=',isub
        stop
      endif
      
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0

c--- Q-Qbar
      if     ((j .gt. 0) .and. (k .lt. 0)) then
        if (j .eq. 5) then
          msq(j,k)=(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)
     .             +Vsq(+4,k)+Vsq(+5,k))*bubar
        else        
          msq(j,k)=Vsq(j,k)*qqb 
        endif        
c--- Qbar-Q
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
        if (k .eq. 5) then
          msq(j,k)=(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)
     .             +Vsq(j,+4)+Vsq(j,+5))*ubarb
        else
          msq(j,k)=Vsq(j,k)*qbq 
        endif
c--- Q-Q
      elseif ((j .eq. 5) .and. (k .gt. 0)) then
        msq(j,k)=(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)
     .           +Vsq(-4,k)+Vsq(-5,k))*bu
      elseif ((j .gt. 0) .and. (k .eq. 5)) then
        msq(j,k)=(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)
     .           +Vsq(j,-4)+Vsq(j,-5))*ub    
c--- g-Q
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
        if (k .eq. 5) then
          msq(j,k)=2d0*gb 
        else
          msq(j,k)=(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)
     .             +Vsq(-4,k)+Vsq(-5,k))*gq 
        endif
c--- g-Qbar
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
        msq(j,k)=(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)
     .           +Vsq(+4,k)+Vsq(+5,k))*gqb 
c--- Q-g
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
        if (j .eq. 5) then
          msq(j,k)=2d0*bg
        else
          msq(j,k)=(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)
     .             +Vsq(j,-4)+Vsq(j,-5))*qg
        endif
c--- Qbar-g
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
        msq(j,k)=(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)
     .           +Vsq(j,+4)+Vsq(j,+5))*qbg 
      endif

c      write(6,*) j,k,msq(j,k)

      enddo
      enddo

c      pause
      elseif (nwz .eq. -1) then
c--- t~ production
      
      if     (isub .eq. 1) then
        qqb=0d0
        qbq=0d0
        qg=0d0
        gq=0d0
        qbg=0d0
        gqb=0d0
        ub=aveqq*fac
     .   *(ubtdg_l(6,2,4,3,5,1,7,p)+ubtdg_h(6,2,4,3,5,1,7,p))
        bu=aveqq*fac
     .   *(ubtdg_l(6,1,4,3,5,2,7,p)+ubtdg_h(6,1,4,3,5,2,7,p))
        ubarb=aveqq*fac
     .   *(ubtdg_l(1,2,4,3,5,6,7,p)+ubtdg_h(1,2,4,3,5,6,7,p))
        bubar=aveqq*fac
     .   *(ubtdg_l(2,1,4,3,5,6,7,p)+ubtdg_h(2,1,4,3,5,6,7,p))
c--- Note that the two '_h' radiation terms below correspond to
c---  diagrams of the form g+b-->W(s+c)+t(->W+b) and are a large
c---  contribution at the LHC that is not included in B. Harris et al.  
        bg=aveqg*fac
     .   *(ubtdg_l(6,1,4,3,5,7,2,p)+zip*ubtdg_h(6,1,4,3,5,7,2,p))
        gb=aveqg*fac
     .   *(ubtdg_l(6,2,4,3,5,7,1,p)+zip*ubtdg_h(6,2,4,3,5,7,1,p))
c In this version of the matrix elements, the b-quark is in
c position 6, for potential merging with the s-channel process
c        qg= aveqg*fac
c     .   *(zip*ubtdg_l(7,6,4,3,5,1,2,p)+ubtdg_h(7,6,4,3,5,1,2,p))
c        gq= aveqg*fac
c     .   *(zip*ubtdg_l(7,6,4,3,5,2,1,p)+ubtdg_h(7,6,4,3,5,2,1,p))
c        qbg=aveqg*fac
c     .   *(zip*ubtdg_l(1,6,4,3,5,7,2,p)+ubtdg_h(1,6,4,3,5,7,2,p))
c        gqb=aveqg*fac
c     .   *(zip*ubtdg_l(2,6,4,3,5,7,1,p)+ubtdg_h(2,6,4,3,5,7,1,p))
c In this version of the matrix elements, the b-quark is in
c position 7 with label "pp", for counting the extra b-quark
c contribution - Z. Sullivan 1/25/05
        if (masslessb) then
        qg= aveqg*fac
     .   *(zip*ubtdg_l(6,7,4,3,5,1,2,p)+ubtdg_h(6,7,4,3,5,1,2,p))
        gq= aveqg*fac
     .   *(zip*ubtdg_l(6,7,4,3,5,2,1,p)+ubtdg_h(6,7,4,3,5,2,1,p))
        qbg=aveqg*fac
     .   *(zip*ubtdg_l(1,7,4,3,5,6,2,p)+ubtdg_h(1,7,4,3,5,6,2,p))
        gqb=aveqg*fac
     .   *(zip*ubtdg_l(2,7,4,3,5,6,1,p)+ubtdg_h(2,7,4,3,5,6,1,p))
        endif
      elseif (isub .eq. 2) then
        qqb=aveqq*fac
     .   *(ubtdg_l(2,6,4,3,5,1,7,p)+ubtdg_h(2,6,4,3,5,1,7,p))
        qbq=aveqq*fac
     .   *(ubtdg_l(1,6,4,3,5,2,7,p)+ubtdg_h(1,6,4,3,5,2,7,p))
        qg= aveqg*fac
     .   *(ubtdg_l(7,6,4,3,5,1,2,p)+zip*ubtdg_h(7,6,4,3,5,1,2,p))
        gq= aveqg*fac
     .   *(ubtdg_l(7,6,4,3,5,2,1,p)+zip*ubtdg_h(7,6,4,3,5,2,1,p))
        qbg=aveqg*fac
     .   *(ubtdg_l(1,6,4,3,5,7,2,p)+zip*ubtdg_h(1,6,4,3,5,7,2,p))
        gqb=aveqg*fac
     .   *(ubtdg_l(2,6,4,3,5,7,1,p)+zip*ubtdg_h(2,6,4,3,5,7,1,p))
        ub=0d0
        bu=0d0
        ubarb=0d0
        bubar=0d0
        bg=0d0
        gb=0d0
      else
        write(6,*) 'Value of isub is wrong in qqb_tbb_g.f: isub=',isub
        stop
      endif
      
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0

c--- Q-Qbar
      if     ((j .gt. 0) .and. (k .lt. 0)) then
        if (k .eq. -5) then
          msq(j,k)=(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)
     .             +Vsq(j,-4)+Vsq(j,-5))*ub
        else        
          msq(j,k)=Vsq(j,k)*qqb 
        endif        
c--- Qbar-Q
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
        if (j .eq. -5) then
          msq(j,k)=(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)
     .             +Vsq(-4,k)+Vsq(-5,k))*bu
        else
          msq(j,k)=Vsq(j,k)*qbq 
        endif
c--- Qbar-Qbar
      elseif ((j .eq. -5) .and. (k .lt. 0)) then
        msq(j,k)=(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)
     .           +Vsq(+4,k)+Vsq(+5,k))*bubar
      elseif ((j .lt. 0) .and. (k .eq. -5)) then
        msq(j,k)=(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)
     .           +Vsq(j,+4)+Vsq(j,+5))*ubarb    
c--- g-Q
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)
     .             +Vsq(-4,k)+Vsq(-5,k))*gq 
c--- g-Qbar
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
        if (k .eq. -5) then
          msq(j,k)=2d0*gb 
        else
        msq(j,k)=(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)
     .           +Vsq(+4,k)+Vsq(+5,k))*gqb 
        endif
c--- Q-g
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)
     .             +Vsq(j,-4)+Vsq(j,-5))*qg
c--- Qbar-g
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
        if (j .eq. -5) then
          msq(j,k)=2d0*bg
        else
        msq(j,k)=(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)
     .           +Vsq(j,+4)+Vsq(j,+5))*qbg 
        endif
      endif

c      write(6,*) j,k,msq(j,k)

      enddo
      enddo

      endif
      
c      pause

      return
      end


      double precision function ubtdg_l(ju,jb,jn,je,jc,jd,jg,p)
C     Matrix element squared for single top production with gluon 
C     radiation in production
C      u(ju) b(jb) -> t(n~(jn)+e+(je)+c(jc))+d(jd)+g(jg)
C     masses of b quarks c.c=b.b=0
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer ju,jb,jn,je,jc,jd,jg,nu
      double precision p(mxpart,4),pt(4),ptDpt
      double precision sne,sdug,prop
      double complex ampi(2)

      do nu=1,4
      pt(nu)=p(je,nu)+p(jn,nu)+p(jc,nu)
      enddo
      ptDpt=pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2

      sne=s(jn,je)
      sdug=s(jd,ju)+s(jd,jg)+s(ju,jg)
      
      if (sdug .lt. 0d0) then
      prop=(sdug-wmass**2)**2
      else
      prop=(sdug-wmass**2)**2+(wmass*wwidth)**2
      endif
      prop=((sne-wmass**2)**2+(wmass*wwidth)**2)
     .    *((ptDpt-mt**2)**2+(mt*twidth)**2)*prop
      
C  -Lefthanded gluon
      ampi(1)=(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
     . *zb(ju,jd)/(zb(jg,ju)*zb(jg,jd))
     . -(zb(je,jc)*za(jc,jg)+zb(je,jn)*za(jn,jg))/zb(jg,jd)
      ampi(1)=ampi(1)*za(jc,jn)*zb(ju,jb)

C  -Righthanded gluon
      ampi(2)=za(ju,jd)/(za(ju,jg)*za(jd,jg))
     . -zb(jg,jb)/(zb(ju,jb)*za(ju,jg))
      ampi(2)=-ampi(2)*za(jc,jn)*zb(ju,jb)
     . *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))

      
      ubtdg_l=(cdabs(ampi(1))**2+cdabs(ampi(2))**2)/prop
      return
      end


      double precision function ubtdg_h(ju,jb,jn,je,jc,jd,jg,p)
C     Matrix element squared for single top production with gluon 
C     radiation in production (radiation from final line)
C      u(ju) b(jb) -> t(n~(jn)+e+(je)+c(jc))+d(jd)+g(jg)
C     masses of b quarks c.c=b.b=0
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer ju,jb,jn,je,jc,jd,jg,nu
      double precision p(mxpart,4),pt(4),ptDpt
      double precision sne,sdu,prop,twoptg
      double complex ampf(2)

      do nu=1,4
      pt(nu)=p(je,nu)+p(jn,nu)+p(jc,nu)
      enddo
      ptDpt=pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2
      twoptg=
     . 2d0*(pt(4)*p(jg,4)-pt(1)*p(jg,1)-pt(2)*p(jg,2)-pt(3)*p(jg,3))

      sne=s(jn,je)
      sdu=s(jd,ju)
      if (sdu .lt. 0d0) then
      prop=(sdu-wmass**2)**2
      else
      prop=(sdu-wmass**2)**2+(wmass*wwidth)**2
      endif
      prop=((sne-wmass**2)**2+(wmass*wwidth)**2)
     .    *((ptDpt-mt**2)**2+(mt*twidth)**2)*prop
      
C  -Lefthanded gluon tb-line
      ampf(1)=
     . -(zb(je,jc)*za(jc,jg)+zb(je,jn)*za(jn,jg))*zb(jb,ju)*za(ju,jd)
     . +mt**2*zb(je,jb)*za(jg,jd)
      ampf(1)=za(jc,jn)*zb(ju,jb)/(twoptg*zb(jg,jb))*ampf(1)

c---eikonal form
c      ampf(1)=
c     . -(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
c     . *(za(jg,ju)*zb(ju,jb)+za(jg,jd)*zb(jd,jb))
c     . /zb(jg,jb)
c     . -(zb(je,jc)*za(jc,jg)+zb(je,jn)*za(jn,jg))
c     . *za(jg,jd)
c      ampf(1)=za(jc,jn)*zb(ju,jb)/twoptg*ampf(1)

C  -Righthanded gluon tb-line
c---eikonal form
      ampf(2)=
     . +(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))*zb(ju,jg)/zb(ju,jb)
     . -(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
     . *(zb(jg,ju)*za(ju,jb)+zb(jg,jd)*za(jd,jb))/twoptg
      ampf(2)=za(jc,jn)*zb(ju,jb)/za(jb,jg)*ampf(2)

C---alternative form
c      ampf(2)=
c     . +(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))*zb(ju,jg)/zb(ju,jb)
c     . +(-(zb(je,jc)*za(jc,jb)+zb(je,jn)*za(jn,jb))
c     .   *(zb(jg,ju)*za(ju,jd)+zb(jg,jb)*za(jb,jd))
c     . +mt**2*zb(je,jg)*za(jb,jd))/twoptg
c      ampf(2)=za(jc,jn)*zb(ju,jb)/za(jb,jg)*ampf(2)


      ubtdg_h=(cdabs(ampf(1))**2+cdabs(ampf(2))**2)/prop

      return
      end
