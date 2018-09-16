      subroutine qqb_ww_g(P,msq)
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  f'(p5)+bar{f'}(p6) + n(p3)+ebar(p4)+ g(p7)
c   for the moment --- radiation only from initial line
      implicit none 
      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'zerowidth.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'anomcoup.f'
      include 'xanomcoup.f'
      include 'plabel.f'
      include 'pchoice.f'
      integer jk,tjk,polg,polq,minus,mplus,jp,kp,jtype
      parameter(minus=1,mplus=2)
      double precision P(mxpart,4),qdks(mxpart,4),msq(-nf:nf,-nf:nf),
     . ave,s127,fac,fac1,offsh,xfac
      double complex ct(2,2),cs_z(2,2),cs_g(2,2),
     . cgamz(2,2),cz(2,2)
      double complex u_ub(5,2,2),d_db(5,2,2),ub_u(5,2,2),db_d(5,2,2),
     .               u_g(5,2,2), d_g(5,2,2), g_ub(5,2,2),g_db(5,2,2),
     .               ub_g(5,2,2),db_g(5,2,2),g_u(5,2,2),g_d(5,2,2),
     .               amp(5),propwp,propwm,propzg,prop12,cprop,A(2,2)
      double precision, parameter :: mp(nf)=(/-1d0,+1d0,-1d0,+1d0,-1d0/)

      do jp=-nf,nf
      do kp=-nf,nf
      msq(jp,kp)=0d0
      enddo
      enddo

      fac=gw**4
      fac1=two*gsq*cf
C---multiply by factor for c-sbar+u-dbar hadronic decay
      if (plabel(5) .eq. 'qj') fac1=2d0*xn*fac1

C----Change the momenta to DKS notation 
C   swapped possibility if we want to swap momenta for hadronic case
c   We have --- f(p1) + f'(p2)-->mu^-(p3)+nubar(p4)+e^+(p6)+nu(p5)+g(p7)
c   DKS have--- ubar(q1)+u(q2)-->mu^-(q3)+nubar(q4)+e^+(q5)+nu(q6)+g(p7)
c----
C   or normal configuration
c   We have --- f(p1) + f'(p2)-->mu^-(p5)+nubar(p6)+e^+(p4)+nu(p3)+g(p7)
c   DKS have--- ubar(q1)+u(q2)-->mu^-(q3)+nubar(q4)+e^+(q5)+nu(q6)+g(p7)

      if ((plabel(5) .eq. 'qj') .and. (plabel(3) .eq. 'el')) then 
c----swapped case
      do j=1,4
      qdks(1,j)=p(1,j)
      qdks(2,j)=p(2,j)
      qdks(3,j)=p(3,j)
      qdks(4,j)=p(4,j)
      qdks(5,j)=p(6,j)
      qdks(6,j)=p(5,j)
      qdks(7,j)=p(7,j)
      enddo
      else
c----normal case
      do j=1,4
      qdks(1,j)=p(1,j)
      qdks(2,j)=p(2,j)
      qdks(3,j)=p(5,j)
      qdks(4,j)=p(6,j)
      qdks(5,j)=p(4,j)
      qdks(6,j)=p(3,j)
      qdks(7,j)=p(7,j)
      enddo
      endif

c--   s returned from sprodx (common block) is 2*dot product
      call spinoru(7,qdks,za,zb)

c--   calculate propagators
      s127=s(1,2)+s(1,7)+s(2,7)
      if     (zerowidth  .eqv. .true.) then
      prop12=s127/dcmplx(s127-zmass**2,zmass*zwidth)  
      cprop=1d0    
      elseif (zerowidth .neqv. .true.) then
      prop12=dcmplx(s127/(s127-zmass**2))
      offsh=s(3,4)-wmass**2
      propwp=dcmplx(offsh)/dcmplx(offsh,wmass*wwidth)
      offsh=s(5,6)-wmass**2
      propwm=dcmplx(offsh)/dcmplx(offsh,wmass*wwidth)
      offsh=s127-zmass**2
      propzg=dcmplx(offsh)/dcmplx(offsh,zmass*zwidth)
      cprop=propwp*propwm*propzg
      endif

c-- couplings according to 3.4 and 3.6
      do j=1,2
      ct(minus,j)=1d0
      ct(mplus,j)=0d0
      cs_z(minus,j)=+mp(j)*l(j)*sin2w*prop12
      cs_z(mplus,j)=-mp(j)*2d0*Q(j)*xw*prop12
      cs_g(minus,j)=+mp(j)*2d0*Q(j)*xw
      cs_g(mplus,j)=+mp(j)*2d0*Q(j)*xw
      cz(minus,j)=0d0
      cz(mplus,j)=0d0
      cgamz(minus,j)=0d0
      cgamz(mplus,j)=0d0
c-- couplings with or without photon pole
      if (zerowidth .neqv. .true.) then
      cz(minus,j)=two*xw*ln*L(j)*prop12
      cz(mplus,j)=two*xw*ln*R(j)*prop12
      cgamz(minus,j)=two*xw*(-Q(j)+le*L(j)*prop12)
      cgamz(mplus,j)=two*xw*(-Q(j)+le*R(j)*prop12)
      endif
      enddo
c-- 
c      l(j)=(tau(j)-two*Q(j)*xw)/sin2w ; r(j)=(-two*Q(j)*xw)/sin2w
c      le=(-1d0-two*(-1d0)*xw)/sin2w ; re=(-two*(-1d0)*xw)/sin2w
c      ln=(+1d0-two*(+0d0)*xw)/sin2w ; rn=0d0
c---

c--- apply a dipole form factor to anomalous couplings (only if tevscale > 0)
      if (tevscale .gt. 0d0) then
        xfac=1d0/(1d0+s127/(tevscale*1d3)**2)**2
      else
        xfac=1d0
      endif
      xdelg1_z=xfac*delg1_z
      xdelg1_g=xfac*delg1_g
      xdelk_z=xfac*delk_z
      xdelk_g=xfac*delk_g
      xlambda_z=xfac*lambda_z
      xlambda_g=xfac*lambda_g
      
c---remember ub-u is the basic process.
c---case ubar-u
      call wwamps(1,2,3,4,5,6,7,za,zb,ub_u)
c---case ubar-g
      call wwamps(1,7,3,4,5,6,2,za,zb,ub_g)
c---case g-ubar
      call wwamps(2,7,3,4,5,6,1,za,zb,g_ub)

c---case u-ubar
      call wwamps(2,1,3,4,5,6,7,za,zb,u_ub)
c---case u-g
      call wwamps(7,1,3,4,5,6,2,za,zb,u_g)
c---case g-u
      call wwamps(7,2,3,4,5,6,1,za,zb,g_u)

c---case dbar-d
      call wwamps(1,2,6,5,4,3,7,za,zb,db_d)
c---case dbar-g
      call wwamps(1,7,6,5,4,3,2,za,zb,db_g)
c---case g-dbar
      call wwamps(2,7,6,5,4,3,1,za,zb,g_db)

c---case d-dbar
      call wwamps(2,1,6,5,4,3,7,za,zb,d_db)
c---case d-g
      call wwamps(7,1,6,5,4,3,2,za,zb,d_g)
c---case g-d
      call wwamps(7,2,6,5,4,3,1,za,zb,g_d)


      do j=-nf,nf
      do k=-nf,nf
c-- skip gluon-gluon case
      if((j .eq. 0).and.(k .eq. 0)) goto 19
c-- skip non-diagonal quark flavors except gluon
      if((j .ne. 0 .and. k .ne. 0) .and. (j .ne. -k)) goto 19
      jk=max(j,k)
      ave=xn*aveqq
      if (j .eq. 0 .or. k .eq. 0) then
          jk=j+k
          ave=xn*aveqg
      endif
      do polg=1,2
      do polq=1,2
c---sum is over diagram type t,s(Z),e,n,s(photon)
      do jtype=1,5
          if    (j .lt. 0 .and. tau(jk) .eq. -1d0 .and. k .ne. 0) then
            amp(jtype)=db_d(jtype,polg,polq)
          elseif(j .lt. 0 .and. tau(jk) .eq.  1d0 .and. k .ne. 0) then
            amp(jtype)=ub_u(jtype,polg,polq)
          elseif(j .gt. 0 .and. tau(jk) .eq. -1d0 .and. k .ne. 0) then
            amp(jtype)=d_db(jtype,polg,polq)
          elseif(j .gt. 0 .and. tau(jk) .eq.  1d0 .and. k .ne. 0) then
            amp(jtype)=u_ub(jtype,polg,polq)
          elseif(j .eq. 0 .and. tau(jk) .eq.  1d0 .and. jk .gt. 0) then
            amp(jtype)=g_u(jtype,polg,polq)
          elseif(j .eq. 0 .and. tau(jk) .eq. -1d0 .and. jk .gt. 0) then
            amp(jtype)=g_d(jtype,polg,polq)
          elseif(j .eq. 0 .and. tau(jk) .eq. -1d0 .and. jk .lt. 0) then
            amp(jtype)=g_ub(jtype,polg,polq)
          elseif(j .eq. 0 .and. tau(jk) .eq.  1d0 .and. jk .lt. 0) then
            amp(jtype)=g_db(jtype,polg,polq)
          elseif(k .eq. 0 .and. tau(jk) .eq.  1d0 .and. jk .gt. 0) then
            amp(jtype)=u_g(jtype,polg,polq)
          elseif(k .eq. 0 .and. tau(jk) .eq. -1d0 .and. jk .gt. 0) then
            amp(jtype)=d_g(jtype,polg,polq)
          elseif(k .eq. 0 .and. tau(jk) .eq. -1d0 .and. jk .lt. 0) then
            amp(jtype)=ub_g(jtype,polg,polq)
          elseif(k .eq. 0 .and. tau(jk) .eq.  1d0 .and. jk .lt. 0) then 
            amp(jtype)=db_g(jtype,polg,polq)
        endif
      enddo

C---tjk is equal to 2 (u,c) or 1 (d,s,b)
        tjk=2-mod(abs(jk),2)

      A(polg,polq)=dcmplx(fac)*cprop
     . *(ct(polq,tjk)*amp(1)+cs_z(polq,tjk)*amp(2)+cs_g(polq,tjk)*amp(5)
     .  +cz(polq,tjk)*amp(3)+cgamz(polq,tjk)*amp(4))
          
      enddo
      enddo

      msq(j,k)=fac1*ave*
     . (cdabs(A(mplus,minus))**2+cdabs(A(minus,minus))**2
     . +cdabs(A(mplus,mplus))**2+cdabs(A(minus,mplus))**2)

        
   19 continue
      enddo
      enddo
      return
      end

      
      




