      subroutine qqb_tbb_v(p,msqv)
      implicit none
c     Virtual matrix element for t-bbar production
C     (nwz=+1)
c      u(-p1)+dbar(-p2)-->n(p3)+e^+(p4)+b(p5)+bbar(p6)
C     or for
C     (nwz=-1)
c      ubar(-p1)+d(-p2)-->e^-(p3)+n(p4)+bbar(p5)+b(p6)
C     averaged(summed) over initial(final) colours and spins
c--NB average over spins only -- colour factors cancel
      integer j,k
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'zprods_com.f'
      include 'scheme.f'
      double precision p(mxpart,4),msqv(-nf:nf,-nf:nf),
     . qqb,qbq,fac,virttop

      scheme='dred'

      call spinoru(6,p,za,zb)
C-----Using our standard practice that virtuals are in units as/4/pi
C-----and by multiplying only by ason2pi we effectively include the
C-----factor of 2 for interference.
      fac=ason2pi*cf
      fac=aveqq*xn**2*gw**8*fac 
              
      if     (nwz .eq. +1) then
        qqb=fac*virttop(1,6,3,4,5,2,za,zb)
        qbq=fac*virttop(2,6,3,4,5,1,za,zb)
c        qqb=fac*virtqqb(1,6,3,4,5,2)
c        qbq=fac*virtqqb(2,6,3,4,5,1)
      elseif (nwz .eq. -1) then
        qqb=fac*virttop(2,6,4,3,5,1,zb,za)
        qbq=fac*virttop(1,6,4,3,5,2,zb,za)
c        qqb=fac*virtqqb(2,6,4,3,5,1)
c        qbq=fac*virtqqb(1,6,4,3,5,2)
      endif

      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0d0
      if     ((j .gt. 0) .and. (k .lt. 0)) then
      msqv(j,k)=Vsq(j,k)*qqb 
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
      msqv(j,k)=Vsq(j,k)*qbq 
      endif
      enddo
      enddo

      return
      end


      double precision function virtqqb(ju,jb,jn,je,jc,jd)
      implicit none

      integer ju,jd,jn,je,jc,jb
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      double precision snec,prop,cv,cv0,mtsq
      double complex c1,amp,ampho


      mtsq=mt**2
      snec=+s(jn,je)+s(je,jc)+s(jc,jn)

      call coefs(s(ju,jd),mtsq,cv0,cv,c1)

      if (s(ju,jd) .lt. 0d0) then
      prop=(s(ju,jd)-wmass**2)**2
      else
      prop=(s(ju,jd)-wmass**2)**2+(wmass*wwidth)**2
      endif
      prop=prop*((snec-mt**2)**2+(mt*twidth)**2)
      prop=prop*((s(jn,je)-wmass**2)**2+(wmass*wwidth)**2)

      amp=za(jc,jn)*zb(ju,jb)
     . *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
      ampho=za(jc,jn)*zb(ju,jb)
     . *(dcmplx(cv0+cv)*(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
     . +c1*dcmplx(0.5d0)*zb(je,jb)*za(jb,jd))

      virtqqb=dble(amp*dconjg(ampho))/prop

      return
      end



      subroutine coefs(s12,mtsq,cv0,cv,c1)
C-----In this routine:-
C-----cv0 is the results for all massless vertex function
C-----cv and c1 are  is the results for one-mass vertex function
      implicit none
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'scheme.f'
      double precision cv,cv0,Li2la
      double precision s12,mtsq,taucs,ddilog,eta,la,oml
      double complex lnrat,logoml,logla,xl12,logsca,Kfun,c1

      if (scheme .eq.'dred') then
C------        eta=0 4d-hel
         eta=0d0
      elseif (scheme .eq. 'tH-V') then
C------       eta=1 t'Hooft Veltman
         eta=1d0
      endif

C**********************************************************************
C   Massless case
C   Taken from
C   %\cite{Altarelli:1979ub}
C   \bibitem{Altarelli:1979ub}
C   G.~Altarelli, R.~K.~Ellis and G.~Martinelli,
C   %``Large Perturbative Corrections To The Drell-Yan Process In QCD,''
C   Nucl.\ Phys.\ B {\bf 157}, 461 (1979).
C   %%CITATION = NUPHA,B157,461;%%
C   Using Eqn(58) with normalization changed to 
C   as/2/pi*cf*(4*pi)^ep/Gamma(1-ep) 
C   Taking account that Gamma(1-ep)^2/Gamma(1-2*ep)=1-ep^2*pi^2/6
C**********************************************************************
      xl12=lnrat(-s12,musq) 

C-----2/22/2012
c-----This appears to be the correction to the vertex in units of as/4/pi*cf
c-----despite the above comment
      cv0=-2d0*epinv*(epinv2-dble(xl12))-dble(xl12**2)
     .           -3d0*(epinv-dble(xl12))-7d0-eta



C---- this routine has been constructed following closely 
C---- the notation of
C---- %\cite{Gottschalk:1980rv}
C---- \bibitem{Gottschalk:1980rv}
C---- T.~Gottschalk,
C---- %``Chromodynamic Corrections To Neutrino Production Of Heavy Quarks,''
C---- Phys.\ Rev.\ D {\bf 23}, 56 (1981).
C---- %%CITATION = PHRVA,D23,56;%%
C----- Adapted from Eqs.(A8,A9)

C NB  s12=-Q^2, -taucs=mtsq+Q^2
      taucs=s12-mtsq
      la=-s12/(mtsq-s12)
      oml=1d0-la
C-----oml=mtsq/(mtsq-s12)
      logoml=-lnrat(-taucs,mtsq)
      logsca=lnrat(-taucs,musq)
      Kfun=dcmplx(oml/la)*logoml

c--- Minus sign relative to Gottschalk since incoming b has momentum
c--- vector reversed for the t-channel process
c--- s-channel process follows by crossing
      c1=-dcmplx(2d0)*Kfun

      if (la .lt. 1d0) then
      Li2la=ddilog(la)
      else
      logla=lnrat(-s12,-taucs)
      Li2la=pisqo6-ddilog(oml)-dble(logla*logoml)
      endif
C-----Again from A8 and A9 these are in units of alpha_s/4/pi*CF
      cv=-epinv*epinv2
     . -epinv*(2.5d0+dble(logoml-logsca))
     . -0.5d0*(11d0+eta)-pisqo6+2d0*Li2la-dble(Kfun)
     .  -0.5d0*dble(logoml*(cone-logoml))
     .  +2.5d0*dble(logsca)+dble(logsca*logoml)-0.5d0*dble(logsca**2)

C answer from gotts.mac
c   ans:
c   -1/ep^2;
c   -2.5/ep -(log(oml)/ep-+log(sca)/ep)
c  -6-pisqo6+2*Li2-Ka
c  -0.5*log(oml)*(1-log(oml))
c  +2.5*log(sca)+log(oml)*log(sca)-log(sca)^2/2
      return
      end

      double precision function virttop(ju,jb,jn,je,jc,jd,za,zb)
      implicit none

      integer ju,jd,jn,je,jc,jb
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      double precision snec,prop,cv,cv0,mtsq
      double complex c1,amp,ampho


      mtsq=mt**2
      snec=+s(jn,je)+s(je,jc)+s(jc,jn)

      call coefs(s(ju,jd),mtsq,cv0,cv,c1)

      if (s(ju,jd) .lt. 0d0) then
      prop=(s(ju,jd)-wmass**2)**2
      else
      prop=(s(ju,jd)-wmass**2)**2+(wmass*wwidth)**2
      endif
      prop=prop*((snec-mt**2)**2+(mt*twidth)**2)
      prop=prop*((s(jn,je)-wmass**2)**2+(wmass*wwidth)**2)

      amp=za(jc,jn)*zb(ju,jb)
     . *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
      ampho=za(jc,jn)*zb(ju,jb)
     . *(dcmplx(cv0+cv)*(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
     . +c1*dcmplx(0.5d0)*zb(je,jb)*za(jb,jd))

      virttop=dble(amp*dconjg(ampho))/prop

      return
      end


