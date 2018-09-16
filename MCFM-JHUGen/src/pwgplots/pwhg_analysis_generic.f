      block data blockwhcprg
      character * 6 whcprg      
      common/cwhcprg/whcprg
      data whcprg/'NLO   '/
      end
            
      function lastnbx(string)
      include 'types.f'
      character *(*) string
      integer:: l,k
      l=len(string)
      do k=l,1,-1
         if(string(k:k).ne.' ') then
            lastnbx=k
            return
         endif
      enddo
      end
      
      function prodvec2(vec1,vec2)
      include 'types.f'
      real(dp) prodvec2,vec1(4),vec2(4)
      prodvec2=vec1(4)*vec2(4)-vec1(1)*vec2(1)
     1 -vec1(2)*vec2(2)-vec1(3)*vec2(3)
      end


      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      include 'types.f'
      include 'constants.f'
      
      real(dp) p(4),y,eta,pt,mass,pv
      real(dp) tiny
      parameter (tiny=1.e-5_dp)
      y=0.5_dp*log((p(4)+p(3))/(p(4)-p(3)))
      pt=sqrt(p(1)**2+p(2)**2)
      pv=sqrt(pt**2+p(3)**2)
      if(pt<tiny)then
         eta=sign(one,p(3))*1.e8_dp
      else
         eta=0.5_dp*log((pv+p(3))/(pv-p(3)))
      endif
      mass=sqrt(abs(p(4)**2-pv**2))
      end


      subroutine yetaptmassplot(p,dsig,prefix)
      implicit none
      include 'types.f'
      
      real(dp) p(4),dsig
      character *(*) prefix
      real(dp) y,eta,pt,m
      call getyetaptmass(p,y,eta,pt,m)
      call filld(prefix//'_y',y,dsig)
      call filld(prefix//'_eta',eta,dsig)
      call filld(prefix//'_pt',pt,dsig)
      call filld(prefix//'_m',m,dsig)
      end


      subroutine pwhgfinalopshist
      implicit none
      include 'types.f'
      
      include 'pwhg_bookhist-new.h'
      integer::  f_idx,b_idx,x_idx,ixx,jxx,j,k,l
      integer::  indexhist
      real(dp) num,numerr,den,denerr

      f_idx=indexhist('dF-dpT')
      b_idx=indexhist('dB-dpT')
      x_idx=indexhist('dAFB-dpT')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('dF-dpT')
      b_idx=indexhist('dB-dpT')
      x_idx=indexhist('F-gt-pT')
      call integratehist(f_idx,x_idx,1)
      x_idx=indexhist('B-gt-pT')
      call integratehist(b_idx,x_idx,1)
      x_idx=indexhist('F-lt-pT')
      call integratehist(f_idx,x_idx,-1)
      x_idx=indexhist('B-lt-pT')
      call integratehist(b_idx,x_idx,-1)

      f_idx=indexhist('F-gt-pT')
      b_idx=indexhist('B-gt-pT')
      x_idx=indexhist('AFB-gt-pT')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('F-lt-pT')
      b_idx=indexhist('B-lt-pT')
      x_idx=indexhist('AFB-lt-pT')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('dF-dm_tt')
      b_idx=indexhist('dB-dm_tt')
      x_idx=indexhist('dAFB-dm_tt')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('dF-dm_tt')
      b_idx=indexhist('dB-dm_tt')
      x_idx=indexhist('F-gt-m_tt')
      call integratehist(f_idx,x_idx,1)
      x_idx=indexhist('B-gt-m_tt')
      call integratehist(b_idx,x_idx,1)
      x_idx=indexhist('F-lt-m_tt')
      call integratehist(f_idx,x_idx,-1)
      x_idx=indexhist('B-lt-m_tt')
      call integratehist(b_idx,x_idx,-1)

      f_idx=indexhist('F-gt-m_tt')
      b_idx=indexhist('B-gt-m_tt')
      x_idx=indexhist('AFB-gt-m_tt')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('F-lt-m_tt')
      b_idx=indexhist('B-lt-m_tt')
      x_idx=indexhist('AFB-lt-m_tt')
      call get_afb_hist(f_idx,b_idx,x_idx)

      end


      subroutine integratehist(integrand,integral,direction)
      implicit none
      include 'types.f'
      
      include 'pwhg_bookhist-new.h'
      integer::  integrand,integral,direction
      integer::  indexhist,ixx,jxx
      external indexhist
      if(direction==1) then
        do ixx=1,nbins(integrand)
          do jxx=ixx,nbins(integrand)
            yhistarr2(ixx,integral)=yhistarr2(ixx,integral)
     1         +yhistarr2(jxx,integrand)
     1         *(xhistarr(jxx+1,integrand)-xhistarr(jxx,integrand))
            errhistarr2(ixx,integral)=errhistarr2(ixx,integral)
     1         +(errhistarr2(jxx,integrand)
     1         *(xhistarr(jxx+1,integrand)-xhistarr(jxx,integrand)))**2
          enddo
          errhistarr2(ixx,integral)=sqrt(errhistarr2(ixx,integral))
        enddo
      else if(direction==-1) then
        do ixx=1,nbins(integrand)
          do jxx=ixx,1,-1
            yhistarr2(ixx,integral)=yhistarr2(ixx,integral)
     1         +yhistarr2(jxx,integrand)
     1         *(xhistarr(jxx+1,integrand)-xhistarr(jxx,integrand))
            errhistarr2(ixx,integral)=errhistarr2(ixx,integral)
     1         +(errhistarr2(jxx,integrand)
     1         *(xhistarr(jxx+1,integrand)-xhistarr(jxx,integrand)))**2
          enddo
          errhistarr2(ixx,integral)=sqrt(errhistarr2(ixx,integral))
        enddo
      else
         write(*,*) 'subroutine integratehist: error!'
         write(*,*) '--------------------------------'
         write(*,*) 'Invalid direction specified for histogram given.'
         write(*,*) 'direction=1/-1 only - integral from lowest bin to'
         write(*,*) 'highest bin or highest bin to lowest bin.'
         call exit(-1)
      endif

      end


      subroutine get_afb_hist(f_idx,b_idx,afb_idx)
      implicit none
      include 'types.f'
      
      include 'pwhg_bookhist-new.h'
      real(dp)   f,ef,b,eb
      integer::  ixx,f_idx,b_idx,afb_idx
      
      do ixx=1,nbins(f_idx)
         f=yhistarr2(ixx,f_idx)
         ef=errhistarr2(ixx,f_idx)
         b=yhistarr2(ixx,b_idx)
         eb=errhistarr2(ixx,b_idx)
         if((f+b)>0._dp) then         ! Guard against division by zero.
            yhistarr2(ixx,afb_idx)=(f-b)/(f+b)
            errhistarr2(ixx,afb_idx)=2*sqrt((f*ef)**2+(b*eb)**2)
     1                              /(f+b)**2
         else
            yhistarr2(ixx,afb_idx)=0._dp
            errhistarr2(ixx,afb_idx)=0._dp
         endif
      enddo

      end


      function islept(j)
      logical:: islept
      integer:: j
      if(abs(j)>=11.and.abs(j)<=16) then
         islept = .true.
      else
         islept = .false.
      endif
      end

      function phepDot(p_A,p_B)
      include 'types.f'
      real(dp)  phepDot
      real(dp)  p_A(4),p_B(4)
      phepDot=p_A(4)*p_B(4)-p_A(1)*p_B(1)
     1       -p_A(2)*p_B(2)-p_A(3)*p_B(3)
      end

c     calculate the separation in the lego plot between the two momenta
c     p1 and p2 in azi and pseudorapidity
      function rsepn_p(p1,p2)
      include 'types.f'      
      include 'constants.f'      
      real(dp) rsepn_p,p1(0:3),p2(0:3)
      real(dp) eta1,phi1,eta2,phi2
      real(dp) delphi
      real(dp) pseudorapidity,azi
      external pseudorapidity,azi

      phi1 = azi(p1)   
      phi2 = azi(p2)
      eta1 = pseudorapidity(p1)
      eta2 = pseudorapidity(p2)

      delphi = abs(phi1-phi2)
      if (delphi>pi) then
         delphi = 2*pi-delphi
      endif
      if (delphi<0 .or. delphi>pi) then
         print*,' problem in rsepn. delphi = ',delphi
      endif
      rsepn_p = sqrt( (eta1-eta2)**2 + delphi**2 )
      end

      function azi(p)
      include 'types.f'      
      include 'constants.f'      
      real(dp) azi,p(0:3)
      azi = atan(p(2)/p(1))
      if (p(1)<0._dp) then
         if (azi>0._dp) then               
            azi = azi - pi
         else
            azi = azi + pi
         endif
      endif    
      end

      function pseudorapidity(p)
      implicit none
      include 'types.f'      
      real(dp) p(0:3),pseudorapidity
      real(dp) mod, costh
      mod = sqrt(p(1)**2+p(2)**2+p(3)**2)
      costh = p(3)/mod
      pseudorapidity=0.5_dp*log((1._dp+costh)/(1._dp-costh))
      end

      subroutine buildjets(mjets,kt,eta,rap,phi,pjet,ibtag,
     &                                               isForClustering)
      implicit none
      include 'types.f'
c     arrays to reconstruct jets
      integer::   mjets
      real(dp)::kt(mjets),eta(mjets),rap(mjets),
     & phi(mjets),pjet(4,mjets)
      integer::   ibtag(mjets)
      logical::   isForClustering(mjets)
      include  'hepevt.h'
      integer::   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=20)
      real(dp)  ptrack(4,maxtrack),pj(4,maxjet)
      integer::   jetvec(maxtrack),itrackhep(maxtrack)
      integer::   ntracks,njets
      integer::   j,k,mu
      real(dp)  r,palg,ptmin,pp,tmp
C - Initialize arrays and counters for output jets
      do j=1,maxtrack
         do mu=1,4
            ptrack(mu,j)=0._dp
         enddo
         jetvec(j)=0
      enddo      
      ntracks=0
      do j=1,mjets
         do mu=1,4
            pjet(mu,j)=0._dp
            pj(mu,j)=0._dp
         enddo
      enddo
C - Extract final state particles to feed to jet finder
      do j=1,nhep
         if(.not.isForClustering(j)) cycle
         if(ntracks==maxtrack) then
            write(*,*) 'analyze: need to increase maxtrack!'
            write(*,*) 'ntracks: ',ntracks
            stop
         endif
         ntracks=ntracks+1
         do mu=1,4
            ptrack(mu,ntracks)=phep(mu,j)
         enddo
         itrackhep(ntracks)=j
      enddo
      if (ntracks==0) then
         return
      endif
C --------------- C
C - Run FastJet - C
C --------------- C
C - R = 0.7   radius parameter
C - f = 0.75  overlapping fraction
      palg  = -1
      r     = 0.5_dp
      ptmin = 10._dp
c      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,jetvec)
      mjets=min(mjets,njets)
      if(njets==0) return
c check consistency
      do k=1,ntracks
         if(jetvec(k)>0) then
            do mu=1,4
               pj(mu,jetvec(k))=pj(mu,jetvec(k))+ptrack(mu,k)
            enddo
         endif
      enddo
      tmp=0
      do j=1,mjets
         do mu=1,4
            tmp=tmp+abs(pj(mu,j)-pjet(mu,j))
         enddo
      enddo
      if(tmp>1d-4) then
         write(*,*) ' bug!'
      endif
C --------------------------------------------------------------------- C
C - Computing arrays of useful kinematics quantities for hardest jets - C
C --------------------------------------------------------------------- C
      do j=1,mjets
         kt(j)=sqrt(pjet(1,j)**2+pjet(2,j)**2)
         pp = sqrt(kt(j)**2+pjet(3,j)**2)
         eta(j)=0.5_dp*log((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
         rap(j)=0.5_dp*log((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
         phi(j)=atan2(pjet(2,j),pjet(1,j))
      enddo
      end
