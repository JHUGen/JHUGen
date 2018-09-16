      function pt(j,p)
      implicit none
      include 'types.f'
      real(dp):: pt
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt      
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j
      real(dp):: p(mxpart,4),getet
      
      pt=getet(p(j,4),p(j,1),p(j,2),p(j,3))

      return
      end

      function pttwo(j,k,p)
      implicit none
      include 'types.f'
      real(dp):: pttwo
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt      
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k
      real(dp):: p(mxpart,4),getet

      pttwo=getet(p(j,4)+p(k,4),p(j,1)+p(k,1),
     &            p(j,2)+p(k,2),p(j,3)+p(k,3))

      return
      end

      function ptthree(j,k,m,p)
      implicit none
      include 'types.f'
      real(dp):: ptthree
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt      
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k,m
      real(dp):: p(mxpart,4),getet
      
      ptthree=getet(p(j,4)+p(k,4)+p(m,4),p(j,1)+p(k,1)+p(m,1),
     &              p(j,2)+p(k,2)+p(m,2),p(j,3)+p(k,3)+p(m,3))

      return
      end

      function ptfour(j,k,m,n,p)
      implicit none
      include 'types.f'
      real(dp):: ptfour
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt      
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k,m,n
      real(dp):: p(mxpart,4),getet
           
      ptfour=getet(p(j,4)+p(k,4)+p(m,4)+p(n,4),
     &             p(j,1)+p(k,1)+p(m,1)+p(n,1),
     &             p(j,2)+p(k,2)+p(m,2)+p(n,2),
     &             p(j,3)+p(k,3)+p(m,3)+p(n,3))

      return
      end

      function ptsix(j1,k1,m1,j2,k2,m2,p)
      implicit none
      include 'types.f'
      real(dp):: ptsix
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt      
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j1,k1,m1,j2,k2,m2
      real(dp):: p(mxpart,4),getet
           
      ptsix=getet(p(j1,4)+p(k1,4)+p(m1,4)+p(j2,4)+p(k2,4)+p(m2,4),
     &            p(j1,1)+p(k1,1)+p(m1,1)+p(j2,1)+p(k2,1)+p(m2,1),
     &            p(j1,2)+p(k1,2)+p(m1,2)+p(j2,2)+p(k2,2)+p(m2,2),
     &            p(j1,3)+p(k1,3)+p(m1,3)+p(j2,3)+p(k2,3)+p(m2,3))

      return
      end

      function etarap(j,p)
      implicit none
      include 'types.f'
      real(dp):: etarap
      
C---returns the value of the pseudorapidity
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j
      real(dp):: p(mxpart,4)
      etarap=sqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)
      etarap=(etarap+p(j,3))/(etarap-p(j,3))

      if (etarap < 1.e-13_dp) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      etarap=100._dp
      else
      etarap=0.5_dp*log(etarap)
      endif
      return
      end

      function aetarap(j,p)
      implicit none
      include 'types.f'
      real(dp):: aetarap
      
C---returns the absolute value of the pseudorapidity
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j
      real(dp):: p(mxpart,4)
      aetarap=sqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)
      aetarap=(aetarap+p(j,3))/(aetarap-p(j,3))
      if (aetarap < 1.e-13_dp) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      aetarap=100._dp
      else
      aetarap=0.5_dp*abs(log(aetarap))
      endif
      return
      end
 
      function yrap(j,p)
      implicit none
      include 'types.f'
      real(dp):: yrap
      
C---returns the value of the rapidity
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j
      real(dp):: p(mxpart,4)
      yrap=(p(j,4)+p(j,3))/(p(j,4)-p(j,3))
      if (yrap < 1.e-13_dp) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrap=100._dp
      else
      yrap=0.5_dp*log(yrap)
      endif
      return
      end

      function ayrap(j,p)
      implicit none
      include 'types.f'
      real(dp):: ayrap
      
C---returns the absolute value of the rapidity
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j
      real(dp):: p(mxpart,4)
      ayrap=(p(j,4)+p(j,3))/(p(j,4)-p(j,3))
      if (ayrap < 1.e-13_dp) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      ayrap=100._dp
      else
      ayrap=0.5_dp*abs(log(ayrap))
      endif
      return
      end
 
c--- this is the rapidity of pair j,k
      function yraptwo(j,k,p)
      implicit none
      include 'types.f'
      real(dp):: yraptwo
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k
      real(dp):: p(mxpart,4)
      yraptwo=(p(j,4)+p(k,4)+p(j,3)+p(k,3))
     &       /(p(j,4)+p(k,4)-p(j,3)-p(k,3))
      if (yraptwo < 1.e-13_dp) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yraptwo=100._dp
      else 
      yraptwo=0.5_dp*log(yraptwo)
      endif
            
      return
      end

c--- this is the pseudo-rapidity of pair j,k
      function etaraptwo(j,k,p)
      implicit none
      include 'types.f'
      real(dp):: etaraptwo
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k
      real(dp):: p(mxpart,4)
      
      etaraptwo=sqrt((p(j,1)+p(k,1))**2+(p(j,2)+p(k,2))**2
     &               +(p(j,3)+p(k,3))**2)
      if (abs(etaraptwo)-abs(p(j,3)+p(k,3)) < 1.e-13_dp) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      etaraptwo=100._dp
      else 
      etaraptwo=(etaraptwo+p(j,3)+p(k,3))
     &         /(etaraptwo-p(j,3)-p(k,3))
      etaraptwo=0.5_dp*log(etaraptwo)
      endif
      
      return
      end

      function yrapthree(j,k,m,p)
      implicit none
      include 'types.f'
      real(dp):: yrapthree
c--- this is the rapidity of the combination j+k+m
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k,m
      real(dp):: p(mxpart,4)
      yrapthree=(p(j,4)+p(k,4)+p(m,4)+p(j,3)+p(k,3)+p(m,3))
     &         /(p(j,4)+p(k,4)+p(m,4)-p(j,3)-p(k,3)-p(m,3))
      if (yrapthree < 1.e-13_dp) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrapthree=100._dp
      else 
      yrapthree=0.5_dp*log(yrapthree)
      endif
            
      return
      end

      function etarapthree(j,k,m,p)
      implicit none
      include 'types.f'
      real(dp):: etarapthree
c--- this is the pseudo-rapidity of the combination j+k+m
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k,m
      real(dp):: p(mxpart,4)
      
      etarapthree=
     &    sqrt((p(j,1)+p(k,1)+p(m,1))**2+(p(j,2)+p(k,2)+p(m,2))**2
     &         +(p(j,3)+p(k,3)+p(m,3))**2)
      if (abs(etarapthree)-abs(p(j,3)+p(k,3)+p(m,3)) < 1.e-13_dp) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      etarapthree=100._dp
      else 
      etarapthree=(etarapthree+p(j,3)+p(k,3)+p(m,3))
     &           /(etarapthree-p(j,3)-p(k,3)-p(m,3))
      etarapthree=0.5_dp*log(etarapthree)
      endif
      
      return
      end

      function yrapfour(j,k,m,n,p)
      implicit none
      include 'types.f'
      real(dp):: yrapfour
c--- this is the rapidity of the combination j+k+m+n
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k,m,n
      real(dp):: p(mxpart,4)
      yrapfour=(p(j,4)+p(k,4)+p(m,4)+p(n,4)+p(j,3)+p(k,3)+p(m,3)+p(n,3))
     &        /(p(j,4)+p(k,4)+p(m,4)+p(n,4)-p(j,3)-p(k,3)-p(m,3)-p(n,3))
      if (yrapfour < 1.e-13_dp) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrapfour=100._dp
      else 
      yrapfour=0.5_dp*log(yrapfour)
      endif
            
      return
      end

      function yrapsix(j1,k1,m1,j2,k2,m2,p)
      implicit none
      include 'types.f'
      real(dp):: yrapsix
c--- this is the rapidity of the combination j1+k1+m1+j2+k2+m2
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j1,k1,m1,j2,k2,m2
      real(dp):: p(mxpart,4)
      yrapsix=(p(j1,4)+p(k1,4)+p(m1,4)+p(j1,3)+p(k1,3)+p(m1,3)
     &        +p(j2,4)+p(k2,4)+p(m2,4)+p(j2,3)+p(k2,3)+p(m2,3))
     &       /(p(j1,4)+p(k1,4)+p(m1,4)-p(j1,3)-p(k1,3)-p(m1,3)
     &        +p(j2,4)+p(k2,4)+p(m2,4)-p(j2,3)-p(k2,3)-p(m2,3))
      if (yrapsix < 1.e-13_dp) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrapsix=100._dp
      else 
      yrapsix=0.5_dp*log(yrapsix)
      endif
            
      return
      end

      function yrapseven(j1,k1,m1,j2,k2,m2,n,p)
      implicit none
      include 'types.f'
      real(dp):: yrapseven
c--- this is the rapidity of the combination j1+k1+m1+j2+k2+m2+n
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j1,k1,m1,j2,k2,m2,n
      real(dp):: p(mxpart,4)
      yrapseven=(p(j1,4)+p(k1,4)+p(m1,4)+p(j1,3)+p(k1,3)+p(m1,3)+p(n,4)
     &          +p(j2,4)+p(k2,4)+p(m2,4)+p(j2,3)+p(k2,3)+p(m2,3)+p(n,3))
     &         /(p(j1,4)+p(k1,4)+p(m1,4)-p(j1,3)-p(k1,3)-p(m1,3)+p(n,4)
     &          +p(j2,4)+p(k2,4)+p(m2,4)-p(j2,3)-p(k2,3)-p(m2,3)-p(n,3))
      if (yrapseven < 1.e-13_dp) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrapseven=100._dp
      else 
      yrapseven=0.5_dp*log(yrapseven)
      endif
       return
       end
c -- RR: mass definitions

       function onemass(j,p)
       implicit none
      include 'types.f'
      real(dp):: onemass
       
       include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
       integer:: j
       real(dp):: p(mxpart,4)
       onemass=p(j,4)**2-p(j,1)**2-p(j,2)**2-p(j,3)**2
       onemass=sqrt(onemass)
       return
       end

      function twomass(j,k1,p)
      implicit none
      include 'types.f'
      real(dp):: twomass
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k1
      real(dp):: p(mxpart,4),pjk(4)
      pjk(:)=p(j,:)+p(k1,:)
      twomass=pjk(4)**2-pjk(1)**2-pjk(2)**2-pjk(3)**2
      twomass=sqrt(twomass)
      return
      end

      function threemass(j,k1,k2,p)
      implicit none
      include 'types.f'
      real(dp):: threemass
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k1,k2
      real(dp):: p(mxpart,4),ptot(4)
      ptot(:)=p(j,:)+p(k1,:)+p(k2,:)
      threemass=ptot(4)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2
      threemass=sqrt(threemass)
      return
      end

      function fourmass(j,k1,k2,k3,p)
      implicit none
      include 'types.f'
      real(dp):: fourmass
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k1,k2,k3
      real(dp):: p(mxpart,4),ptot(4)
      ptot(:)=p(j,:)+p(k1,:)+p(k2,:)+p(k3,:)
      fourmass=ptot(4)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2
      fourmass=sqrt(fourmass)
      return
      end

      

      


