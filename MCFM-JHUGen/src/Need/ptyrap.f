      double precision function pt(j,p)
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt      
      implicit none
      include 'constants.f'
      integer j
      double precision p(mxpart,4),getet
      
      pt=getet(p(j,4),p(j,1),p(j,2),p(j,3))

      return
      end

      double precision function pttwo(j,k,p)
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt      
      implicit none
      include 'constants.f'
      integer j,k
      double precision p(mxpart,4),getet

      pttwo=getet(p(j,4)+p(k,4),p(j,1)+p(k,1),
     .            p(j,2)+p(k,2),p(j,3)+p(k,3))

      return
      end

      double precision function ptthree(j,k,m,p)
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt      
      implicit none
      include 'constants.f'
      integer j,k,m
      double precision p(mxpart,4),getet
      
      ptthree=getet(p(j,4)+p(k,4)+p(m,4),p(j,1)+p(k,1)+p(m,1),
     .              p(j,2)+p(k,2)+p(m,2),p(j,3)+p(k,3)+p(m,3))

      return
      end

      double precision function ptfour(j,k,m,n,p)
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt      
      implicit none
      include 'constants.f'
      integer j,k,m,n
      double precision p(mxpart,4),getet
           
      ptfour=getet(p(j,4)+p(k,4)+p(m,4)+p(n,4),
     &             p(j,1)+p(k,1)+p(m,1)+p(n,1),
     &             p(j,2)+p(k,2)+p(m,2)+p(n,2),
     &             p(j,3)+p(k,3)+p(m,3)+p(n,3))

      return
      end

      double precision function ptsix(j1,k1,m1,j2,k2,m2,p)
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt      
      implicit none
      include 'constants.f'
      integer j1,k1,m1,j2,k2,m2
      double precision p(mxpart,4),getet
           
      ptsix=getet(p(j1,4)+p(k1,4)+p(m1,4)+p(j2,4)+p(k2,4)+p(m2,4),
     &            p(j1,1)+p(k1,1)+p(m1,1)+p(j2,1)+p(k2,1)+p(m2,1),
     &            p(j1,2)+p(k1,2)+p(m1,2)+p(j2,2)+p(k2,2)+p(m2,2),
     &            p(j1,3)+p(k1,3)+p(m1,3)+p(j2,3)+p(k2,3)+p(m2,3))

      return
      end

      double precision function etarap(j,p)
      implicit none
C---returns the value of the pseudorapidity
      include 'constants.f'
      integer j
      double precision p(mxpart,4)
      etarap=dsqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)
      etarap=(etarap+p(j,3))/(etarap-p(j,3))

      if (etarap .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      etarap=100d0
      else
      etarap=0.5d0*dlog(etarap)
      endif
      return
      end

      double precision function aetarap(j,p)
      implicit none
C---returns the absolute value of the pseudorapidity
      include 'constants.f'
      integer j
      double precision p(mxpart,4)
      aetarap=dsqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)
      aetarap=(aetarap+p(j,3))/(aetarap-p(j,3))
      if (aetarap .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      aetarap=100d0
      else
      aetarap=0.5d0*abs(dlog(aetarap))
      endif
      return
      end
 
      double precision function yrap(j,p)
      implicit none
C---returns the value of the rapidity
      include 'constants.f'
      integer j
      double precision p(mxpart,4)
      yrap=(p(j,4)+p(j,3))/(p(j,4)-p(j,3))
      if (yrap .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrap=100d0
      else
      yrap=0.5d0*dlog(yrap)
      endif
      return
      end

      double precision function ayrap(j,p)
      implicit none
C---returns the absolute value of the rapidity
      include 'constants.f'
      integer j
      double precision p(mxpart,4)
      ayrap=(p(j,4)+p(j,3))/(p(j,4)-p(j,3))
      if (ayrap .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      ayrap=100d0
      else
      ayrap=0.5d0*dabs(dlog(ayrap))
      endif
      return
      end
 
c--- this is the rapidity of pair j,k
      double precision function yraptwo(j,k,p)
      implicit none
      include 'constants.f'
      integer j,k
      double precision p(mxpart,4)
      yraptwo=(p(j,4)+p(k,4)+p(j,3)+p(k,3))
     .       /(p(j,4)+p(k,4)-p(j,3)-p(k,3))
      if (yraptwo .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yraptwo=100d0
      else 
      yraptwo=0.5d0*dlog(yraptwo)
      endif
            
      return
      end

c--- this is the pseudo-rapidity of pair j,k
      double precision function etaraptwo(j,k,p)
      implicit none
      include 'constants.f'
      integer j,k
      double precision p(mxpart,4)
      
      etaraptwo=dsqrt((p(j,1)+p(k,1))**2+(p(j,2)+p(k,2))**2
     .               +(p(j,3)+p(k,3))**2)
      if (abs(etaraptwo)-abs(p(j,3)+p(k,3)) .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      etaraptwo=100d0
      else 
      etaraptwo=(etaraptwo+p(j,3)+p(k,3))
     .         /(etaraptwo-p(j,3)-p(k,3))
      etaraptwo=0.5d0*dlog(etaraptwo)
      endif
      
      return
      end

      double precision function yrapthree(j,k,m,p)
c--- this is the rapidity of the combination j+k+m
      implicit none
      include 'constants.f'
      integer j,k,m
      double precision p(mxpart,4)
      yrapthree=(p(j,4)+p(k,4)+p(m,4)+p(j,3)+p(k,3)+p(m,3))
     .         /(p(j,4)+p(k,4)+p(m,4)-p(j,3)-p(k,3)-p(m,3))
      if (yrapthree .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrapthree=100d0
      else 
      yrapthree=0.5d0*dlog(yrapthree)
      endif
            
      return
      end

      double precision function etarapthree(j,k,m,p)
c--- this is the pseudo-rapidity of the combination j+k+m
      implicit none
      include 'constants.f'
      integer j,k,m
      double precision p(mxpart,4)
      
      etarapthree=
     .    dsqrt((p(j,1)+p(k,1)+p(m,1))**2+(p(j,2)+p(k,2)+p(m,2))**2
     .         +(p(j,3)+p(k,3)+p(m,3))**2)
      if (abs(etarapthree)-abs(p(j,3)+p(k,3)+p(m,3)) .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      etarapthree=100d0
      else 
      etarapthree=(etarapthree+p(j,3)+p(k,3)+p(m,3))
     .           /(etarapthree-p(j,3)-p(k,3)-p(m,3))
      etarapthree=0.5d0*dlog(etarapthree)
      endif
      
      return
      end

      double precision function yrapfour(j,k,m,n,p)
c--- this is the rapidity of the combination j+k+m+n
      implicit none
      include 'constants.f'
      integer j,k,m,n
      double precision p(mxpart,4)
      yrapfour=(p(j,4)+p(k,4)+p(m,4)+p(n,4)+p(j,3)+p(k,3)+p(m,3)+p(n,3))
     &        /(p(j,4)+p(k,4)+p(m,4)+p(n,4)-p(j,3)-p(k,3)-p(m,3)-p(n,3))
      if (yrapfour .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrapfour=100d0
      else 
      yrapfour=0.5d0*dlog(yrapfour)
      endif
            
      return
      end

      double precision function yrapsix(j1,k1,m1,j2,k2,m2,p)
c--- this is the rapidity of the combination j1+k1+m1+j2+k2+m2
      implicit none
      include 'constants.f'
      integer j1,k1,m1,j2,k2,m2
      double precision p(mxpart,4)
      yrapsix=(p(j1,4)+p(k1,4)+p(m1,4)+p(j1,3)+p(k1,3)+p(m1,3)
     &        +p(j2,4)+p(k2,4)+p(m2,4)+p(j2,3)+p(k2,3)+p(m2,3))
     &       /(p(j1,4)+p(k1,4)+p(m1,4)-p(j1,3)-p(k1,3)-p(m1,3)
     &        +p(j2,4)+p(k2,4)+p(m2,4)-p(j2,3)-p(k2,3)-p(m2,3))
      if (yrapsix .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrapsix=100d0
      else 
      yrapsix=0.5d0*dlog(yrapsix)
      endif
            
      return
      end

      double precision function yrapseven(j1,k1,m1,j2,k2,m2,n,p)
c--- this is the rapidity of the combination j1+k1+m1+j2+k2+m2+n
      implicit none
      include 'constants.f'
      integer j1,k1,m1,j2,k2,m2,n
      double precision p(mxpart,4)
      yrapseven=(p(j1,4)+p(k1,4)+p(m1,4)+p(j1,3)+p(k1,3)+p(m1,3)+p(n,4)
     &          +p(j2,4)+p(k2,4)+p(m2,4)+p(j2,3)+p(k2,3)+p(m2,3)+p(n,3))
     &         /(p(j1,4)+p(k1,4)+p(m1,4)-p(j1,3)-p(k1,3)-p(m1,3)+p(n,4)
     &          +p(j2,4)+p(k2,4)+p(m2,4)-p(j2,3)-p(k2,3)-p(m2,3)-p(n,3))
      if (yrapseven .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrapseven=100d0
      else 
      yrapseven=0.5d0*dlog(yrapseven)
      endif
       return
       end
c -- RR: mass definitions

       double precision function onemass(j,p)
       implicit none
       include 'constants.f'
       integer j
       double precision p(mxpart,4)
       onemass=p(j,4)**2-p(j,1)**2-p(j,2)**2-p(j,3)**2
       onemass=dsqrt(onemass)
       return
       end

      double precision function twomass(j,k1,p)
      implicit none
      include 'constants.f'
      integer j,k1
      double precision p(mxpart,4),pjk(4)
      pjk(:)=p(j,:)+p(k1,:)
      twomass=pjk(4)**2-pjk(1)**2-pjk(2)**2-pjk(3)**2
      twomass=dsqrt(twomass)
      return
      end

      double precision function threemass(j,k1,k2,p)
      implicit none
      include 'constants.f'
      integer j,k1,k2
      double precision p(mxpart,4),ptot(4)
      ptot(:)=p(j,:)+p(k1,:)+p(k2,:)
      threemass=ptot(4)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2
      threemass=dsqrt(threemass)
      return
      end

      double precision function fourmass(j,k1,k2,k3,p)
      implicit none
      include 'constants.f'
      integer j,k1,k2,k3
      double precision p(mxpart,4),ptot(4)
      ptot(:)=p(j,:)+p(k1,:)+p(k2,:)+p(k3,:)
      fourmass=ptot(4)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2
      fourmass=dsqrt(fourmass)
      return
      end

      

      


