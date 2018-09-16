      subroutine genrff(p,i1,i2,i7,r1,r2,phit,wt5_4,*)
c----i1 is an initial state vector.
c---final,final--i1 is the emitter
      implicit none
      include 'constants.f'
      include 'debug.f'
      integer i1,i2,i7,j
      double precision p(mxpart,4),rtalbe,c(4),d(4),phi,phit,wt5_4
      double precision dot,r1,r2,jacbit
      double precision beta,alpha
      double precision oma,y,omy,wt0
      parameter(wt0=1d0/eight/pisq)

c----final-final case
      
c      call writeout(p)
      phi=twopi*phit

C---for p soft set r1=0.9999d0 and r2=0.4999999999d0

c      write(6,*) 'setting r1=1d-7'
c      r1=9.999999999999d0
c      write(6,*) 'setting r2=0.49999999'
c      r2=1d-7


      y=r1**2
      omy=one-y
      if (r2 .le. 0.5d0) oma=(two*r2)**2
      if (r2 .gt. 0.5d0) oma=one-(two*r2-one)**2


      alpha=one-oma
      beta=y*oma
      rtalbe=sqrt(beta*alpha)
c      write(6,*) 'i1',i1
c      write(6,*) 'i2',i2
c      write(6,*) 'i7',i7
c      write(6,*) 'alpha,beta,rtalbe',alpha,beta,rtalbe
      jacbit=four*sqrt(y)/(half/sqrt(oma)+half/sqrt(alpha))

      wt5_4=wt0*omy*dot(p,4,5)*jacbit

      if (debug) write(6,*) 'wt5_4 in genrff',wt5_4

c---generate n-1 momenta from n momenta
c---sudakov in terms of original vectors
      call gtperp(rtalbe,p,i1,i2,3,c,d)
      do j=1,4
c--define p(i7)
      p(i7,j)=alpha*p(i1,j)+beta*p(i2,j)+cos(phi)*c(j)+sin(phi)*d(j)
c--now-define modified p4 and p5
      p(i1,j)=p(i1,j)+y*p(i2,j)-p(i7,j)
      p(i2,j)=p(i2,j)*omy
      enddo
      call writeout(p)

      end



