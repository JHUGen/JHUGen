      subroutine genrif(p,i1,i2,i7,r1,r2,phit,wt5_4,*)
      implicit none
      include 'types.f'
c----i1 is an initial state vector.
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'debug.f'
      include 'xmin.f'
      integer:: i1,i2,i7,j
      real(dp):: p(mxpart,4),rtalbe,c(4),d(4),phi,phit,wt5_4,
     & jacbit,dot,r1,r2
      real(dp):: beta,alpha
      real(dp):: x,omx,omxmin,wt0
      parameter(wt0=1._dp/eight/pisq)

c------initial final

      phi=twopi*phit
      omxmin=one-xmin

c -xx(i1)
      omx=omxmin*r1**2
      x=one-omx
      if (r2 < 0.5_dp) alpha=(two*r2)**2
      if (r2 > 0.5_dp) alpha=one-(two*r2-one)**2

      jacbit=four*r1/(half/sqrt(alpha)+half/sqrt(one-alpha))

      wt5_4=wt0*dot(p,i1,i2)/x**2*omxmin*jacbit
      if (debug) write(6,*) 'i1 in genrif',i1
      if (debug) write(6,*) 'i2 in genrif',i2
      if (debug) write(6,*) 'x in genrif',x
      if (debug) write(6,*) 'omxmin in genrif',omxmin
      if (debug) write(6,*) 'wt5_4 in genrif',wt5_4
      beta=omx*(1._dp-alpha)/x
      rtalbe=sqrt(beta*alpha)
c---generate transverse vectors c and d with length^2=rtalbe^2*2*p1Dp2
c-- with direction in transverse plane picked by 6
      call gtperp(rtalbe,p,i1,i2,6,c,d)

      do j=1,4
c---sudakov in terms of original variables
      p(i7,j)=alpha*p(i2,j)+beta*p(i1,j)+cos(phi)*c(j)+sin(phi)*d(j)
      p(i1,j)=p(i1,j)/x
      p(i2,j)=p(i2,j)-p(i7,j)+omx*p(i1,j)
      enddo

c---we have now finished generating momenta and can return
      return
      end

