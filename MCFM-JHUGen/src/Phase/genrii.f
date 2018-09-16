      subroutine genrii(p,i1,i2,i7,r1,r2,phit,wt5_4,*)
      implicit none
      include 'types.f'
c----i1,i2 initial state vectors.
c----i7 label of generated vector
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'debug.f'
      include 'impsample.f'
      include 'x1x2.f'
      integer:: i1,i2,i7,j,k
      real(dp):: p(mxpart,4),rtalbe,c(4),d(4),phi,phit,jacbit
      real(dp):: qtDp(mxpart),qsDp(mxpart),dot,r1,r2
      real(dp):: q(4),qs(4),qt(4),beta,alpha,qsDqs,qtDqt
      real(dp):: a,oma,x,omx,omxmin,wt0,wt5_4
      parameter(wt0=1._dp/8._dp/pisq)

      phi=twopi*phit
c      omxmin=one-xmin
      omxmin=one-xx(i1)
      r1=0.5_dp
c      write(6,*) 'Enter r1'
c      read(5,*) r1
      write(6,*) 'Enter r2'
      read(5,*) r2
      if (impsample) then
        omx=omxmin*r1**2
        x=one-omx
        if (r2 < 0.5_dp) oma=(two*r2)**2
        if (r2 >= 0.5_dp) oma=one-(two*r2-one)**2
        a=one-oma
        jacbit=four*sqrt(omx*omxmin)/(half/sqrt(oma)+half/sqrt(a))
      else
        omx=omxmin*r1
        x=one-omx
        oma=r2
        a=one-oma
        jacbit=omxmin
      endif

      alpha=omx*a
      beta=omx-alpha
      rtalbe=sqrt(beta*alpha)

      wt5_4=wt0*dot(p,i1,i2)/x**2*omx*jacbit

c      write(6,*) 'wt5_4 in genrii',wt5_4

      if (debug) write(6,*) 'jacbit in genrii',jacbit
      if (debug) write(6,*) 'wt0 in genrii',wt0
      if (debug) write(6,*) 'omxmin in genrii',omxmin
      if (debug) write(6,*) 'omx in genrii',omx
      if (debug) write(6,*) 'x in genrii',x
      if (debug) write(6,*) 'i1 in genrii',i1
      if (debug) write(6,*) 'i2 in genrii',i2
      if (debug) write(6,*) 'omxmin in genrii',omxmin
      if (debug) write(6,*) 'wt5_4 in genrii',wt5_4
      if (debug) write(6,*)

c---rescale p(i1)
      do j=1,4
      p(i1,j)=p(i1,j)/x
      enddo

c---Sudakov wrt new vectors
c---generate transverse vectors c and d with length^2=rtalbe^2*2*p1Dp2
c-- with direction in transverse plane picked by 3
      call gtperp(rtalbe,p,i1,i2,3,c,d)

c---generate p7 and auxiliary vectors
      do j=1,4
      p(i7,j)=alpha*p(i1,j)+beta*p(i2,j)+cos(phi)*c(j)+sin(phi)*d(j)
      q(j) =p(i1,j)+p(i2,j)-p(i7,j)
      qt(j)=x*p(i1,j)+p(i2,j)
      qs(j)=q(j)+qt(j)
      enddo

      qtDqt=qt(4)**2-qt(1)**2-qt(2)**2-qt(3)**2
      qsDqs=qs(4)**2-qs(1)**2-qs(2)**2-qs(3)**2

C--generate the remaining vectors 3 through i7-1
      do k=3,i7-1
      qtDp(k)=qt(4)*p(k,4)-qt(1)*p(k,1)-qt(2)*p(k,2)-qt(3)*p(k,3)
      qsDp(k)=qs(4)*p(k,4)-qs(1)*p(k,1)-qs(2)*p(k,2)-qs(3)*p(k,3)
      do j=1,4
      p(k,j)=p(k,j)+two*(qtDp(k)*q(j)/qtDqt-qsDp(k)*qs(j)/qsDqs)
      enddo
      enddo
c----this completes the generation of the new momenta;
c----we can now return

      return
      end

