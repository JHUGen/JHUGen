      function I3me(k1sq,k2sq,k3sq)
      implicit none
      include 'types.f'
      complex(dp):: I3me

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
c     Triangle function with one internal mass and two spacelike
c     offshell lines and one timelike of mass k1sq
c
c                         | k2
c                         |
c                         /\
c                        /__\____
c               k3  *****-------- k1,(k1sq=msq)
c
c   Mass of the one massive internal line is equal to k1sq
c
c   Adapted from
c   %\cite{Denner:kt}
c   \bibitem{Denner:kt}
c   A.~Denner,
c   %``Techniques For Calculation Of Electroweak Radiative
c   Corrections At The One Loop Level And Results For W Physics At Lep-200,''
c   Fortsch.\ Phys.\  {\bf 41}, 307 (1993).
c   %%CITATION = FPYKA,41,307;%%

      integer:: i,j,k
      real(dp):: lambda,alpha,y0(0:2),xp(0:2),xm(0:2),
     & yip(0:2),yim(0:2),alphai(0:2),x,y,z,k1sq,k2sq,k3sq
      real(dp):: psq(0:2,0:2),msq(0:2),I3mer,I3mei

      real(dp):: ddilog,polylog
      real(dp):: aa,bb

      integer,parameter::jj(0:2)=(/1,2,0/)
      integer,parameter::kk(0:2)=(/2,0,1/)

      polylog(y,z)=ddilog((y-1._dp)/z)-ddilog(y/z)
      lambda(x,y,z)=sqrt(x**2+y**2+z**2-2._dp*(x*y+y*z+z*x))

      psq(0,1)=k1sq
      psq(1,2)=k2sq
      psq(2,0)=k3sq
      msq(0)=k1sq
      msq(1)=0._dp
      msq(2)=0._dp
      alpha=psq(0,1)**2+psq(1,2)**2+psq(2,0)**2
     & -2._dp*(psq(0,1)*psq(1,2)+psq(1,2)*psq(2,0)+psq(2,0)*psq(0,1))
      if (alpha < 0._dp) then
      write(6,*) 'I3me:Formula not implemented for alpha imaginary'
      stop
      endif
      alpha=sqrt(alpha)
      I3mer=0._dp
      I3mei=0._dp
      do i=0,2
      j=jj(i)
      k=kk(i)
      alphai(i)=lambda(psq(j,k),msq(j),msq(k))
      y0(i)=0.5_dp/alpha/psq(j,k)
     & *(psq(j,k)*(psq(j,k)-psq(k,i)-psq(i,j)
     & +2._dp*msq(i)-msq(j)-msq(k))
     & -(psq(k,i)-psq(i,j))*(msq(j)-msq(k))
     & +alpha*(psq(j,k)-msq(j)+msq(k)))
      xp(i)=0.5_dp/psq(j,k)*(psq(j,k)-msq(j)+msq(k)+alphai(i))
      xm(i)=0.5_dp/psq(j,k)*(psq(j,k)-msq(j)+msq(k)-alphai(i))
      yip(i)=y0(i)-xp(i)
      yim(i)=y0(i)-xm(i)
      I3mer=I3mer+polylog(y0(i),yip(i))+polylog(y0(i),yim(i))

      enddo

      aa=(y0(1)-1._dp)/yim(1)
      bb=y0(0)/yip(0)
      I3mei=pi*log(aa/bb)

C---- Minus sign for (-1)^n
      I3me=-cplx2(I3mer,I3mei)/alpha

      return
      end

