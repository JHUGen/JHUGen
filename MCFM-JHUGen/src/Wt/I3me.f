      double complex function I3me(k1sq,k2sq,k3sq)
      implicit none
      include 'constants.f'
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

      integer i,j,k
      double precision lambda,alpha,y0(0:2),xp(0:2),xm(0:2),
     . yip(0:2),yim(0:2),alphai(0:2),x,y,z,k1sq,k2sq,k3sq
      double precision psq(0:2,0:2),msq(0:2),I3mer,I3mei

      double precision ddilog,polylog
      double precision aa,bb

      integer,parameter::jj(0:2)=(/1,2,0/)
      integer,parameter::kk(0:2)=(/2,0,1/)

      polylog(y,z)=ddilog((y-1d0)/z)-ddilog(y/z)
      lambda(x,y,z)=sqrt(x**2+y**2+z**2-2d0*(x*y+y*z+z*x))

      psq(0,1)=k1sq
      psq(1,2)=k2sq
      psq(2,0)=k3sq
      msq(0)=k1sq
      msq(1)=0d0
      msq(2)=0d0
      alpha=psq(0,1)**2+psq(1,2)**2+psq(2,0)**2
     . -2d0*(psq(0,1)*psq(1,2)+psq(1,2)*psq(2,0)+psq(2,0)*psq(0,1))
      if (alpha .lt. 0d0) then
      write(6,*) 'I3me:Formula not implemented for alpha imaginary'
      stop
      endif
      alpha=sqrt(alpha)
      I3mer=0d0
      I3mei=0d0
      do i=0,2
      j=jj(i)
      k=kk(i)
      alphai(i)=lambda(psq(j,k),msq(j),msq(k))
      y0(i)=0.5d0/alpha/psq(j,k)
     . *(psq(j,k)*(psq(j,k)-psq(k,i)-psq(i,j)
     . +2d0*msq(i)-msq(j)-msq(k))
     . -(psq(k,i)-psq(i,j))*(msq(j)-msq(k))
     . +alpha*(psq(j,k)-msq(j)+msq(k)))
      xp(i)=0.5d0/psq(j,k)*(psq(j,k)-msq(j)+msq(k)+alphai(i))
      xm(i)=0.5d0/psq(j,k)*(psq(j,k)-msq(j)+msq(k)-alphai(i))
      yip(i)=y0(i)-xp(i)
      yim(i)=y0(i)-xm(i)
      I3mer=I3mer+polylog(y0(i),yip(i))+polylog(y0(i),yim(i))

      enddo 

      aa=(y0(1)-1d0)/yim(1)
      bb=y0(0)/yip(0)
      I3mei=pi*dlog(aa/bb)

C---- Minus sign for (-1)^n
      I3me=-Dcmplx(I3mer,I3mei)/alpha

      return
      end

