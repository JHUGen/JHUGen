      subroutine genii(nperms,p,wt,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'debug.f'
      include 'impsample.f'
      include 'npart.f'
      include 'x1x2.f'
      logical:: justjac
      integer:: i1,i2,j,k,nperms
      real(dp):: p(mxpart,4),x,dot,q(mxpart,4),alpha,
     & msq(-nf:nf,-nf:nf)
      real(dp):: omx,Pqq,Pqg,facq,facg,s13,omxmin,a,oma,jacbit
      real(dp):: wt4,wt,wt5_4,wt0
      parameter(wt0=1._dp/eight/pisq)
      common/justjac/justjac

      integer,parameter:: j1(2)=(/1,2/)
      integer,parameter:: j2(2)=(/2,1/)

      i1=j1(nperms)
      i2=j2(nperms)

c first of all calculate the variables with which one started
c---NB all incoming

      s13=2._dp*dot(p,i1,3)
      x=(dot(p,i1,i2)+dot(p,i1,3)+dot(p,i2,3))/dot(p,i1,i2)
c      write(6,*) 'impsample',impsample
c      omxmin=one-xmin
      omxmin=one-xx(i1)
      alpha=-dot(p,i2,3)/dot(p,i1,i2)
      omx=one-x
      a=alpha/omx
      oma=1._dp-a
      if (impsample) then
      jacbit=four*sqrt(omx*omxmin)/(half/sqrt(a)+half/sqrt(oma))
      else
      jacbit=omxmin
      endif

c---at this stage the p are momenta including radiation
      wt5_4=wt0*dot(p,i1,i2)*omx/x*jacbit

      call itransform(p,q,x,i1,3,i2)
      call wtgen(npart,q,wt4)
c---calculate total weight
      wt=wt5_4*wt4


c      write(6,*) 'wt2 in genii.f',wt4
c      write(6,*) 'wt3_2 in genii.f',wt5_4
c      write(6,*) 'wt in genii.f',wt

      if (debug) then
      write(6,*) 'jacbit in genii',jacbit
      write(6,*) 'omxmin in genii',omxmin
      write(6,*) 'omx in genii',omx
      write(6,*) 'wt5_4 in genii',wt5_4
      write(6,*) 'wt4 in genii',wt4
      write(6,*) 'wt in genii',wt
      endif

c---q are in Born level four momenta

      if (justjac) return

      call qqb_WH(q,msq)

      Pqq=CF*(one+x**2)/omx
      Pqg=TR*(one-two*x*omx)
      facq=-2*gsq/x*Pqq
      facg=-2*gsq/x*Pqg

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp


      if     ((j > 0) .and. (k < 0)) then
      msq(j,k)=facq/s13*msq(j,k)
      elseif ((j < 0) .and. (k > 0)) then
      msq(j,k)=facq/s13*msq(j,k)
      elseif ((j == 0) .and. (k > 0)) then
      msq(j,k)=facg/s13*
     &(msq(-1,k)+msq(-2,k)+msq(-3,k)+msq(-4,k)+msq(-5,k))
      elseif ((j == 0) .and. (k < 0)) then
      msq(j,k)=facg/s13*
     &(msq(+1,k)+msq(+2,k)+msq(+3,k)+msq(+4,k)+msq(+5,k))

      endif

      enddo
      enddo

      return
      end
