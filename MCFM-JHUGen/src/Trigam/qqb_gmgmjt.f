      subroutine qqb_gmgmjt(p,msq)
      implicit none
      include 'types.f'
c--- matrix element squared for the process
c---    q(p1) + q~(p2) --> gam(p3) + gam(p4) + g(p5)
c--- (and all crossings)
c---
c--- C. Williams, March 2013

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      complex(dp):: qqbg(2,2,2,2),qbqg(2,2,2,2)
      complex(dp):: qgqb(2,2,2,2),qbgq(2,2,2,2)
      complex(dp):: gqqb(2,2,2,2),gqbq(2,2,2,2)
      real(dp):: qqbg_sum,qbqg_sum
      real(dp):: qgqb_sum,qbgq_sum
      real(dp):: gqbq_sum,gqqb_sum
      integer:: h1,h2,h3,h4,j,k
      real(dp):: fac,statfac
      parameter(statfac=0.5_dp)

      qqbg_sum=0._dp
      qgqb_sum=0._dp
      gqqb_sum=0._dp
c      qbqg_sum=0._dp
c      qbgq_sum=0._dp
c      gqbq_sum=0._dp
      msq(:,:)=0._dp

      fac=8._dp*cf*xn*gsq*esq**2*statfac

      call spinoru(5,p,za,zb)
      call amp_lord_gmgmjt(1,2,5,3,4,za,zb,qqbg)
      call amp_lord_gmgmjt(1,5,2,3,4,za,zb,qgqb)
      call amp_lord_gmgmjt(2,5,1,3,4,za,zb,gqqb)
c      call amp_lord_gmgmjt(2,1,5,3,4,za,zb,qbqg)
c      call amp_lord_gmgmjt(5,1,2,3,4,za,zb,qbgq)
c      call amp_lord_gmgmjt(5,2,1,3,4,za,zb,gqbq)

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2

        qqbg_sum=qqbg_sum+abs(qqbg(h1,h2,h3,h4))**2
        qgqb_sum=qgqb_sum+abs(qgqb(h1,h2,h3,h4))**2
        gqqb_sum=gqqb_sum+abs(gqqb(h1,h2,h3,h4))**2
c        qbqg_sum=qbqg_sum+abs(qbqg(h1,h2,h3,h4))**2
c        qbgq_sum=qbgq_sum+abs(qbgq(h1,h2,h3,h4))**2
c        gqbq_sum=gqbq_sum+abs(gqbq(h1,h2,h3,h4))**2

      enddo
      enddo
      enddo
      enddo

c--- use symmetry to avoid calculating half the matrix elements
      qbqg_sum=qqbg_sum
      gqbq_sum=gqqb_sum
      qbgq_sum=qgqb_sum

      do j=1,nf
        msq(j,-j)=fac*aveqq*qqbg_sum*Q(j)**4
        msq(-j,j)=fac*aveqq*qbqg_sum*Q(j)**4
        msq(0,+j)=fac*aveqg*gqqb_sum*Q(j)**4
        msq(+j,0)=fac*aveqg*qgqb_sum*Q(j)**4
        msq(-j,0)=fac*aveqg*qbgq_sum*Q(j)**4
        msq(0,-j)=fac*aveqg*gqbq_sum*Q(j)**4
      enddo

      return
      end


      subroutine amp_lord_gmgmjt(i1,i2,i3,i4,i5,za,zb,amp)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5
      complex(dp):: amp(2,2,2,2)
      complex(dp):: amp_2gam1g
!===== default is for gluon MHV amplitude

      amp(:,:,:,:)=czip
      amp(1,1,1,1)=czip
      amp(2,2,2,2)=czip
      amp(2,1,1,1)=czip
      amp(1,2,2,2)=czip

!===== 3,4,5 symmetry
      amp(1,1,2,2)=amp_2gam1g(i1,i2,i3,i4,i5,za,zb)
      amp(1,2,1,2)=amp_2gam1g(i1,i2,i4,i3,i5,za,zb)
      amp(1,2,2,1)=amp_2gam1g(i1,i2,i5,i3,i4,za,zb)

!===== line reversal
      amp(2,1,2,2)=amp_2gam1g(i2,i1,i3,i4,i5,za,zb)
      amp(2,2,1,2)=amp_2gam1g(i2,i1,i4,i3,i5,za,zb)
      amp(2,2,2,1)=amp_2gam1g(i2,i1,i5,i3,i4,za,zb)

!====== conjugation
      amp(2,2,1,1)=-amp_2gam1g(i1,i2,i3,i4,i5,zb,za)
      amp(2,1,2,1)=-amp_2gam1g(i1,i2,i4,i3,i5,zb,za)
      amp(2,1,1,2)=-amp_2gam1g(i1,i2,i5,i3,i4,zb,za)
!====== conjugation + line reversal
      amp(1,2,1,1)=-amp_2gam1g(i2,i1,i3,i4,i5,zb,za)
      amp(1,1,2,1)=-amp_2gam1g(i2,i1,i4,i3,i5,zb,za)
      amp(1,1,1,2)=-amp_2gam1g(i2,i1,i5,i3,i4,zb,za)

      return
      end

