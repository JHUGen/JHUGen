      subroutine qqb_gmgmjt(p,msq) 
c--- matrix element squared for the process
c---    q(p1) + q~(p2) --> gam(p3) + gam(p4) + g(p5)
c--- (and all crossings)
c---
c--- C. Williams, March 2013
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'ewcharge.f' 
      include 'ewcouple.f'
      include 'qcdcouple.f' 
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf)
      double complex qqbg(2,2,2,2),qbqg(2,2,2,2)    
      double complex qgqb(2,2,2,2),qbgq(2,2,2,2)
      double complex gqqb(2,2,2,2),gqbq(2,2,2,2)
      double precision qqbg_sum,qbqg_sum
      double precision qgqb_sum,qbgq_sum
      double precision gqbq_sum,gqqb_sum
      integer h1,h2,h3,h4,j,k 
      double precision fac,statfac
      parameter(statfac=0.5d0)

      qqbg_sum=0d0 
      qgqb_sum=0d0 
      gqqb_sum=0d0 
c      qbqg_sum=0d0 
c      qbgq_sum=0d0 
c      gqbq_sum=0d0 
      msq(:,:)=0d0 

      fac=8d0*cf*xn*gsq*esq**2*statfac

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

        qqbg_sum=qqbg_sum+cdabs(qqbg(h1,h2,h3,h4))**2
        qgqb_sum=qgqb_sum+cdabs(qgqb(h1,h2,h3,h4))**2
        gqqb_sum=gqqb_sum+cdabs(gqqb(h1,h2,h3,h4))**2
c        qbqg_sum=qbqg_sum+cdabs(qbqg(h1,h2,h3,h4))**2
c        qbgq_sum=qbgq_sum+cdabs(qbgq(h1,h2,h3,h4))**2
c        gqbq_sum=gqbq_sum+cdabs(gqbq(h1,h2,h3,h4))**2
            
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
      include 'constants.f' 
      include 'zprods_decl.f' 
      integer i1,i2,i3,i4,i5
      double complex amp(2,2,2,2) 
      double complex amp_2gam1g
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

