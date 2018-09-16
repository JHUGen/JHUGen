      subroutine qqb_wgam_frag(p,msq)
      implicit none
      include 'types.f'
      
c----Matrix element for W production
C----averaged over initial colours and spins
C for nwz=+1
c     u(-p1)+dbar(-p2)--> W^+(n(p3)+e^+(p4))   + (g(p5)-->gamma(z*p5))
C For nwz=-1
c     d(-p1)+ubar(-p2)--> W^-(e^-(p3)+nbar(p4))+ (g(p5)-->gamma(z*p5)) 
c---  qqb_w_g routine modified for fragmentation C. Williams Dec 2010
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ckm.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'frag.f'
      integer:: j,k,i
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      real(dp):: qqbWg,qbqWg,qgWq,qbgWqb,gqbWqb,gqWq,w1jet
      real(dp):: fsq,D(0:5)
      common/D/D
!$omp threadprivate(/D/)

c      fsq=wmass**2+half*((pttwo(3,4,p))**2+z_frag**2*pt(5,p)**2)
      fsq=frag_scale**2
c---- Generate array D(j) corresponding to MCFM notation 0=gluon 1=down 2=up ....
      do i=0,5
         D(i)=0._dp
       if     (fragset == 'BFGset_I') then
            call get_frag(z_frag,fsq,1,i,D(i))   
       elseif (fragset == 'BFGsetII') then  
            call get_frag(z_frag,fsq,2,i,D(i))   
         elseif (fragset == 'GdRG__LO') then 
            call GGdR_frag(z_frag,i,D(i),0) 
         else
            write(6,*) 'Unrecognized fragmentation set name: ',fragset
            stop
         endif
      enddo

     

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call dotem(5,p,s)
c---calculate the propagator
      fac=gwsq**2*gsq*V

      qqbWg= +aveqq*fac*w1jet(1,2,3,4,5)
      gqbWqb=-aveqg*fac*w1jet(5,2,3,4,1)
      qgWq=  -aveqg*fac*w1jet(1,5,3,4,2)
      
      qbqWg= +aveqq*fac*w1jet(2,1,3,4,5)
      qbgWqb=-aveqg*fac*w1jet(5,1,3,4,2)
      gqWq=  -aveqg*fac*w1jet(2,5,3,4,1)

      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=Vsq(j,k)*qqbWg*D(0)
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=Vsq(j,k)*qbqWg*D(0)
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=
     &   (Vsq(j,-1)*D(1)+Vsq(j,-2)*D(2)+Vsq(j,-3)*D(3)
     &    +Vsq(j,-4)*D(4)+Vsq(j,-5)*D(5))*qgWq
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=
     &    (Vsq(j,+1)*D(1)+Vsq(j,+2)*D(2)+Vsq(j,+3)*D(3)
     &    +Vsq(j,+4)*D(4)+Vsq(j,+5)*D(5))*qbgWqb
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=
     &    (Vsq(-1,k)*D(1)+Vsq(-2,k)*D(2)+Vsq(-3,k)*D(3)
     &    +Vsq(-4,k)*D(4)+Vsq(-5,k)*D(5))*gqWq
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=
     &    (Vsq(+1,k)*D(1)+Vsq(+2,k)*D(2)+Vsq(+3,k)*D(3)
     &    +Vsq(+4,k)*D(4)+Vsq(+5,k)*D(5))*gqbWqb
      endif

      enddo
      enddo
      return
      end
 
     

