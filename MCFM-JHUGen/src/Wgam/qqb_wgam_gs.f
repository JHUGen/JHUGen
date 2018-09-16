      subroutine qqb_wgam_gs(p,msq)
      implicit none
      include 'types.f'
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  W^+(nu(p3)+e+(p4))+a(p5)+g(p6)
c     q(-p1)+qbar(-p2) -->  W^-(e-(p3)+nubar(p4))+a(p5)+g(p6)
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'frag.f'
      include 'ewcharge.f'
      include 'phot_dip.f'
      integer:: j,k,nd
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp):: msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     & sub16_2(4),sub26_1(4),dummyv(-nf:nf,-nf:nf),dsubv
      real(dp):: msq56_1(-nf:nf,-nf:nf),msq56_2(-nf:nf,-nf:nf)
      real(dp):: sub56_1,sub56_2
      external qqb_wgam,donothing_gvec
      external qqb_w_g

c--- this flag controls whether or not to include dipoles that
c--- subtract the photon-quark divergence   
      if (frag) then 
        ndmax=3              
      else
        ndmax=2
      endif

      do j=1,mxpart
         phot_dip(j)=.false.
      enddo      

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,6,2,sub16_2,dsubv,msq16_2,dummyv,
     & qqb_wgam,donothing_gvec)
      call dips(2,p,2,6,1,sub26_1,dsubv,msq26_1,dummyv,
     & qqb_wgam,donothing_gvec)

c--- additional fragmentation dipoles, if necessary      
      if (frag) then 
      call dipsfrag(3,p,5,6,2,sub56_2,msq56_2,qqb_w_g)
      phot_dip(3)=.true.
      endif
    
      do j=-nf,nf
      do k=-nf,nf

      do nd=1,ndmax
      msq(nd,j,k)=0._dp
      enddo

      if  ((j == 0) .and. (k == 0)) then
         goto 20

      elseif  ((j > 0) .and. (k < 0)
     &     .or.(j < 0) .and. (k > 0)) then
         msq(1,j,k)=2._dp*cf*sub16_2(qq)*msq16_2(j,k)
         msq(2,j,k)=2._dp*cf*sub26_1(qq)*msq26_1(j,k)        

      elseif ((j .ne. 0) .and. (k == 0)) then
         msq(2,j,k)=2._dp*tr*sub26_1(qg)*(msq26_1(j,+1)+msq26_1(j,+2)
     &   +msq26_1(j,+3)+msq26_1(j,+4)+msq26_1(j,+5)
     &                                 +msq26_1(j,-1)+msq26_1(j,-2)
     &   +msq26_1(j,-3)+msq26_1(j,-4)+msq26_1(j,-5))

c--- additional fragmentation dipole         
         if (frag) then
           if(mod(abs(j),2) == 1) then 
             msq(3,j,k)=Q(2)**2*msq56_2(j,k)*sub56_2
           else
             msq(3,j,k)=Q(1)**2*msq56_2(j,k)*sub56_2
           endif            
         endif

      elseif ((j == 0) .and. (k .ne. 0)) then
         msq(1,j,k)=2._dp*tr*sub16_2(qg)*(msq16_2(+1,k)+msq16_2(+2,k)
     &   +msq16_2(+3,k)+msq16_2(+4,k)+msq16_2(+5,k)
     &                                 +msq16_2(-1,k)+msq16_2(-2,k)
     &   +msq16_2(-3,k)+msq16_2(-4,k)+msq16_2(-5,k))
         
c--- additional fragmentation dipole         
         if (frag) then 
           if(mod(abs(k),2) == 1) then  
             msq(3,j,k)=Q(2)**2*msq56_2(j,k)*sub56_2              
           else
             msq(3,j,k)=Q(1)**2*msq56_2(j,k)*sub56_2              
           endif
         endif

      endif      
 20   continue
      
      enddo
      enddo

      return      
      end


