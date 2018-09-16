!====== fragmentation routine for fourgam
!==== based of trigam tree C.Williams Oct 2014
      subroutine qqb_fourgam_frag(p,msq) 
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      include 'frag.f' 
      integer:: j
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),
     & ampsq_3gam1g,qqb,qbq,qg,gq,qbg,gqb,fac,cfac
      real(dp):: fsq,D(0:5)
      integer:: i 

      fsq=frag_scale**2
c---- Generate array D(j) corresponding to MCFM notation 0=gluon 1=up 2=down ....
        do i=0,5
           D(i)=0d0
           if     (fragset == 'BFGset_I') then
              call get_frag(z_frag,fsq,1,i,D(i))   
           elseif (fragset == 'BFGsetII') then  
              call get_frag(z_frag,fsq,2,i,D(i))   
           elseif (fragset == 'GdRG__LO') then 
            call GGdR_frag(z_frag,i,D(i),0)
           elseif (fragset == 'GdRG_NLO') then 
            call GGdR_frag(z_frag,i,D(i),1)
           else
              write(6,*) 'Unrecognized fragmentation set name: ',fragset
              stop        
           endif
        enddo
      call spinoru(6,p,za,zb)

      fac=16d0*esq**3*gsq*xn*Cf/3d0
      
      qqb=aveqq*fac*ampsq_3gam1g(6,3,4,5,1,2,za,zb)
      qg=aveqg*fac*ampsq_3gam1g(2,3,4,5,1,6,za,zb)
      gq=aveqg*fac*ampsq_3gam1g(1,3,4,5,2,6,za,zb)
      qbq=qqb
      qbg=qg
      gqb=gq
      
c      qbq=aveqq*fac*ampsq_3gam1g(6,3,4,5,2,1,za,zb)
c      qbg=aveqg*fac*ampsq_3gam1g(2,3,4,5,6,1,za,zb)
c      gqb=aveqg*fac*ampsq_3gam1g(1,3,4,5,6,2,za,zb)

      msq(:,:)=0d0
      do j=1,nf
        cfac=Q(j)**6
        msq(j,-j)=cfac*qqb*D(0)
        msq(-j,j)=cfac*qbq*D(0)
        msq(j,0)=cfac*qg*D(j)
        msq(0,j)=cfac*gq*D(j)
        msq(-j,0)=cfac*qbg*D(j)
        msq(0,-j)=cfac*gqb*D(j)
      enddo
      
      return
      end
      
