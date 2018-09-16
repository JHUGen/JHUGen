!----- fills amplitudes for gg=>HZ process
!===== these are the pieces which come from box topologies 
      subroutine gg_HZ_box(p,amp,mt2) 
      implicit none 
      include 'types.f'
      include 'constants.f' 
      include 'mxpart.f'
      include 'masses.f' 
      include 'zprods_decl.f' 
      include 'sprods_com.f' 
      include 'scale.f'
      include 'cplx.h'
      real(kind=dp):: p(mxpart,4),mt2
      complex(kind=dp):: amp(2,2,2),test
      complex(kind=dp):: ggHZ_pp_box
      complex(kind=dp):: ggHZ_mp_box
      external ggHZ_pp_box 
      external ggHZ_mp_box 
      complex(kind=dp):: prop_34 



!----- debug
!----- KC phase space point check 
!      p(1,4)=  -3.0000000000000000_dp
!      p(1,1)=   2.1213203435596424_dp
!      p(1,2)=   1.0606601717798212_dp
!      p(1,3)=   1.8371173070873839_dp
!      p(2,4)=  -3.0000000000000000_dp
!      p(2,1)=  -1.3016512173526746_dp
!      p(2,2)=   2.6346884688170253_dp
!      p(2,3)= -0.60342421284441206_dp
!      p(3,4)=  0.85714285714285710_dp
!      p(3,1)= -0.31578947368421051_dp
!      p(3,2)=  0.79685060448070799_dp
!      p(3,3)=   0.0000000000000000_dp
!      p(4,4)=  -1.5000000000000000_dp
!      p(4,1)= -0.88167787843870971_dp
!      p(4,2)=   1.1895840590472970_dp
!      p(4,3)= -0.23986222114450001_dp


!      scale=one 
!      musq=one
!      mt=0.4255266775_dp
!----- end debug

      call spinoru(4,p,za,zb)
!===== z propogator 
      prop_34=s(3,4)/cplx2(s(3,4)-zmass**2,zmass*zwidth) 

!------- fill amplitudes which come from box diagrams

!------ left handed lepton coupling (- sign from line reversal)
      amp(2,2,1)=ggHZ_pp_box(1,2,3,4,za,zb,mt2)
      amp(1,2,1)=ggHZ_mp_box(1,2,3,4,za,zb,mt2)
      amp(2,1,1)=-ggHZ_mp_box(1,2,4,3,zb,za,mt2)
      amp(1,1,1)=-ggHZ_pp_box(1,2,4,3,zb,za,mt2)
      
!------- right handed lepton coupling 
      amp(2,2,2)=ggHZ_pp_box(1,2,4,3,za,zb,mt2)
      amp(1,2,2)=ggHZ_mp_box(1,2,4,3,za,zb,mt2)
      amp(2,1,2)=-ggHZ_mp_box(1,2,3,4,zb,za,mt2)
      amp(1,1,2)=-ggHZ_pp_box(1,2,3,4,zb,za,mt2)


!----- debug, all checked against numeric D-dim code CW 2/10/15
!      write(6,*) '******** L ********* '
!      write(6,*) '2,2,1',amp(2,2,1)*im
!      write(6,*) '1,2,1',amp(1,2,1)*im
!      write(6,*) '2,1,1',amp(2,1,1)*im
!      write(6,*) '1,1,1',amp(1,1,1)*im
!      write(6,*) '******************** '
!      write(6,*) 
!      write(6,*) '******** R ********* '
!      write(6,*) '2,2,2',amp(2,2,2)*im
!      write(6,*) '1,2,2',amp(1,2,2)*im
!      write(6,*) '2,1,2',amp(2,1,2)*im
!      write(6,*) '1,1,2',amp(1,1,2)*im
!      write(6,*) '******************** '
!      pause

!---- dress with propagators  
      amp(:,:,:)=prop_34*amp(:,:,:)
      
      return 
      end

