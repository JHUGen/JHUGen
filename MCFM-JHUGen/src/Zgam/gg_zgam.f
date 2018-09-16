      subroutine gg_zgam(p,msq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'

      real(dp):: p(mxpart,4),msq
      complex(dp):: ggZgam_vec,prop
!     arg 1 = gluon mom p1, arg 2 = gluon mom p2, arg 3 = photon mom p5, arg 4 = helicity of lepton line 34,
      complex(dp):: amp(2,2,2,2),cu(2),cd(2)
      integer:: h1,h2,h3,h4
      real(dp):: fac

      call spinoru(5,p,za,zb)
      amp(2,2,2,2)=+ggZgam_vec('+++',p,1,2,5,4,3,za,zb)
      amp(1,1,1,1)=-ggZgam_vec('+++',p,1,2,5,4,3,zb,za)
      amp(2,1,1,2)=+ggZgam_vec('+--',p,1,2,5,4,3,za,zb)
      amp(1,2,2,1)=-ggZgam_vec('+--',p,1,2,5,4,3,zb,za)
      
      amp(2,2,2,1)=+ggZgam_vec('+++',p,1,2,5,3,4,za,zb)
      amp(1,1,1,2)=-ggZgam_vec('+++',p,1,2,5,3,4,zb,za)
      amp(2,1,1,1)=+ggZgam_vec('+--',p,1,2,5,3,4,za,zb)
      amp(1,2,2,2)=-ggZgam_vec('+--',p,1,2,5,3,4,zb,za)   

      amp(2,2,1,2)=+ggZgam_vec('++-',p,1,2,5,4,3,za,zb)
      amp(1,1,2,1)=-ggZgam_vec('++-',p,1,2,5,4,3,zb,za)
      amp(2,2,1,1)=+ggZgam_vec('++-',p,1,2,5,3,4,za,zb)
      amp(1,1,2,2)=-ggZgam_vec('++-',p,1,2,5,3,4,zb,za)

      amp(1,2,1,1)=+ggZgam_vec('+--',p,2,1,5,3,4,za,zb)
      amp(2,1,2,2)=-ggZgam_vec('+--',p,2,1,5,3,4,zb,za)
      amp(2,1,2,1)=-ggZgam_vec('+--',p,2,1,5,4,3,zb,za)
      amp(1,2,1,2)=+ggZgam_vec('+--',p,2,1,5,4,3,za,zb)

      prop=s(3,4)/cplx2(s(3,4)-zmass**2,zmass*zwidth) 
      
!     argument is helicity of lepton line
      cu(1)=Qu*(Qu*q1+0.5_dp*(L(2)+R(2))*l1*prop)
      cd(1)=Qd*(Qd*q1+0.5_dp*(L(1)+R(1))*l1*prop)
      
      cu(2)=Qu*(Qu*q1+0.5_dp*(L(2)+R(2))*r1*prop)
      cd(2)=Qd*(Qd*q1+0.5_dp*(L(1)+R(1))*r1*prop)
      
    
      msq=0._dp

      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2

                  msq=msq+abs(2._dp*amp(h1,h2,h3,h4)*cu(h4)
     &                         +3._dp*amp(h1,h2,h3,h4)*cd(h4))**2
                  
               enddo
            enddo
         enddo
      enddo      

!     colour prefactor
      fac=avegg*V*(2._dp*rt2*esq*gsq/(16._dp*pisq))**2*esq      
      msq=msq*fac

      return 
      end
      
      

      
      
