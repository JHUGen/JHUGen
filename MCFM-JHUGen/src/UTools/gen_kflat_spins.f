!===== Subroutine for generating a flatted vector for the 
!===== two mass triangle topology, here P1^2 ne 0 and P2^2=0 
      subroutine gen_kflat_2mtri(p,n,P1,P2,za,zb) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f' 
      include 'utoolscc.f' 
      integer:: n,nu 
      real(dp):: p(mxpart,4),P1(4),P2(4) 
      real(dp):: po(mxpart,4) 
      real(dp):: P1sq,P1dP2

      po(:,:)=p(:,:) 
      P1sq=P1(4)**2-P1(3)**2-P1(2)**2-P1(1)**2
      
      P1dP2=2._dp*(P1(4)*P2(4)-P1(3)*P2(3)-P1(2)*P2(2)-P1(1)*P2(1))

      do nu=1,4
         po(n+1,nu)=P1(nu)-P1sq/P1dP2*P2(nu) 
      enddo
      if (utoolscc) then
        call spinoruorz(n+1,po,zb,za) 
      else      
        call spinoruorz(n+1,po,za,zb) 
      endif
      
      return 
      end
         
      

!===== Subroutine for generating a flatted vector for the 
!===== three mass triangle topology, here P1^2 ne 0 and P2^2 ne 0  
!===== also return gamma (with the positive sign)
      subroutine gen_kflat_3mtri(p,n,P1,P2,za,zb,gamma) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f' 
      include 'utoolscc.f' 
      integer:: n,nu 
      real(dp):: p(mxpart,4),P1(4),P2(4) 
      real(dp):: po(mxpart,4) 
      real(dp):: P1sq,P1dP2,P2sq,gamma

      po(:,:)=p(:,:) 
      P1sq=P1(4)**2-P1(3)**2-P1(2)**2-P1(1)**2
      P2sq=P2(4)**2-P2(3)**2-P2(2)**2-P2(1)**2
      
      P1dP2=(P1(4)*P2(4)-P1(3)*P2(3)-P1(2)*P2(2)-P1(1)*P2(1))

      gamma=P1dP2+sqrt(P1dP2**2-P1sq*P2sq)

      do nu=1,4
         po(n+1,nu)=P1(nu)-P1sq/gamma*P2(nu) 
         po(n+1,nu)=po(n+1,nu)/(one-P1sq*P2sq/gamma**2) 
         po(n+2,nu)=P2(nu)-P2sq/gamma*P1(nu) 
         po(n+2,nu)=po(n+2,nu)/(one-P1sq*P2sq/gamma**2) 
      enddo
      if (utoolscc) then
        call spinoruorz(n+2,po,zb,za) 
      else      
        call spinoruorz(n+2,po,za,zb) 
      endif
      
      
      return 
      end
         
      
