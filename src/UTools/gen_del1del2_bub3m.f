
      subroutine gen_del1del2_bub3m(p,n,P1,P2,za,zb) 
      implicit none 
      include 'constants.f' 
      include 'zprods_decl.f' 
      include 'utoolscc.f' 
      double precision p(mxpart,4),P1(4),P2(4) 
      integer n,nu
      double precision g1,DEL,a,b 
      double precision P1dP2,P1dP1,P2dP2
      double precision del1(4),del2(4) 
      double precision pmo(mxpart,4) 

      pmo(:,:)=p(:,:) 

      P1dP2=0d0 
      P1dP1=0d0 
      P2dP2=0d0
      do nu=1,3 
         P1dP2=P1dP2-P1(nu)*P2(nu) 
         P1dP1=P1dP1-P1(nu)*P1(nu) 
         P2dP2=P2dP2-P2(nu)*P2(nu)
      enddo
      nu=4
      P1dP2=P1dP2+P1(nu)*P2(nu)      
      P1dP1=P1dP1+P1(nu)*P1(nu)     
      P2dP2=P2dP2+P2(nu)*P2(nu) 
      
      DEL=P1dP2**2-P1dP1*P2dP2 
      g1=P1dP2+dsqrt(DEL) 
      a=P1dP1/g1
      b=P2dP2/g1
      
      do nu=1,4
         del1(nu)=one/(one-a*b)*(P1(nu)-a*P2(nu))
         del2(nu)=one/(one-a*b)*(P2(nu)-b*P1(nu))
      enddo
      
      do nu=1,4
         pmo(n+1,nu)=del1(nu)
         pmo(n+2,nu)=del2(nu)
      enddo

      if (utoolscc) then
        call spinoruorz(n+2,pmo,zb,za) 
      else      
        call spinoruorz(n+2,pmo,za,zb) 
      endif
      
      return 
      end 
