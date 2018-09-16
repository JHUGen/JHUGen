      subroutine singcoeffg(p,j1,j2,j3,j4,cab,cba,c00ab,c00ba)          
      implicit none
      include 'types.f'
                                                           
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'                                             
      include 'masses.f'                                                
      include 'scale.f'                                                 
      real(dp):: p(mxpart,4),v12,htheta,x,s12,s13,s14,s23,s24,s34,
     & xlog12,xl13,xl14,xl23,xl24,xl34,xlm13,xlm14,xlm23,xlm24,xlm34    
 
      complex(dp):: cab(-2:-1),cba(-2:-1),c00ab(-2:-1),c00ba(-2:-1)    
      integer:: j1,j2,j3,j4                                               
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)     
      htheta(x)=half+half*sign(one,x)                                   
                                                                        
      s12=2d0
     &*(p(j1,4)*p(j2,4)-p(j1,1)*p(j2,1)-p(j1,2)*p(j2,2)-p(j1,3)*p(j2,3))
      s13=2d0
     &*(p(j1,4)*p(j3,4)-p(j1,1)*p(j3,1)-p(j1,2)*p(j3,2)-p(j1,3)*p(j3,3))
      s14=2d0
     &*(p(j1,4)*p(j4,4)-p(j1,1)*p(j4,1)-p(j1,2)*p(j4,2)-p(j1,3)*p(j4,3))
      s23=2d0
     &*(p(j2,4)*p(j3,4)-p(j2,1)*p(j3,1)-p(j2,2)*p(j3,2)-p(j2,3)*p(j3,3))
      s24=2d0
     &*(p(j2,4)*p(j4,4)-p(j2,1)*p(j4,1)-p(j2,2)*p(j4,2)-p(j2,3)*p(j4,3))
      s34=2d0
     &*(p(j3,4)*p(j4,4)-p(j3,1)*p(j4,1)-p(j3,2)*p(j4,2)-p(j3,3)*p(j4,3))
                                                                        
      v12=sqrt(1d0-4d0*mt**4/s12**2)                                    
      xlog12=log((1d0-v12)/(1d0+v12))/v12                               
 
      xlm23=log(mt**2/abs(s23))                                         
      xlm14=log(mt**2/abs(s14))*xn                                      
                                                                        
      xl14=log(musq/abs(s14))                                           
      xl13=log(musq/abs(s13))                                           
      xl23=log(musq/abs(s23))                                           
      xl24=log(musq/abs(s24))                                           
      xl34=log(musq/abs(s34))                                           
                                                                        
      cab(-2)= - 2.D0*xn

      cba(-2)= - 2.D0*xn

      c00ab(-2)= 0

      c00ba(-2)= 0

      cab(-1)= + xn**(-1) + 1.D0/2.D0*xn**(-1)*xlog12 + 1.D0/4.D0*xlm23
     &     + 1.D0/4.D0*xlm14 - 1.D0/2.D0*xl34 + 1.D0/4.D0*xl23 + 1.D0/4.
     &    D0*xl14 - 1.D0/4.D0*xlog12 - xn - 1.D0/2.D0*xn*xlm23 - 1.D0/2.
     &    D0*xn*xlm14 - xn*xl34 - 1.D0/2.D0*xn*xl23 - 1.D0/2.D0*xn*xl14
     &     + 2.D0*htheta(s12)*impi*v12**(-1) + 2.D0*htheta(s13)*impi + 
     &    2.D0*htheta(s14)*impi + 2.D0*htheta(s23)*impi + 2.D0*htheta(
     &    s24)*impi + 2.D0*htheta(s34)*impi

      cba(-1)= + xn**(-1) + 1.D0/2.D0*xn**(-1)*xlog12 + 1.D0/4.D0*xlm23
     &     + 1.D0/4.D0*xlm14 - 1.D0/2.D0*xl34 + 1.D0/4.D0*xl23 + 1.D0/4.
     &    D0*xl14 - 1.D0/4.D0*xlog12 - xn - 1.D0/2.D0*xn*xlm24 - 1.D0/2.
     &    D0*xn*xlm13 - xn*xl34 - 1.D0/2.D0*xn*xl24 - 1.D0/2.D0*xn*xl13
     &     + 2.D0*htheta(s12)*impi*v12**(-1) + 2.D0*htheta(s13)*impi + 
     &    2.D0*htheta(s14)*impi + 2.D0*htheta(s23)*impi + 2.D0*htheta(
     &    s24)*impi + 2.D0*htheta(s34)*impi

      c00ab(-1)= + 1.D0/4.D0*xlm24 + 1.D0/4.D0*xlm13 - 1.D0/2.D0*xl34
     &     + 1.D0/4.D0*xl24 + 1.D0/4.D0*xl13 - 1.D0/4.D0*xlog12

      c00ba(-1)= + 1.D0/4.D0*xlm23 + 1.D0/4.D0*xlm14 - 1.D0/2.D0*xl34
     &     + 1.D0/4.D0*xl23 + 1.D0/4.D0*xl14 - 1.D0/4.D0*xlog12


      return                                                            
      end                                                               
