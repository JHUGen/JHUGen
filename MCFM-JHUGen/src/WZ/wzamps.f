      subroutine wzamps(j1,j2,j3,j4,j5,j6,j7,za,zb,f)
      implicit none
      include 'types.f'
c- first label of fs,ft,fu, is gluon polarization, second is zdecay line

      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'anomcoup.f'
      include 'xanomcoup.f'
      include 'masses.f'
      include 'srdiags.f'
      integer:: j1,j2,j3,j4,j5,j6,j7,mplus,minus,jtype,j,k
      complex(dp):: A7treea,B7treea,B7treeb,A7b_1,A7b_2,A7b_3,A7b_4
      complex(dp):: f(10,2,2)
      real(dp):: s127
      parameter(minus=1,mplus=2)

c----initialize to zero
c--- added amplitude 10 on 4/25/2002, to include anomalous couplings
      do jtype=4,9
      do j=1,2
      do k=1,2
      f(jtype,j,k)=czip
      enddo
      enddo
      enddo
      
      s127=real(za(1,7)*zb(7,1)+za(2,7)*zb(7,2)+za(1,2)*zb(2,1))

c      f(1,mplus,minus)=+A7treeb_anom_wz(j1,j2,j3,j4,j5,j6,j7,za,zb) !fs
      call A7treeb_anom_wz(j1,j2,j3,j4,j5,j6,j7,za,zb,
     & A7b_1,A7b_2,A7b_3,A7b_4)                                     !fs
      f(1,mplus,minus)=
     &           A7b_1*(2._dp+xdelg1_z+xdelk_z+xlambda_z*s127/wmass**2)
     &          +A7b_2*(2._dp+xdelg1_z+xdelk_z+xlambda_z)
     &          +A7b_3*2._dp*(1._dp+xdelg1_z) 
     &          +A7b_4*xlambda_z/wmass**2
      f(10,mplus,minus)=
     &           A7b_1*(2._dp+xdelg1_g+xdelk_g+xlambda_g*s127/wmass**2)
     &          +A7b_2*(2._dp+xdelg1_g+xdelk_g+xlambda_g)
     &          +A7b_3*2._dp*(1._dp+xdelg1_g) 
     &          +A7b_4*xlambda_g/wmass**2
      f(2,mplus,minus)=+A7treea(j1,j2,j3,j4,j5,j6,j7,za,zb)         !ft
      f(3,mplus,minus)=+A7treea(j1,j2,j6,j5,j4,j3,j7,za,zb)         !fu

c      f(1,mplus,mplus)=+A7treeb_anom_wz(j1,j2,j3,j4,j6,j5,j7,za,zb)
      call A7treeb_anom_wz(j1,j2,j3,j4,j6,j5,j7,za,zb,
     & A7b_1,A7b_2,A7b_3,A7b_4)
      f(1,mplus,mplus)=
     &           A7b_1*(2._dp+xdelg1_z+xdelk_z+xlambda_z*s127/wmass**2)
     &          +A7b_2*(2._dp+xdelg1_z+xdelk_z+xlambda_z)
     &          +A7b_3*2._dp*(1._dp+xdelg1_z) 
     &          +A7b_4*xlambda_z/wmass**2
      f(10,mplus,mplus)=
     &           A7b_1*(2._dp+xdelg1_g+xdelk_g+xlambda_g*s127/wmass**2)
     &          +A7b_2*(2._dp+xdelg1_g+xdelk_g+xlambda_g)
     &          +A7b_3*2._dp*(1._dp+xdelg1_g) 
     &          +A7b_4*xlambda_g/wmass**2
      f(2,mplus,mplus)=+A7treea(j1,j2,j3,j4,j6,j5,j7,za,zb)
      f(3,mplus,mplus)=+A7treea(j1,j2,j5,j6,j4,j3,j7,za,zb)

c-- old version
c      f(1,minus,minus)=-A7treeb_anom_wz(j2,j1,j5,j6,j3,j4,j7,zb,za)
c-- new version
c      f(1,minus,minus)=+A7treeb_anom_wz(j2,j1,j4,j3,j6,j5,j7,zb,za)
      call A7treeb_anom_wz(j2,j1,j4,j3,j6,j5,j7,zb,za,
     & A7b_1,A7b_2,A7b_3,A7b_4)
      f(1,minus,minus)=
     &           A7b_1*(2._dp+xdelg1_z+xdelk_z+xlambda_z*s127/wmass**2)
     &          +A7b_2*(2._dp+xdelg1_z+xdelk_z+xlambda_z)
     &          +A7b_3*2._dp*(1._dp+xdelg1_z) 
     &          +A7b_4*xlambda_z/wmass**2
      f(10,minus,minus)=
     &           A7b_1*(2._dp+xdelg1_g+xdelk_g+xlambda_g*s127/wmass**2)
     &          +A7b_2*(2._dp+xdelg1_g+xdelk_g+xlambda_g)
     &          +A7b_3*2._dp*(1._dp+xdelg1_g) 
     &          +A7b_4*xlambda_g/wmass**2
      f(2,minus,minus)=-A7treea(j2,j1,j5,j6,j3,j4,j7,zb,za)
      f(3,minus,minus)=-A7treea(j2,j1,j4,j3,j6,j5,j7,zb,za)

c-- old version
c      f(1,minus,mplus)=-A7treeb_anom_wz(j2,j1,j6,j5,j3,j4,j7,zb,za)
c-- new version
c      f(1,minus,mplus)=+A7treeb_anom_wz(j2,j1,j4,j3,j5,j6,j7,zb,za)
      call A7treeb_anom_wz(j2,j1,j4,j3,j5,j6,j7,zb,za,
     & A7b_1,A7b_2,A7b_3,A7b_4)
      f(1,minus,mplus)=
     &           A7b_1*(2._dp+xdelg1_z+xdelk_z+xlambda_z*s127/wmass**2)
     &          +A7b_2*(2._dp+xdelg1_z+xdelk_z+xlambda_z)
     &          +A7b_3*2._dp*(1._dp+xdelg1_z) 
     &          +A7b_4*xlambda_z/wmass**2
      f(10,minus,mplus)=
     &           A7b_1*(2._dp+xdelg1_g+xdelk_g+xlambda_g*s127/wmass**2)
     &          +A7b_2*(2._dp+xdelg1_g+xdelk_g+xlambda_g)
     &          +A7b_3*2._dp*(1._dp+xdelg1_g) 
     &          +A7b_4*xlambda_g/wmass**2
      f(2,minus,mplus)=-A7treea(j2,j1,j6,j5,j3,j4,j7,zb,za)
      f(3,minus,mplus)=-A7treea(j2,j1,j4,j3,j5,j6,j7,zb,za)
            
      if (srdiags .eqv. .false.) return
      
c+++ the extra amplitude 7 accounts for W+ production instead of W-
c+++ amplitudes 8 and 9 are needed when considering ub+d instead of
c+++ .e+_dpub and are obtained from 6 and 7 by simply swapping 1 <+> 2

c+++ further review (12/05/01) shows that 8 and 9 are not needed,
c+++ hence they are now commented out

      f(4,mplus,minus)=+B7treeb(j2,j1,j6,j5,j4,j3,j7,za,zb)    
      f(5,mplus,minus)=+B7treea(j1,j2,j3,j4,j5,j6,j7,za,zb)    
      f(6,mplus,minus)=+B7treeb(j2,j1,j3,j4,j5,j6,j7,za,zb)    
      f(7,mplus,minus)=+B7treea(j1,j2,j6,j5,j4,j3,j7,za,zb)    
c      f(8,mplus,minus)=+B7treeb(j1,j2,j3,j4,j5,j6,j7,za,zb)    
c      f(9,mplus,minus)=+B7treea(j2,j1,j6,j5,j4,j3,j7,za,zb)    

      f(4,mplus,mplus)=+B7treeb(j2,j1,j5,j6,j4,j3,j7,za,zb)
      f(5,mplus,mplus)=+B7treea(j1,j2,j3,j4,j6,j5,j7,za,zb)
         
      f(4,minus,minus)=-B7treea(j2,j1,j4,j3,j6,j5,j7,zb,za)
      f(5,minus,minus)=-B7treeb(j1,j2,j5,j6,j3,j4,j7,zb,za)
      f(6,minus,minus)=-B7treea(j2,j1,j5,j6,j3,j4,j7,zb,za)
      f(7,minus,minus)=-B7treeb(j1,j2,j4,j3,j6,j5,j7,zb,za)
c      f(8,minus,minus)=-B7treea(j1,j2,j5,j6,j3,j4,j7,zb,za)
c      f(9,minus,minus)=-B7treeb(j2,j1,j4,j3,j6,j5,j7,zb,za)

      f(4,minus,mplus)=-B7treea(j2,j1,j3,j4,j6,j5,j7,zb,za)
      f(5,minus,mplus)=-B7treeb(j1,j2,j5,j6,j4,j3,j7,zb,za)
      f(4,minus,mplus)=-B7treea(j2,j1,j4,j3,j5,j6,j7,zb,za)
      f(5,minus,mplus)=-B7treeb(j1,j2,j6,j5,j3,j4,j7,zb,za)
     
      return
      end
