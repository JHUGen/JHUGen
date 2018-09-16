      subroutine wwamps(j1,j2,j3,j4,j5,j6,j7,za,zb,f)
      implicit none
      include 'types.f'
c  -first label of fs,ft is gluon polarization, second is qqb line
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'zerowidth.f'
      include 'masses.f'
      include 'anomcoup.f'
      include 'xanomcoup.f'
      include 'kprocess.f'
      include 'srdiags.f'      
      integer:: j,k,jtype,j1,j2,j3,j4,j5,j6,j7,mplus,minus
      complex(dp):: A7treea,B7treea,B7treeb
      complex(dp):: f(5,2,2),A7b_1,A7b_2,A7b_3
      complex(dp):: prop34,prop56,propboth
      parameter (minus=1,mplus=2)
      
c----initialize to zero
      do jtype=3,4
      do j=1,2
      do k=1,2
      f(jtype,j,k)=czip
      enddo
      enddo
      enddo
            
      if     (zerowidth  .eqv. .true.) then
      prop34=s(3,4)/cplx2(s(3,4)-wmass**2,wmass*wwidth)
      prop56=s(5,6)/cplx2(s(5,6)-wmass**2,wmass*wwidth)
      elseif (zerowidth .neqv. .true.) then
      prop34=cplx1(s(3,4)/(s(3,4)-wmass**2))
      prop56=cplx1(s(5,6)/(s(5,6)-wmass**2)) 
      endif
      propboth=prop34*prop56

      f(1,mplus,mplus)= czip
c      f(2,mplus,mplus)=-A7treeb(j2,j1,j3,j4,j5,j6,j7,za,zb)*propboth
      call A7treeb_anom(j2,j1,j3,j4,j5,j6,j7,za,zb,A7b_1,A7b_2,A7b_3)
      f(2,mplus,mplus)=-(A7b_1*(2._dp+xdelg1_z+xdelk_z+xlambda_z)
     &                 +A7b_2*(2._dp*(1._dp+xdelg1_z))
     &                 +A7b_3*(xlambda_z/wmass**2))*propboth
      f(5,mplus,mplus)=-(A7b_1*(2._dp+xdelg1_g+xdelk_g+xlambda_g)
     &                 +A7b_2*(2._dp*(1._dp+xdelg1_g))
     &                 +A7b_3*(xlambda_g/wmass**2))*propboth

      f(1,mplus,minus)=+A7treea(j1,j2,j3,j4,j5,j6,j7,za,zb)*propboth      
c      f(2,mplus,minus)=+A7treeb(j1,j2,j3,j4,j5,j6,j7,za,zb)*propboth
      call A7treeb_anom(j1,j2,j3,j4,j5,j6,j7,za,zb,A7b_1,A7b_2,A7b_3)
      f(2,mplus,minus)=(A7b_1*(2._dp+xdelg1_z+xdelk_z+xlambda_z)
     &                 +A7b_2*(2._dp*(1._dp+xdelg1_z))
     &                 +A7b_3*(xlambda_z/wmass**2))*propboth
      f(5,mplus,minus)=(A7b_1*(2._dp+xdelg1_g+xdelk_g+xlambda_g)
     &                 +A7b_2*(2._dp*(1._dp+xdelg1_g))
     &                 +A7b_3*(xlambda_g/wmass**2))*propboth

      f(1,minus,mplus)= czip
c      f(2,minus,mplus)=+A7treeb(j1,j2,j5,j6,j3,j4,j7,zb,za)*propboth
      call A7treeb_anom(j1,j2,j5,j6,j3,j4,j7,zb,za,A7b_1,A7b_2,A7b_3)
      f(2,minus,mplus)=(A7b_1*(2._dp+xdelg1_z+xdelk_z+xlambda_z)
     &                 +A7b_2*(2._dp*(1._dp+xdelg1_z))
     &                 +A7b_3*(xlambda_z/wmass**2))*propboth
      f(5,minus,mplus)=(A7b_1*(2._dp+xdelg1_g+xdelk_g+xlambda_g)
     &                 +A7b_2*(2._dp*(1._dp+xdelg1_g))
     &                 +A7b_3*(xlambda_g/wmass**2))*propboth

      f(1,minus,minus)=-A7treea(j2,j1,j5,j6,j3,j4,j7,zb,za)*propboth
c      f(2,minus,minus)=-A7treeb(j2,j1,j5,j6,j3,j4,j7,zb,za)*propboth
      call A7treeb_anom(j2,j1,j5,j6,j3,j4,j7,zb,za,A7b_1,A7b_2,A7b_3)
      f(2,minus,minus)=-(A7b_1*(2._dp+xdelg1_z+xdelk_z+xlambda_z)
     &                 +A7b_2*(2._dp*(1._dp+xdelg1_z))
     &                 +A7b_3*(xlambda_z/wmass**2))*propboth
      f(5,minus,minus)=-(A7b_1*(2._dp+xdelg1_g+xdelk_g+xlambda_g)
     &                 +A7b_2*(2._dp*(1._dp+xdelg1_g))
     &                 +A7b_3*(xlambda_g/wmass**2))*propboth

      if (srdiags .eqv. .false.) return  
c--- Done all non-singly resonant amplitudes

c--- also return here for WW+jet process since no singly-resonant
c---  diagrams are included in the real contribution
      if (kcase==kWW_jet) return    

      f(3,mplus,mplus)=-B7treeb(j1,j2,3,4,5,6,j7,za,zb)*prop34
     &                 -B7treea(j2,j1,3,4,5,6,j7,za,zb)*prop56
      f(4,mplus,mplus)=-B7treea(j2,j1,6,5,4,3,j7,za,zb)*prop34
     &                 -B7treeb(j1,j2,6,5,4,3,j7,za,zb)*prop56

      f(3,mplus,minus)=+B7treeb(j2,j1,3,4,5,6,j7,za,zb)*prop34
     &                 +B7treea(j1,j2,3,4,5,6,j7,za,zb)*prop56
      f(4,mplus,minus)=+B7treea(j1,j2,6,5,4,3,j7,za,zb)*prop34
     &                 +B7treeb(j2,j1,6,5,4,3,j7,za,zb)*prop56

      f(3,minus,mplus)=+B7treea(j1,j2,5,6,3,4,j7,zb,za)*prop34
     &                 +B7treeb(j2,j1,5,6,3,4,j7,zb,za)*prop56
      f(4,minus,mplus)=+B7treeb(j2,j1,4,3,6,5,j7,zb,za)*prop34
     &                 +B7treea(j1,j2,4,3,6,5,j7,zb,za)*prop56

      f(3,minus,minus)=-B7treea(j2,j1,5,6,3,4,j7,zb,za)*prop34
     &                 -B7treeb(j1,j2,5,6,3,4,j7,zb,za)*prop56
      f(4,minus,minus)=-B7treeb(j1,j2,4,3,6,5,j7,zb,za)*prop34
     &                 -B7treea(j2,j1,4,3,6,5,j7,zb,za)*prop56

      return
      end

