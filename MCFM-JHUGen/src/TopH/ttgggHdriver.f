      subroutine ttgggHdriver(q,ampsq)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: q(mxpart,4),ampsq
      integer:: h1,h2,h3,h4,h6
      complex(dp):: m(6),mqed,
     & m126(2,2,2,2,2),m261(2,2,2,2,2),m612(2,2,2,2,2),
     & m621(2,2,2,2,2),m216(2,2,2,2,2),m162(2,2,2,2,2)

      call ttgggHamp(q,3,4,1,2,6,1,1,m126)
      call ttgggHamp(q,3,4,2,6,1,1,1,m261)
      call ttgggHamp(q,3,4,6,1,2,1,1,m612)
      call ttgggHamp(q,3,4,6,2,1,1,1,m621)
      call ttgggHamp(q,3,4,2,1,6,1,1,m216)
      call ttgggHamp(q,3,4,1,6,2,1,1,m162)
      ampsq=0d0
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do h6=1,2
      m(1)=m126(h1,h2,h3,h4,h6)
      m(2)=m261(h1,h2,h4,h6,h3)
      m(3)=m612(h1,h2,h6,h3,h4)
      m(4)=m621(h1,h2,h6,h4,h3)
      m(5)=m216(h1,h2,h4,h3,h6)
      m(6)=m162(h1,h2,h3,h6,h4)
      mqed=m(1)+m(2)+m(3)+m(4)+m(5)+m(6)
C-----Overall factor of V removed from color sum
      ampsq=ampsq+(xnsq-2d0)*(
     & +real(m(1)*conjg(m(1)))+real(m(2)*conjg(m(2)))
     & +real(m(3)*conjg(m(3)))+real(m(4)*conjg(m(4)))
     & +real(m(5)*conjg(m(5)))+real(m(6)*conjg(m(6))))
       ampsq=ampsq
     & +real(m(1)*conjg(m(4)-m(6)-m(5)))
     & +real(m(2)*conjg(m(6)-m(5)-m(4)))
     & +real(m(3)*conjg(m(5)-m(4)-m(6)))
     & +real(m(4)*conjg(m(1)-m(2)-m(3)))
     & +real(m(5)*conjg(m(3)-m(1)-m(2)))
     & +real(m(6)*conjg(m(2)-m(3)-m(1)))
      ampsq=ampsq+real(mqed*conjg(mqed))/xn**2
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
