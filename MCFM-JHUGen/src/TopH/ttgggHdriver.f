      subroutine ttgggHdriver(q,ampsq)
      implicit none
      include 'constants.f'
      double precision q(mxpart,4),ampsq
      integer h1,h2,h3,h4,h6
      double complex m(6),mqed,
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
     & +dble(m(1)*Dconjg(m(1)))+dble(m(2)*Dconjg(m(2)))
     & +dble(m(3)*Dconjg(m(3)))+dble(m(4)*Dconjg(m(4)))
     & +dble(m(5)*Dconjg(m(5)))+dble(m(6)*Dconjg(m(6))))
       ampsq=ampsq
     & +dble(m(1)*dconjg(m(4)-m(6)-m(5)))
     & +dble(m(2)*dconjg(m(6)-m(5)-m(4)))
     & +dble(m(3)*dconjg(m(5)-m(4)-m(6)))
     & +dble(m(4)*dconjg(m(1)-m(2)-m(3)))
     & +dble(m(5)*dconjg(m(3)-m(1)-m(2)))
     & +dble(m(6)*dconjg(m(2)-m(3)-m(1)))
      ampsq=ampsq+dble(mqed*dconjg(mqed))/xn**2
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
