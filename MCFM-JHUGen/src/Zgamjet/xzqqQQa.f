      subroutine xzqqQQa_qq(i1,i2,i3,i4,i5,i6,i7,ai,bi)
      implicit none
      include 'types.f'
      
************************************************************************
*     0 ---> q(p1)+qb(p2)+Q(p3)+Qb(p4)+ph(p5)+lb(p6)+l(p7)              *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcharge.f'
      integer:: i1,i2,i3,i4,i5,i6,i7
      complex(dp)::a1(2,2,2,2),a2(2,2,2,2),a3(2,2,2,2),a4(2,2,2,2),
     &               b1(2,2,2,2),b2(2,2,2,2),b3(2,2,2,2),b4(2,2,2,2)
      complex(dp):: ai(2,2,2,2,2,2),bi(2,2,2,2,2,2)
      integer:: hq,Qh,ha,lh,q12,q34

c-----'a' corresponds to the perm 1,2,3,4,5
      call nagyqqQQg(i1,i2,i3,i4,i5,i6,i7,a1,a2,a3,a4)
c-----'b' corresponds to the perm 3,4,1,2,5
      call nagyqqQQg(i3,i4,i1,i2,i5,i6,i7,b1,b2,b3,b4)
c-----sum ai and bi
      do q12=1,2
      do q34=1,2
      do hq=1,2
      do Qh=1,2
      do ha=1,2
      do lh=1,2
      ai(q12,q34,hq,Qh,ha,lh)= Q(q12)*(a1(hq,Qh,ha,lh)+a2(hq,Qh,ha,lh))
     &                        +Q(q34)*(a3(hq,Qh,ha,lh)+a4(hq,Qh,ha,lh))
      bi(q12,q34,hq,Qh,ha,lh)= Q(q34)*(b1(hq,Qh,ha,lh)+b2(hq,Qh,ha,lh))
     &                        +Q(q12)*(b3(hq,Qh,ha,lh)+b4(hq,Qh,ha,lh))
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c-----
      return
      end


      subroutine xzqqQQa_ql(i1,i2,i3,i4,i5,i6,i7,ai,bi)
      implicit none
      include 'types.f'
      
************************************************************************
*     0 ---> q(p1)+qb(p2)+Q(p3)+Qb(p4)+ph(p5)+lb(p6)+l(p7)              *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: i1,i2,i3,i4,i5,i6,i7
      integer:: w1,w2,w3,w4
      complex(dp):: ai(2,2,2,2),bi(2,2,2,2)
      integer:: hq,Qh,ha,lh
c-----'a' corresponds to the perm 1,2,3,4,5
      call amp_zqqQQa_ql(i1,i2,i3,i4,i5,i6,i7,ai)
c-----'b' corresponds to the perm 3,4,1,2,5
      call amp_zqqQQa_ql(i3,i4,i1,i2,i5,i6,i7,bi)
c-----
      return
      end

