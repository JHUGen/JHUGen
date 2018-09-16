      subroutine pgen2(p,jcount,i3,i4,pout)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      integer:: jcount,i3,i4
      real(dp):: p(mxpart,4),pout(mxpart,4),p3out(4),p4out(4)

C----now setup pout
      pout(:,:)=p(:,:)


      p3out(4)=0.5d0*zmass
      p4out(4)=0.5d0*zmass

      if ((jcount == 1) .or. (jcount == 2)) then
             p3out(1:3)=0d0
             p4out(1:3)=0d0
             p3out(1)=+0.5d0*zmass
             p4out(1)=-0.5d0*zmass

             pout(i3,:)=p3out(:)
             pout(i4,:)=p4out(:)
             if (jcount == 2) p(i3,:)=p4out(:)
             if (jcount == 2) p(i4,:)=p3out(:)

      elseif ((jcount == 3) .or. (jcount == 4)) then
             p3out(1:3)=0d0
             p4out(1:3)=0d0
             p3out(2)=+0.5d0*zmass
             p4out(2)=-0.5d0*zmass

             pout(i3,:)=p3out(:)
             pout(i4,:)=p4out(:)
             if (jcount == 4) pout(i3,:)=p4out(:)
             if (jcount == 4) pout(i4,:)=p3out(:)

      elseif ((jcount == 5) .or. (jcount == 6)) then
             p3out(1:3)=0d0
             p4out(1:3)=0d0
             p3out(3)=+0.5d0*zmass
             p4out(3)=-0.5d0*zmass

             pout(i3,:)=p3out(:)
             pout(i4,:)=p4out(:)
             if (jcount == 6) pout(i3,:)=p4out(:)
             if (jcount == 6) pout(i4,:)=p3out(:)
      else
      write(6,*) 'Unimplemnted value of j'
      stop
      endif
c      write(6,*)'j',jcount
c      write(6,*)'p3',pout(3,4)**2-pout(3,1)**2-pout(3,2)**2-pout(3,3)**2
c      write(6,*)'p4',pout(4,4)**2-pout(4,1)**2-pout(4,2)**2-pout(4,3)**2
c      write(6,*)'p5',pout(5,4)**2-pout(5,1)**2-pout(5,2)**2-pout(5,3)**2
c      write(6,*)'p6',pout(6,4)**2-pout(6,1)**2-pout(6,2)**2-pout(6,3)**2
c      write(6,*)'p7',pout(7,4)**2-pout(7,1)**2-pout(7,2)**2-pout(7,3)**2
c      write(6,*)'p8',pout(8,4)**2-pout(8,1)**2-pout(8,2)**2-pout(8,3)**2
c      pause
      return
      end
