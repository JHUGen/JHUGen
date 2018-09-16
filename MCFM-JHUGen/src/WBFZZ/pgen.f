      subroutine pgen(p,jcount,i3,i4,pout)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      integer:: jcount,i3,i4
      real(dp):: p(mxpart,4),pout(mxpart,4),
     & p3in(4),p3out(4),p4out(4),p34(4)

      pout(:,:)=p(:,:)
      p34(:)=p(i3,:)+p(i4,:)

      p3in(4)=+0.5d0*zmass

      if ((jcount == 1) .or. (jcount == 2)) then
             p3in(1:3)=0d0
             p3in(1)=+0.5d0*zmass

             call boost(zmass,p34,p3in,p3out)
             p4out(:)=p34(:)-p3out(:)

             if (jcount == 1) then
             pout(i3,:)=p3out(:)
             pout(i4,:)=p4out(:)
             elseif (jcount == 2) then
             pout(i3,:)=p4out(:)
             pout(i4,:)=p3out(:)
             endif

      elseif ((jcount == 3) .or. (jcount == 4)) then
             p3in(1:3)=0d0
             p3in(2)=+0.5d0*zmass

             call boost(zmass,p34,p3in,p3out)
             p4out(:)=p34(:)-p3out(:)


             if (jcount == 3) then
             pout(i3,:)=p3out(:)
             pout(i4,:)=p4out(:)
             elseif (jcount == 4) then
             pout(i3,:)=p4out(:)
             pout(i4,:)=p3out(:)
             endif

      elseif ((jcount == 5) .or. (jcount == 6)) then
             p3in(1:3)=0d0
             p3in(3)=+0.5d0*zmass

             call boost(zmass,p34,p3in,p3out)
             p4out(:)=p34(:)-p3out(:)

             if (jcount == 5) then
             pout(i3,:)=p3out(:)
             pout(i4,:)=p4out(:)
             elseif (jcount == 6) then
             pout(i3,:)=p4out(:)
             pout(i4,:)=p3out(:)
             endif

      else
      write(6,*) 'Unimplemnted value of j'
      stop
      endif

      return
      end
