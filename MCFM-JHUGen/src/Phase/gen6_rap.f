      subroutine gen6_rap(r,p,wt6,*)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'mxdim.f'
      include 'breit.f'

      integer:: nu
      real(dp):: r(mxdim)
      real(dp):: wt6,p(mxpart,4),tp(4),tm(4),bp(4),bm(4)
      real(dp):: wtepnn,wtnbem,ep(4),em(4),nn(4),nb(4),wp(4),wm(4)
      real(dp):: wtttb,wtwp,wtwm,s3min,wt0
      parameter(wt0=1._dp/twopi**4)


*     q(-p1) +qbar(-p2)=t(nu(p3)+e^+(p4)+b(p5))                        *
*                       +t~(b~(p6)+e^-(p7)+nu(p8))+g(p9)               *
*                                                                      * 


      wt6=0._dp
      mass2=mt
      mass3=mt
      call gen2m(r,p,wtttb,*999)
      wtttb=(mt*twidth*pi)**2*wtttb
      do nu=1,4
      tp(nu)=p(3,nu)
      tm(nu)=p(4,nu)
      enddo
      s3min=0._dp
      n3=1
      mass3=wmass      
      width3=wwidth      

      call phi1_2m(mb,r(5),r(6),r(7),s3min,tp,bm,wp,wtwp,*999)
      call phi1_2m(mb,r(8),r(9),r(10),s3min,tm,bp,wm,wtwm,*999)
      call phi3m0(r(11),r(12),wp,nn,ep,wtepnn,*999)
      call phi3m0(r(13),r(14),wm,em,nb,wtnbem,*999)
      wt6=wt0*wtepnn*wtnbem*wtwp*wtwm*wtttb
*     q(-p1) +qbar(-p2)=t(nu(p3)+e^+(p4)+b(p5))                        *
*                       +t~(b~(p6)+e^-(p7)+nu(p8))               *
*                                                                      * 
      do nu=1,4
      p(3,nu)=nn(nu)
      p(4,nu)=ep(nu)
      p(5,nu)=bm(nu)

      p(6,nu)=bp(nu)
      p(7,nu)=em(nu)
      p(8,nu)=nb(nu)

      enddo

      return
 999  wt6=0._dp
      p(:,:)=0._dp
      return 1
      
      end

