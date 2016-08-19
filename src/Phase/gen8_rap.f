      subroutine gen8_rap(r,p,wt8,*)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'breit.f'
      integer nu

      double precision r(mxdim)
      double precision wt8,p(mxpart,4),tp(4),tm(4),bp(4),bm(4),n(4),e(4)
      double precision wtepnn,wtnbem,ep(4),em(4),nn(4),nb(4),wp(4),wm(4)
      double precision wtttw,wtwp,wtwm,s3min,wt0
      parameter(wt0=1d0/twopi**4)

*     q(-p1) +qbar(-p2)=t(nu(p3)+e^+(p4)+b(p5))                        *
*                       +t~(b~(p6)+e^-(p7)+nu(p8))+g(p9)               *
*                                                                      * 

      wt8=0d0
C---call gen4 that uses r(1)....r(10)
      call gen4(r,p,wtttw,*999)
      wtttw=(mt*twidth*pi)**2*wtttw
      do nu=1,4
      tp(nu)=p(5,nu)
      tm(nu)=p(6,nu)
      n(nu)=p(3,nu)
      e(nu)=p(4,nu)
      
      enddo

      s3min=0d0
      n3=1
      mass3=wmass      
      width3=wwidth      
      
      call phi1_2m(mb,r(11),r(12),r(13),s3min,tp,bm,wp,wtwp,*999)
      call phi1_2m(mb,r(14),r(15),r(16),s3min,tm,bp,wm,wtwm,*999)
      call phi3m0(r(17),r(18),wp,nn,ep,wtepnn,*999)
      call phi3m0(r(19),r(20),wm,em,nb,wtnbem,*999)
      wt8=wt0*wtepnn*wtnbem*wtwp*wtwm*wtttw
*     q(-p1) +qbar(-p2)=t(nu(p3)+e^+(p4)+b(p5))                        *
*                       +t~(b~(p6)+e^-(p7)+nu(p8)+W(nu(p9)+e+(p10)     *
*                                                                      * 
      do nu=1,4
      p(3,nu)=nn(nu)
      p(4,nu)=ep(nu)
      p(5,nu)=bm(nu)

      p(6,nu)=bp(nu)
      p(7,nu)=em(nu)
      p(8,nu)=nb(nu)
      p(9,nu)=n(nu)
      p(10,nu)=e(nu)
      enddo

      return
 999  wt8=0d0
      p(:,:)=0d0
      return 1
      
      end

