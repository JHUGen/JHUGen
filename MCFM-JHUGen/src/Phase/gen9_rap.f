      subroutine gen9_rap(r,p,wt9,*)
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
      real(dp):: wt9,p(mxpart,4),tp(4),tm(4),bp(4),bm(4),
     & n(4),e(4),g(4),ep(4),em(4),nn(4),nb(4),wp(4),wm(4),
     & wtepnn,wtnbem
      real(dp):: wtttwg,wtwp,wtwm,s3min,wt0
      parameter(wt0=1._dp/twopi**4)

*     q(-p1) +qbar(-p2)=t(nu(p3)+e^+(p4)+b(p5))                        *
*                       +t~(b~(p6)+e^-(p7)+nu(p8))+g(p9)               *
*                                                                      * 


      wt9=0._dp
C---call gen5 that uses r(1)....r(13)
      call gen5(r,p,wtttwg,*999)
      wtttwg=(mt*twidth*pi)**2*wtttwg
      do nu=1,4
      n(nu)=p(3,nu)
      e(nu)=p(4,nu)
      tp(nu)=p(5,nu)
      tm(nu)=p(6,nu)
      g(nu)=p(7,nu)
      enddo
      

      s3min=0._dp
      n3=1
      mass3=wmass
      width3=wwidth

      call phi1_2m(mb,r(14),r(15),r(16),s3min,tp,bm,wp,wtwp,*999)
      call phi1_2m(mb,r(17),r(18),r(19),s3min,tm,bp,wm,wtwm,*999)
      call phi3m0(r(20),r(21),wp,nn,ep,wtepnn,*999)
      call phi3m0(r(22),r(23),wm,em,nb,wtnbem,*999)
      wt9=wt0*wtepnn*wtnbem*wtwp*wtwm*wtttwg
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
      p(11,nu)=g(nu)
      enddo

      return
 999  wt9=0._dp
      p(:,:)=0._dp
      return 1
      
      end

