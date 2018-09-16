      subroutine gen7_rap(r,p,wt7,*)
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
      real(dp):: wt7,p(mxpart,4),tp(4),tm(4),bp(4),bm(4),
     & nn(4),nb(4)
      real(dp):: wtepnn,wtnbem,ep(4),em(4),pg(4)
      real(dp):: wp(4),wm(4),m3,m4
      real(dp):: wt0,wtttg,wtwp,wtwm,s3min
      parameter(wt0=1._dp/twopi**4)

*     q(-p1) +qbar(-p2)=t(nu(p3)+e^+(p4)+b(p5))                        *
*                       +t~(b~(p6)+e^-(p7)+nu(p8))+g(p9)               *
*                                                                      * 
      wt7=0._dp
      m3=mt
      m4=mt
      call gen3m_rap(r,p,m3,m4,wtttg,*999)
      wtttg=(pi*mt*twidth)**2*wtttg
      do nu=1,4
      tp(nu)=p(3,nu)
      tm(nu)=p(4,nu)
      pg(nu)=p(5,nu)
      enddo
      s3min=0._dp
      n3=1
      mass3=wmass      
      width3=wwidth      

      call phi1_2m(mb,r(8),r(9),r(10),s3min,tp,bm,wp,wtwp,*999)
      call phi1_2m(mb,r(11),r(12),r(13),s3min,tm,bp,wm,wtwm,*999)
      call phi3m0(r(14),r(15),wp,nn,ep,wtepnn,*999)
      call phi3m0(r(16),r(17),wm,em,nb,wtnbem,*999)

      wt7=wt0*wtepnn*wtnbem*wtwp*wtwm*wtttg
*     q(-p1) +qbar(-p2)=t(nu(p3)+e^+(p4)+b(p5))                        *
*                       +t~(b~(p6)+e^-(p7)+nu(p8))+g(p9)               *
*                                                                      * 
      do nu=1,4
      p(3,nu)=nn(nu)
      p(4,nu)=ep(nu)
      p(5,nu)=bm(nu)

      p(6,nu)=bp(nu)
      p(7,nu)=em(nu)
      p(8,nu)=nb(nu)

      p(9,nu)=pg(nu)
      enddo
      return
 999  wt7=0._dp
      p(:,:)=0._dp
      return 1
      
      end

