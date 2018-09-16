      subroutine qqb_totttH_gs(p,msqc)
      implicit none
      include 'types.f'
************************************************************************
*     Authors: J.M. Campbell and R.K. Ellis                            *
*     August, 2008.                                                    *
*     calculate the subtraction terms                                  *
*     averaged over initial colors and spins                           *
*     for the process                                                  *
*                                                                      *
*     q(-p1) +qbar(-p2)=t(p3)+t~(p4)+H(p5)                             * 
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
************************************************************************
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
     
      integer:: j,k
c --- remember: nd will count the dipoles
      integer:: nd
c--- slightly obtuse notation, to simplify declaration lines      
      real(dp):: p(mxpart,4),msqc(maxd,fn:nf,fn:nf)
      real(dp)::
     & msq16_2(fn:nf,fn:nf),msq26_1(fn:nf,fn:nf),
     & msq16_3(fn:nf,fn:nf),msq16_4(fn:nf,fn:nf),
     & msq36_1(fn:nf,fn:nf),msq46_1(fn:nf,fn:nf),
     & msq26_3(fn:nf,fn:nf),msq26_4(fn:nf,fn:nf),
     & msq36_2(fn:nf,fn:nf),msq46_2(fn:nf,fn:nf),
     & msq36_4(fn:nf,fn:nf),msq46_3(fn:nf,fn:nf),

     & msq16_2v(fn:nf,fn:nf),msq26_1v(fn:nf,fn:nf),
     & msq16_3v(fn:nf,fn:nf),msq16_4v(fn:nf,fn:nf),
     & msq36_1v(fn:nf,fn:nf),msq46_1v(fn:nf,fn:nf),
     & msq26_3v(fn:nf,fn:nf),msq26_4v(fn:nf,fn:nf),
     & msq36_2v(fn:nf,fn:nf),msq46_2v(fn:nf,fn:nf),
     & msq36_4v(fn:nf,fn:nf),msq46_3v(fn:nf,fn:nf),

     & sub16_2(4),sub26_1(4),
     & sub16_3(4),sub16_4(4),
     & sub36_1(4),sub46_1(4),
     & sub26_3(4),sub26_4(4),
     & sub36_2(4),sub46_2(4),
     & sub36_4(4),sub46_3(4),

     & sub16_2v,sub26_1v,
     & sub16_3v,sub16_4v,
     & sub36_1v,sub46_1v,
     & sub26_3v,sub26_4v,
     & sub36_2v,sub46_2v,
     & sub36_4v,sub46_3v

      real(dp)::
     & m16_2(0:2,fn:nf,fn:nf),m26_1(0:2,fn:nf,fn:nf),
     & m16_3(0:2,fn:nf,fn:nf),m36_1(0:2,fn:nf,fn:nf),
     & m16_4(0:2,fn:nf,fn:nf),m46_1(0:2,fn:nf,fn:nf),
     & m26_3(0:2,fn:nf,fn:nf),m36_2(0:2,fn:nf,fn:nf),
     & m26_4(0:2,fn:nf,fn:nf),m46_2(0:2,fn:nf,fn:nf),
     & m36_4(0:2,fn:nf,fn:nf),m46_3(0:2,fn:nf,fn:nf)

      real(dp)::
     & m16_2v(0:2,fn:nf,fn:nf),m26_1v(0:2,fn:nf,fn:nf),
     & m16_3v(0:2,fn:nf,fn:nf),m36_1v(0:2,fn:nf,fn:nf),
     & m16_4v(0:2,fn:nf,fn:nf),m46_1v(0:2,fn:nf,fn:nf),
     & m26_3v(0:2,fn:nf,fn:nf),m36_2v(0:2,fn:nf,fn:nf),
     & m26_4v(0:2,fn:nf,fn:nf),m46_2v(0:2,fn:nf,fn:nf),
     & m36_4v(0:2,fn:nf,fn:nf),m46_3v(0:2,fn:nf,fn:nf)


      external qqb_totttH,qqb_totttH_gvec,donothing_gvec

      
      qqproc=.true.
      qgproc=.false.
      gqproc=.false.
      ggproc=.false.

      ndmax=12

c-- initialize the matrix elements to zero
      do j=-nf,nf
      do k=-nf,nf      
      do nd=1,ndmax
        msqc(nd,j,k)=0d0
      enddo
      enddo
      enddo
      

c--initial-initial
      call dips_mass(1,p,1,6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     & qqb_totttH,qqb_totttH_gvec)
      call storedip_mass(m16_2,m16_2v)
      call dips_mass(2,p,2,6,1,sub26_1,sub26_1v,msq26_1,msq26_1v,
     & qqb_totttH,qqb_totttH_gvec)
      call storedip_mass(m26_1,m26_1v)

c--initial final
      call dips_mass(3,p,1,6,3,sub16_3,sub16_3v,msq16_3,msq16_3v,
     & qqb_totttH,qqb_totttH_gvec)
      call storedip_mass(m16_3,m16_3v)
      call dips_mass(4,p,2,6,3,sub26_3,sub26_3v,msq26_3,msq26_3v,
     & qqb_totttH,qqb_totttH_gvec)
      call storedip_mass(m26_3,m26_3v)
      call dips_mass(5,p,1,6,4,sub16_4,sub16_4v,msq16_4,msq16_4v,
     & qqb_totttH,qqb_totttH_gvec)
      call storedip_mass(m16_4,m16_4v)
      call dips_mass(6,p,2,6,4,sub26_4,sub26_4v,msq26_4,msq26_4v,
     & qqb_totttH,qqb_totttH_gvec)
      call storedip_mass(m26_4,m26_4v)

c--final-initial 
      call dips_mass(7,p,3,6,1,sub36_1,sub36_1v,msq36_1,msq36_1v,
     & qqb_totttH,donothing_gvec)
      call storedip_mass(m36_1,m36_1v)
      call dips_mass(8,p,3,6,2,sub36_2,sub36_2v,msq36_2,msq36_2v,
     & qqb_totttH,donothing_gvec)
      call storedip_mass(m36_2,m36_2v)
      call dips_mass(9,p,4,6,1,sub46_1,sub46_1v,msq46_1,msq46_1v,
     & qqb_totttH,donothing_gvec)
      call storedip_mass(m46_1,m46_1v)
      call dips_mass(10,p,4,6,2,sub46_2,sub46_2v,msq46_2,msq46_2v,
     & qqb_totttH,donothing_gvec)
      call storedip_mass(m46_2,m46_2v)

c--final-final
      call dips_mass(11,p,3,6,4,sub36_4,sub36_4v,msq36_4,msq36_4v,
     & qqb_totttH,donothing_gvec)
      call storedip_mass(m36_4,m36_4v)
      call dips_mass(12,p,4,6,3,sub46_3,sub46_3v,msq46_3,msq46_3v,
     & qqb_totttH,donothing_gvec)
      call storedip_mass(m46_3,m46_3v)


c      write(6,*) 'gs:sub16_2(gg)',sub16_2(gg)
c      write(6,*) 'gs:sub26_1(gg)',sub26_1(gg)

c      write(6,*) 'gs:sub16_3(gg)',sub16_3(gg)
c      write(6,*) 'gs:sub16_4(gg)',sub16_4(gg)
c      write(6,*) 'gs:sub26_3(gg)',sub26_3(gg)
c      write(6,*) 'gs:sub26_4(gg)',sub26_4(gg)

c      write(6,*) 'gs:sub36_1(qq)',sub36_1(qq)
c      write(6,*) 'gs:sub36_2(qq)',sub36_2(qq)
c      write(6,*) 'gs:sub46_1(qq)',sub46_1(qq)
c      write(6,*) 'gs:sub46_2(qq)',sub46_2(qq)

c      write(6,*) 'gs:sub36_4(qq)',sub36_4(qq)
c      write(6,*) 'gs:sub46_3(qq)',sub46_3(qq)


c--- fill the dipole contributions
      do j=-nf,nf
      do k=-nf,nf

      if ((j == 0) .and. (k==0)) then
      msqc(1,j,k)= (+m16_2(1,j,k)+m16_2(2,j,k))*sub16_2(gg)*xn
     &            +(+m16_2v(1,j,k)+m16_2v(2,j,k))*sub16_2v*xn
      msqc(2,j,k)= (+m26_1(1,j,k)+m26_1(2,j,k))*sub26_1(gg)*xn
     &            +(+m26_1v(1,j,k)+m26_1v(2,j,k))*sub26_1v*xn
      
      msqc(3,j,k)= (+m16_3(2,j,k)+m16_3(0,j,k))*sub16_3(gg)*xn
     &            +(+m16_3v(2,j,k)+m16_3v(0,j,k))*sub16_3v*xn
      msqc(4,j,k)= (+m26_3(1,j,k)+m26_3(0,j,k))*sub26_3(gg)*xn
     &            +(+m26_3v(1,j,k)+m26_3v(0,j,k))*sub26_3v*xn

      msqc(5,j,k)= (+m16_4(1,j,k)+m16_4(0,j,k))*sub16_4(gg)*xn
     &            +(+m16_4v(1,j,k)+m16_4v(0,j,k))*sub16_4v*xn
      msqc(6,j,k)= (+m26_4(2,j,k)+m26_4(0,j,k))*sub26_4(gg)*xn
     &            +(+m26_4v(2,j,k)+m26_4v(0,j,k))*sub26_4v*xn

      msqc(7,j,k)= (+m36_1(2,j,k)+m36_1(0,j,k))*sub36_1(qq)*xn
      msqc(8,j,k)= (+m36_2(1,j,k)+m36_2(0,j,k))*sub36_2(qq)*xn

      msqc(9,j,k)= (+m46_1(1,j,k)+m46_1(0,j,k))*sub46_1(qq)*xn
      msqc(10,j,k)=(+m46_2(2,j,k)+m46_2(0,j,k))*sub46_2(qq)*xn

      msqc(11,j,k)=-((m36_4(1,j,k)+m36_4(2,j,k))/xn
     &              +(xn+1d0/xn)*m36_4(0,j,k))*sub36_4(qq)
      msqc(12,j,k)=-((m46_3(1,j,k)+m46_3(2,j,k))/xn
     &              +(xn+1d0/xn)*m46_3(0,j,k))*sub46_3(qq)

      elseif ((j > 0) .and. (k == -j)) then
      msqc(1,j,k)= -msq16_2(j,k)*sub16_2(qq)/xn
      msqc(2,j,k)= -msq26_1(j,k)*sub26_1(qq)/xn

      msqc(3,j,k)= +msq16_3(j,k)*sub16_3(qq)*(xn-2d0/xn)
      msqc(4,j,k)= +msq26_3(j,k)*sub26_3(qq)*2d0/xn
      msqc(5,j,k)= +msq16_4(j,k)*sub16_4(qq)*2d0/xn
      msqc(6,j,k)= +msq26_4(j,k)*sub26_4(qq)*(xn-2d0/xn)

      msqc(7,j,k)= +msq36_1(j,k)*sub36_1(qq)*(xn-2d0/xn)
      msqc(8,j,k)= +msq36_2(j,k)*sub36_2(qq)*2d0/xn
      msqc(9,j,k)= +msq46_1(j,k)*sub46_1(qq)*2d0/xn
      msqc(10,j,k)=+msq46_2(j,k)*sub46_2(qq)*(xn-2d0/xn)

      msqc(11,j,k)=-msq36_4(j,k)*sub36_4(qq)/xn
      msqc(12,j,k)=-msq46_3(j,k)*sub46_3(qq)/xn

      elseif ((j < 0) .and. (k == -j)) then
      msqc(1,j,k)= -msq16_2(j,k)*sub16_2(qq)/xn
      msqc(2,j,k)= -msq26_1(j,k)*sub26_1(qq)/xn

      msqc(3,j,k)= +msq16_3(j,k)*sub16_3(qq)*two/xn
      msqc(4,j,k)= +msq26_3(j,k)*sub26_3(qq)*(xn-two/xn)
      msqc(5,j,k)= +msq16_4(j,k)*sub16_4(qq)*(xn-two/xn)
      msqc(6,j,k)= +msq26_4(j,k)*sub26_4(qq)*two/xn

      msqc(7,j,k)= +msq36_1(j,k)*sub36_1(qq)*two/xn
      msqc(8,j,k)= +msq36_2(j,k)*sub36_2(qq)*(xn-two/xn)
      msqc(9,j,k)= +msq46_1(j,k)*sub46_1(qq)*(xn-two/xn)
      msqc(10,j,k)=+msq46_2(j,k)*sub46_2(qq)*two/xn

      msqc(11,j,k)=-msq36_4(j,k)*sub36_4(qq)/xn
      msqc(12,j,k)=-msq46_3(j,k)*sub46_3(qq)/xn

      elseif ((j == 0) .and. (k .ne. 0)) then
      msqc(1,j,k)= 2d0*tr*msq16_2(-k,k)*sub16_2(qg)
      msqc(2,j,k)= (aveqg/avegg)*
     &    (sub26_1(gq)*msq26_1(0,0)+sub26_1v*msq26_1v(0,0))

      elseif ((j .ne. 0) .and. (k == 0)) then
      msqc(1,j,k)= (aveqg/avegg)*
     &    (sub16_2(gq)*msq16_2(0,0)+sub16_2v*msq16_2v(0,0))
      msqc(2,j,k)= 2d0*tr*msq26_1(j,-j)*sub26_1(qg)
      endif

      enddo
      enddo

      return
      end            
      
