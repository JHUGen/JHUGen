      subroutine qqb_QQbdku_gs(p,msqc)
      implicit none
      include 'types.f'
************************************************************************
*     Authors: J.M. Campbell and R.K. Ellis                            *
*     August, 2008.                                                    *
*     calculate the subtraction terms                                  *
*     averaged over initial colors and spins                           *
*     for the process                                                  *
*                                                                      *
*     q(-p1) +qbar(-p2)=bbar(p6)+e-(p7)+nubar(p8)+nu(p3)+e+(p4)+b(p5)  *
*                                                                      *
*     Top is kept strictly on-shell; spin correlations are dropped     *
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
     & msq15_2(fn:nf,fn:nf),msq25_1(fn:nf,fn:nf),
     & msq15_3(fn:nf,fn:nf),msq15_4(fn:nf,fn:nf),
     & msq35_1(fn:nf,fn:nf),msq45_1(fn:nf,fn:nf),
     & msq25_3(fn:nf,fn:nf),msq25_4(fn:nf,fn:nf),
     & msq35_2(fn:nf,fn:nf),msq45_2(fn:nf,fn:nf),
     & msq35_4(fn:nf,fn:nf),msq45_3(fn:nf,fn:nf),

     & msq15_2v(fn:nf,fn:nf),msq25_1v(fn:nf,fn:nf),
     & msq15_3v(fn:nf,fn:nf),msq15_4v(fn:nf,fn:nf),
     & msq35_1v(fn:nf,fn:nf),msq45_1v(fn:nf,fn:nf),
     & msq25_3v(fn:nf,fn:nf),msq25_4v(fn:nf,fn:nf),
     & msq35_2v(fn:nf,fn:nf),msq45_2v(fn:nf,fn:nf),
     & msq35_4v(fn:nf,fn:nf),msq45_3v(fn:nf,fn:nf),

     & sub15_2(4),sub25_1(4),
     & sub15_3(4),sub15_4(4),
     & sub35_1(4),sub45_1(4),
     & sub25_3(4),sub25_4(4),
     & sub35_2(4),sub45_2(4),
     & sub35_4(4),sub45_3(4),

     & sub15_2v,sub25_1v,
     & sub15_3v,sub15_4v,
     & sub35_1v,sub45_1v,
     & sub25_3v,sub25_4v,
     & sub35_2v,sub45_2v,
     & sub35_4v,sub45_3v

      real(dp)::
     & m15_2(0:2,fn:nf,fn:nf),m25_1(0:2,fn:nf,fn:nf),
     & m15_3(0:2,fn:nf,fn:nf),m35_1(0:2,fn:nf,fn:nf),
     & m15_4(0:2,fn:nf,fn:nf),m45_1(0:2,fn:nf,fn:nf),
     & m25_3(0:2,fn:nf,fn:nf),m35_2(0:2,fn:nf,fn:nf),
     & m25_4(0:2,fn:nf,fn:nf),m45_2(0:2,fn:nf,fn:nf),
     & m35_4(0:2,fn:nf,fn:nf),m45_3(0:2,fn:nf,fn:nf)

      real(dp)::
     & m15_2v(0:2,fn:nf,fn:nf),m25_1v(0:2,fn:nf,fn:nf),
     & m15_3v(0:2,fn:nf,fn:nf),m35_1v(0:2,fn:nf,fn:nf),
     & m15_4v(0:2,fn:nf,fn:nf),m45_1v(0:2,fn:nf,fn:nf),
     & m25_3v(0:2,fn:nf,fn:nf),m35_2v(0:2,fn:nf,fn:nf),
     & m25_4v(0:2,fn:nf,fn:nf),m45_2v(0:2,fn:nf,fn:nf),
     & m35_4v(0:2,fn:nf,fn:nf),m45_3v(0:2,fn:nf,fn:nf)


      external qqb_QQbdku,qqb_QQbdku_gvec,donothing_gvec


      qqproc=.true.
      qgproc=.false.
      gqproc=.false.
      ggproc=.false.

      ndmax=12

c-- initialize the matrix elements to zero
      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msqc(nd,j,k)=0._dp
      enddo
      enddo
      enddo


c--initial-initial
      call dips_mass(1,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & qqb_QQbdku,qqb_QQbdku_gvec)
      call storedip_mass(m15_2,m15_2v)
      call dips_mass(2,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     & qqb_QQbdku,qqb_QQbdku_gvec)
      call storedip_mass(m25_1,m25_1v)

c--initial final
      call dips_mass(3,p,1,5,3,sub15_3,sub15_3v,msq15_3,msq15_3v,
     & qqb_QQbdku,qqb_QQbdku_gvec)
      call storedip_mass(m15_3,m15_3v)
      call dips_mass(4,p,2,5,3,sub25_3,sub25_3v,msq25_3,msq25_3v,
     & qqb_QQbdku,qqb_QQbdku_gvec)
      call storedip_mass(m25_3,m25_3v)
      call dips_mass(5,p,1,5,4,sub15_4,sub15_4v,msq15_4,msq15_4v,
     & qqb_QQbdku,qqb_QQbdku_gvec)
      call storedip_mass(m15_4,m15_4v)
      call dips_mass(6,p,2,5,4,sub25_4,sub25_4v,msq25_4,msq25_4v,
     & qqb_QQbdku,qqb_QQbdku_gvec)
      call storedip_mass(m25_4,m25_4v)

c--final-initial
      call dips_mass(7,p,3,5,1,sub35_1,sub35_1v,msq35_1,msq35_1v,
     & qqb_QQbdku,donothing_gvec)
      call storedip_mass(m35_1,m35_1v)
      call dips_mass(8,p,3,5,2,sub35_2,sub35_2v,msq35_2,msq35_2v,
     & qqb_QQbdku,donothing_gvec)
      call storedip_mass(m35_2,m35_2v)
      call dips_mass(9,p,4,5,1,sub45_1,sub45_1v,msq45_1,msq45_1v,
     & qqb_QQbdku,donothing_gvec)
      call storedip_mass(m45_1,m45_1v)
      call dips_mass(10,p,4,5,2,sub45_2,sub45_2v,msq45_2,msq45_2v,
     & qqb_QQbdku,donothing_gvec)
      call storedip_mass(m45_2,m45_2v)

c--final-final
      call dips_mass(11,p,3,5,4,sub35_4,sub35_4v,msq35_4,msq35_4v,
     & qqb_QQbdku,donothing_gvec)
      call storedip_mass(m35_4,m35_4v)
      call dips_mass(12,p,4,5,3,sub45_3,sub45_3v,msq45_3,msq45_3v,
     & qqb_QQbdku,donothing_gvec)
      call storedip_mass(m45_3,m45_3v)


c      write(6,*) 'gs:sub15_2(gg)',sub15_2(gg)
c      write(6,*) 'gs:sub25_1(gg)',sub25_1(gg)

c      write(6,*) 'gs:sub15_3(gg)',sub15_3(gg)
c      write(6,*) 'gs:sub15_4(gg)',sub15_4(gg)
c      write(6,*) 'gs:sub25_3(gg)',sub25_3(gg)
c      write(6,*) 'gs:sub25_4(gg)',sub25_4(gg)

c      write(6,*) 'gs:sub35_1(qq)',sub35_1(qq)
c      write(6,*) 'gs:sub35_2(qq)',sub35_2(qq)
c      write(6,*) 'gs:sub45_1(qq)',sub45_1(qq)
c      write(6,*) 'gs:sub45_2(qq)',sub45_2(qq)

c      write(6,*) 'gs:sub35_4(qq)',sub35_4(qq)
c      write(6,*) 'gs:sub45_3(qq)',sub45_3(qq)


c--- fill the dipole contributions
      do j=-nf,nf
      do k=-nf,nf

      if ((j == 0) .and. (k==0)) then
      msqc(1,j,k)= (+m15_2(1,j,k)+m15_2(2,j,k))*sub15_2(gg)*xn
     &            +(+m15_2v(1,j,k)+m15_2v(2,j,k))*sub15_2v*xn
      msqc(2,j,k)= (+m25_1(1,j,k)+m25_1(2,j,k))*sub25_1(gg)*xn
     &            +(+m25_1v(1,j,k)+m25_1v(2,j,k))*sub25_1v*xn

      msqc(3,j,k)= (+m15_3(2,j,k)+m15_3(0,j,k))*sub15_3(gg)*xn
     &            +(+m15_3v(2,j,k)+m15_3v(0,j,k))*sub15_3v*xn
      msqc(4,j,k)= (+m25_3(1,j,k)+m25_3(0,j,k))*sub25_3(gg)*xn
     &            +(+m25_3v(1,j,k)+m25_3v(0,j,k))*sub25_3v*xn

      msqc(5,j,k)= (+m15_4(1,j,k)+m15_4(0,j,k))*sub15_4(gg)*xn
     &            +(+m15_4v(1,j,k)+m15_4v(0,j,k))*sub15_4v*xn
      msqc(6,j,k)= (+m25_4(2,j,k)+m25_4(0,j,k))*sub25_4(gg)*xn
     &            +(+m25_4v(2,j,k)+m25_4v(0,j,k))*sub25_4v*xn

      msqc(7,j,k)= (+m35_1(2,j,k)+m35_1(0,j,k))*sub35_1(qq)*xn
      msqc(8,j,k)= (+m35_2(1,j,k)+m35_2(0,j,k))*sub35_2(qq)*xn

      msqc(9,j,k)= (+m45_1(1,j,k)+m45_1(0,j,k))*sub45_1(qq)*xn
      msqc(10,j,k)=(+m45_2(2,j,k)+m45_2(0,j,k))*sub45_2(qq)*xn

      msqc(11,j,k)=-((m35_4(1,j,k)+m35_4(2,j,k))/xn
     &              +(xn+1._dp/xn)*m35_4(0,j,k))*sub35_4(qq)
      msqc(12,j,k)=-((m45_3(1,j,k)+m45_3(2,j,k))/xn
     &              +(xn+1._dp/xn)*m45_3(0,j,k))*sub45_3(qq)

      elseif ((j > 0) .and. (k == -j)) then
      msqc(1,j,k)= -msq15_2(j,k)*sub15_2(qq)/xn
      msqc(2,j,k)= -msq25_1(j,k)*sub25_1(qq)/xn

      msqc(3,j,k)= +msq15_3(j,k)*sub15_3(qq)*(xn-2._dp/xn)
      msqc(4,j,k)= +msq25_3(j,k)*sub25_3(qq)*2._dp/xn
      msqc(5,j,k)= +msq15_4(j,k)*sub15_4(qq)*2._dp/xn
      msqc(6,j,k)= +msq25_4(j,k)*sub25_4(qq)*(xn-2._dp/xn)

      msqc(7,j,k)= +msq35_1(j,k)*sub35_1(qq)*(xn-2._dp/xn)
      msqc(8,j,k)= +msq35_2(j,k)*sub35_2(qq)*2._dp/xn
      msqc(9,j,k)= +msq45_1(j,k)*sub45_1(qq)*2._dp/xn
      msqc(10,j,k)=+msq45_2(j,k)*sub45_2(qq)*(xn-2._dp/xn)

      msqc(11,j,k)=-msq35_4(j,k)*sub35_4(qq)/xn
      msqc(12,j,k)=-msq45_3(j,k)*sub45_3(qq)/xn

      elseif ((j < 0) .and. (k == -j)) then
      msqc(1,j,k)= -msq15_2(j,k)*sub15_2(qq)/xn
      msqc(2,j,k)= -msq25_1(j,k)*sub25_1(qq)/xn

      msqc(3,j,k)= +msq15_3(j,k)*sub15_3(qq)*two/xn
      msqc(4,j,k)= +msq25_3(j,k)*sub25_3(qq)*(xn-two/xn)
      msqc(5,j,k)= +msq15_4(j,k)*sub15_4(qq)*(xn-two/xn)
      msqc(6,j,k)= +msq25_4(j,k)*sub25_4(qq)*two/xn

      msqc(7,j,k)= +msq35_1(j,k)*sub35_1(qq)*two/xn
      msqc(8,j,k)= +msq35_2(j,k)*sub35_2(qq)*(xn-two/xn)
      msqc(9,j,k)= +msq45_1(j,k)*sub45_1(qq)*(xn-two/xn)
      msqc(10,j,k)=+msq45_2(j,k)*sub45_2(qq)*two/xn

      msqc(11,j,k)=-msq35_4(j,k)*sub35_4(qq)/xn
      msqc(12,j,k)=-msq45_3(j,k)*sub45_3(qq)/xn

      elseif ((j == 0) .and. (k .ne. 0)) then
      msqc(1,j,k)= 2._dp*tr*msq15_2(-k,k)*sub15_2(qg)
      msqc(2,j,k)= (aveqg/avegg)*
     &    (sub25_1(gq)*msq25_1(0,0)+sub25_1v*msq25_1v(0,0))

      elseif ((j .ne. 0) .and. (k == 0)) then
      msqc(1,j,k)= (aveqg/avegg)*
     &    (sub15_2(gq)*msq15_2(0,0)+sub15_2v*msq15_2v(0,0))
      msqc(2,j,k)= 2._dp*tr*msq25_1(j,-j)*sub25_1(qg)
      endif

      enddo
      enddo

      return
      end

