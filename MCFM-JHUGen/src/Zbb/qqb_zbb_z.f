      subroutine qqb_zbb_z(p,z)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'agq.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'PR_cs_new.f'
      integer:: is
      real(dp):: z,p(mxpart,4),dot,tempqg,tempgq
      real(dp):: xl12,xl15,xl16,xl25,xl26,xl56
      real(dp):: ii_qq,ii_qg,ii_gq,if_qq,fi_qq,
     &                 ff_qq,ii_gg,if_gg

      xl12=log(+two*dot(p,1,2)/musq)
      xl15=log(-two*dot(p,1,5)/musq)
      xl16=log(-two*dot(p,1,6)/musq)
      xl25=log(-two*dot(p,2,5)/musq)
      xl26=log(-two*dot(p,2,6)/musq)
      xl56=log(+two*dot(p,5,6)/musq)

      do is=1,3
c--- pieces for the 4Q matrix elements
      Q1(q,q,a,is)=ason4pi
     &     *((xn-two/xn)*(if_qq(z,xl16,is)+fi_qq(z,xl16,is))
     &          +two/xn *(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
     &          -one/xn *(ii_qq(z,xl12,is)+ff_qq(z,xl56,is)))
      Q1(a,a,q,is)=ason4pi
     &     *((xn-two/xn)*(if_qq(z,xl26,is)+fi_qq(z,xl26,is))
     &          +two/xn *(if_qq(z,xl25,is)+fi_qq(z,xl25,is))
     &          -one/xn *(ii_qq(z,xl12,is)+ff_qq(z,xl56,is)))
      Q2(a,a,q,is)=ason4pi
     &     *((xn-two/xn)*(if_qq(z,xl25,is)+fi_qq(z,xl25,is))
     &          +two/xn *(if_qq(z,xl26,is)+fi_qq(z,xl26,is))
     &          -one/xn *(ii_qq(z,xl12,is)+ff_qq(z,xl56,is)))
      Q2(q,q,a,is)=ason4pi
     &     *((xn-two/xn)*(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
     &          +two/xn *(if_qq(z,xl16,is)+fi_qq(z,xl16,is))
     &          -one/xn *(ii_qq(z,xl12,is)+ff_qq(z,xl56,is)))

      tempqg=ason2pi*tr*ii_qg(z,xl12,is)
      Q1(a,g,q,is)=tempqg
      Q1(q,g,a,is)=tempqg
      Q2(a,g,q,is)=tempqg
      Q2(q,g,a,is)=tempqg

      tempgq=ason4pi*two*cf*ii_gq(z,xl12,is)
      Q1(g,q,g,is)=tempgq
      Q1(g,a,g,is)=tempgq
      Q2(g,q,g,is)=tempgq
      Q2(g,a,g,is)=tempgq

c--- pieces for the 2Q matrix elements (i.e. gg sub-process)
c--- note that these are separated by colour structure
      R1(g,g,g,0,is)=ason4pi*xn*(if_gg(z,xl16,is)+if_gg(z,xl15,is)
     &                          +fi_qq(z,xl16,is)+fi_qq(z,xl15,is)
     &                          -ff_qq(z,xl56,is))
     &              -ason4pi/xn*(ff_qq(z,xl56,is))
      R1(g,g,g,1,is)=ason4pi*xn*(if_gg(z,xl16,is)+fi_qq(z,xl16,is)
     &                          +ii_gg(z,xl12,is))
     &              -ason4pi/xn*(ff_qq(z,xl56,is))
      R1(g,g,g,2,is)=ason4pi*xn*(if_gg(z,xl15,is)+fi_qq(z,xl15,is)
     &                          +ii_gg(z,xl12,is))
     &              -ason4pi/xn*(ff_qq(z,xl56,is))
      R2(g,g,g,0,is)=ason4pi*xn*(if_gg(z,xl26,is)+if_gg(z,xl25,is)
     &                          +fi_qq(z,xl26,is)+fi_qq(z,xl25,is)
     &                          -ff_qq(z,xl56,is))
     &              -ason4pi/xn*(ff_qq(z,xl56,is))
      R2(g,g,g,1,is)=ason4pi*xn*(if_gg(z,xl25,is)+fi_qq(z,xl25,is)
     &                          +ii_gg(z,xl12,is))
     &              -ason4pi/xn*(ff_qq(z,xl56,is))
      R2(g,g,g,2,is)=ason4pi*xn*(if_gg(z,xl26,is)+fi_qq(z,xl26,is)
     &                          +ii_gg(z,xl12,is))
     &              -ason4pi/xn*(ff_qq(z,xl56,is))
      enddo

      return
      end

