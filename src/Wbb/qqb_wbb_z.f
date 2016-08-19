      subroutine qqb_wbb_z(p,z)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      integer is
      double precision z,p(mxpart,4),dot
      double precision xl12,xl15,xl16,xl25,xl26,xl56
      double precision ii_qq,ii_qg,if_qq,fi_qq,ff_qq,tempqg

      xl12=dlog(+two*dot(p,1,2)/musq)
      xl15=dlog(-two*dot(p,1,5)/musq)
      xl16=dlog(-two*dot(p,1,6)/musq)
      xl25=dlog(-two*dot(p,2,5)/musq)
      xl26=dlog(-two*dot(p,2,6)/musq)
      xl56=dlog(+two*dot(p,5,6)/musq)

      do is=1,3
      Q1(q,q,a,is)=ason4pi
     .             *((xn-two/xn)*(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
     .                  +two/xn *(if_qq(z,xl16,is)+fi_qq(z,xl16,is))
     .                  -one/xn *(ii_qq(z,xl12,is)+ff_qq(z,xl56,is)))
      Q1(a,a,q,is)=ason4pi
     .             *((xn-two/xn)*(if_qq(z,xl25,is)+fi_qq(z,xl25,is))
     .                  +two/xn *(if_qq(z,xl26,is)+fi_qq(z,xl26,is))
     .                  -one/xn *(ii_qq(z,xl12,is)+ff_qq(z,xl56,is)))
      Q2(a,a,q,is)=ason4pi
     .             *((xn-two/xn)*(if_qq(z,xl26,is)+fi_qq(z,xl26,is))
     .                  +two/xn *(if_qq(z,xl25,is)+fi_qq(z,xl25,is))
     .                  -one/xn *(ii_qq(z,xl12,is)+ff_qq(z,xl56,is)))
      Q2(q,q,a,is)=ason4pi
     .             *((xn-two/xn)*(if_qq(z,xl16,is)+fi_qq(z,xl16,is))
     .                  +two/xn *(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
     .                  -one/xn *(ii_qq(z,xl12,is)+ff_qq(z,xl56,is)))

      tempqg=ason2pi*tr*ii_qg(z,xl12,is)
      Q1(a,g,q,is)=tempqg
      Q1(q,g,a,is)=tempqg
      Q2(a,g,q,is)=tempqg
      Q2(q,g,a,is)=tempqg
      enddo

      return
      end

