      subroutine bq_tpq_z(p,z)
      implicit none
      include 'types.f'

C     "_z pieces"
c     Matrix element for t-bbar production
c      b(-p1)+u(-p2)-->n(p3)+e^+(p4)+b(p5)+d(p6)
C     averaged(summed) over initial(final) colours and spins

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'PR_stop.f'
      include 'agq.f'
      include 'nwz.f'
      integer:: is
      real(dp):: z,p(mxpart,4),dot,if_qq,fi_qq,if_mqq,fi_mqq,
     & xl16,xl25,mbar25,xl26,xl15,mbar15,xl12,
     & ii_qg,tempqg


      xl12=log(+two*dot(p,1,2)/musq)

      xl16=log(-two*dot(p,1,6)/musq)
      xl25=log((-two*dot(p,1,6)+mt**2)/musq)
      mbar25=mt/sqrt(-two*dot(p,1,6)+mt**2)

      xl26=log(-two*dot(p,2,6)/musq)
      xl15=log((-two*dot(p,2,6)+mt**2)/musq)
      mbar15=mt/sqrt(-two*dot(p,2,6)+mt**2)

c----contributions for one leg

      do is=1,3

      if     (nwz == +1) then
c--- for the case of a b in the initial state (nwz=+1)

c--- ub
      B1(q,q,b,is)=+ason2pi*cf*(
     & if_qq(z,xl16,is)+fi_qq(z,xl16,is))
      B2(b,b,q,is)=+ason2pi*cf*(
     & if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))

c--- ubarb
      Q1(a,a,q,is)=+ason2pi*cf*(
     & if_qq(z,xl16,is)+fi_qq(z,xl16,is))
      Q2(q,q,a,is)=+ason2pi*cf*(
     & if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))

c--- bu
      B1(b,b,q,is)=+ason2pi*cf*(
     & if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
      B2(q,q,b,is)=+ason2pi*cf*(
     & if_qq(z,xl26,is)+fi_qq(z,xl26,is))

c--- bubar
      Q1(q,q,a,is)=+ason2pi*cf*(
     & if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
      Q2(a,a,q,is)=+ason2pi*cf*(
     & if_qq(z,xl26,is)+fi_qq(z,xl26,is))

      tempqg=+ason2pi*tr*ii_qg(z,xl12,is)

c--- ug - will be set by bg below
c      Q2(q,g,q,is)=tempqg

c--- bg
      Q2(q,g,q,is)=tempqg
      Q2(a,g,q,is)=tempqg

c--- ubarg
      Q2(q,g,a,is)=tempqg

c--- gb
      Q1(q,g,q,is)=tempqg
      Q1(a,g,q,is)=tempqg

c--- gu - already set by gb above
c      Q1(q,g,q,is)=tempqg

c--- gubar
      Q1(q,g,a,is)=tempqg

      elseif (nwz == -1) then
c--- for the case of a b~ in the initial state (nwz=-1)

c--- ubbar
      Q1(q,q,a,is)=+ason2pi*cf*(
     & if_qq(z,xl16,is)+fi_qq(z,xl16,is))
      Q2(a,a,q,is)=+ason2pi*cf*(
     & if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))

c--- ubarbbar
      B1(a,a,b,is)=+ason2pi*cf*(
     & if_qq(z,xl16,is)+fi_qq(z,xl16,is))
      B2(b,b,a,is)=+ason2pi*cf*(
     & if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))

c--- bbaru
      Q1(a,a,q,is)=+ason2pi*cf*(
     & if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
      Q2(q,q,a,is)=+ason2pi*cf*(
     & if_qq(z,xl26,is)+fi_qq(z,xl26,is))

c--- bbarubar
      B1(b,b,a,is)=+ason2pi*cf*(
     & if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
      B2(a,a,b,is)=+ason2pi*cf*(
     & if_qq(z,xl26,is)+fi_qq(z,xl26,is))

      tempqg=+ason2pi*tr*ii_qg(z,xl12,is)

c--- ug
      Q2(a,g,q,is)=tempqg

c--- bbarg
      Q2(q,g,a,is)=tempqg
      Q2(a,g,a,is)=tempqg

c--- ubarg - aleady set by bg above
c      Q2(q,g,a,is)=tempqg

c--- gbbar
      Q1(q,g,a,is)=tempqg
      Q1(a,g,a,is)=tempqg

c--- gu
      Q1(a,g,q,is)=tempqg

c--- gubar - already set by gbbar above
c      Q1(q,g,a,is)=tempqg
      endif

      enddo

      return
      end
