      subroutine qq_tchan_ztq_dk_z(p,z)
      implicit none
      include 'types.f'
      
C     "_z pieces"
c     u(-j1)+b(j2)->Z(j3,j4)+t(nu(j5)+e+(j6)+b(j7))+d(j8)+g(j9)

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
     & ii_qg,tempqg, ps(mxpart,4)



      ps(1,:)=p(1,:)
      ps(2,:)=p(2,:)
      ps(3,:)=p(3,:)
      ps(4,:)=p(4,:)
      ps(5,:)=p(5,:)+p(6,:)+p(7,:)
      ps(6,:)=p(8,:)
      
      xl12=log(+two*dot(ps,1,2)/musq)            
      xl16=log(-two*dot(ps,1,6)/musq)
      xl25=log(-two*dot(ps,2,5)/musq)
      mbar25=mt/sqrt(-two*dot(ps,2,5))

      xl26=log(-two*dot(ps,2,6)/musq)
      xl15=log(-two*dot(ps,1,5)/musq)
      mbar15=mt/sqrt(-two*dot(ps,1,5))
      
c----contributions for one leg

c--- NEED B1, B2 COMMON BLOCKS TO DISTINGUISH BETWEEN ub AND bu CASES
      do is=1,3

      if (nwz == +1) then
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

      elseif (nwz == -1) then

c--- ubarbbar
      B1(a,a,b,is)=+ason2pi*cf*(
     & if_qq(z,xl16,is)+fi_qq(z,xl16,is))
      B2(b,b,a,is)=+ason2pi*cf*(
     & if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
      
c--- ubbar
      Q1(q,q,a,is)=+ason2pi*cf*(
     & if_qq(z,xl16,is)+fi_qq(z,xl16,is))
      Q2(a,a,q,is)=+ason2pi*cf*(
     & if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
      
c--- bbarubar
      B1(b,b,a,is)=+ason2pi*cf*(
     & if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
      B2(a,a,b,is)=+ason2pi*cf*(
     & if_qq(z,xl26,is)+fi_qq(z,xl26,is))
      
c--- bbaru
      Q1(a,a,q,is)=+ason2pi*cf*(
     & if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
      Q2(q,q,a,is)=+ason2pi*cf*(
     & if_qq(z,xl26,is)+fi_qq(z,xl26,is))

      endif

      tempqg=+ason2pi*tr*ii_qg(z,xl12,is)


         
c--- ug/bg
      Q2(q,g,q,is)=tempqg
      Q2(a,g,q,is)=tempqg
      
c--- ubarg/bbarg
      Q2(q,g,a,is)=tempqg
      Q2(a,g,a,is)=tempqg

c--- gu/gb
      Q1(q,g,q,is)=tempqg
      Q1(a,g,q,is)=tempqg
     
c--- gubar/gbbar
      Q1(q,g,a,is)=tempqg
      Q1(a,g,a,is)=tempqg

      enddo

      return
      end
