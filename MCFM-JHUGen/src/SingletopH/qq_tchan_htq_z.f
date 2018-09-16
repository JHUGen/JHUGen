      subroutine qq_tchan_htq_z(p,z)
      implicit none
C     "_z pieces"
c     u(-j1)+b(j2)->H(j3,j4)+t(5)+d(j6)+g(j7)

      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'PR_stop.f'
      include 'agq.f'
      include 'nwz.f'
      integer is
      double precision z,p(mxpart,4),dot,if_qq,fi_qq,if_mqq,fi_mqq,
     . xl16,xl25,mbar25,xl26,xl15,mbar15,xl12,
     . ii_qg,tempqg

      
      xl12=dlog(+two*dot(p,1,2)/musq)
            
      xl16=dlog(-two*dot(p,1,6)/musq)
      xl25=dlog(-two*dot(p,2,5)/musq)
      mbar25=mt/dsqrt(-two*dot(p,2,5))

      xl26=dlog(-two*dot(p,2,6)/musq)
      xl15=dlog(-two*dot(p,1,5)/musq)
      mbar15=mt/dsqrt(-two*dot(p,1,5))
      
c----contributions for one leg

      do is=1,3

      if (nwz .eq. +1) then
c--- Counterterms for t
c--- NEED B1, B2 COMMON BLOCKS TO DISTINGUISH BETWEEN ub AND bu CASES
c--- ub
      B1(q,q,b,is)=+ason2pi*cf*(
     . if_qq(z,xl16,is)+fi_qq(z,xl16,is))
      B2(b,b,q,is)=+ason2pi*cf*(
     . if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
      
c--- dbarb
      Q1(a,a,q,is)=+ason2pi*cf*(
     . if_qq(z,xl16,is)+fi_qq(z,xl16,is))
      Q2(q,q,a,is)=+ason2pi*cf*(
     . if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
      
c--- bu
      B1(b,b,q,is)=+ason2pi*cf*(
     . if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
      B2(q,q,b,is)=+ason2pi*cf*(
     . if_qq(z,xl26,is)+fi_qq(z,xl26,is))
      
c--- bdbar
      Q1(q,q,a,is)=+ason2pi*cf*(
     . if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
      Q2(a,a,q,is)=+ason2pi*cf*(
     . if_qq(z,xl26,is)+fi_qq(z,xl26,is))
      
      else      
c--- Counterterms for t~
c--- NEED B1, B2 COMMON BLOCKS TO DISTINGUISH ubarbbar AND bbarubar CASES
c--- dbbar
      Q1(q,q,a,is)=+ason2pi*cf*(
     . if_qq(z,xl16,is)+fi_qq(z,xl16,is))
      Q2(a,a,q,is)=+ason2pi*cf*(
     . if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
      
c--- ubarbbar
      B1(a,a,b,is)=+ason2pi*cf*(
     . if_qq(z,xl16,is)+fi_qq(z,xl16,is))
      B2(b,b,a,is)=+ason2pi*cf*(
     . if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
      
c--- bbard
      Q1(a,a,q,is)=+ason2pi*cf*(
     . if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
      Q2(q,q,a,is)=+ason2pi*cf*(
     . if_qq(z,xl26,is)+fi_qq(z,xl26,is))
      
c--- bbarubar
      B1(b,b,a,is)=+ason2pi*cf*(
     . if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
      B2(a,a,b,is)=+ason2pi*cf*(
     . if_qq(z,xl26,is)+fi_qq(z,xl26,is))
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
