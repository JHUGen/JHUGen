      subroutine w2jetnx(i1,i2,i3,i4,i5,i6,p,n,za,zb,zab,zba)
C----matrix element squared with p5 line contracted with n(mu)
C----nDp6 should be equal to zero
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'mmsqv_cs.f'
      include 'flags.f'
      include 'lc.f'
      double complex qcdabn(2,2,2),qcdban(2,2,2),qedn(2,2,2)
      double complex zab(mxpart,mxpart),zba(mxpart,mxpart)
      double precision msq1,msq2,msqq,n(4),p(mxpart,4)
      double precision nDp5
      integer i1,i2,i3,i4,i5,i6

      mmsqv_cs(:,:,:)=0d0

      nDp5=n(4)*p(i5,4)-n(3)*p(i5,3)-n(2)*p(i5,2)-n(1)*p(i5,1)

      call checkndotp(p,n,i6)

      call subqcdn(i1,i2,i3,i4,i5,i6,nDp5,za,zb,zab,zba,qcdabn,qcdban)
            
C--first argument is quark line
C--second argument is polarization of i5 line
C  1=L,2=R
      qedn(1,1,1)=qcdabn(1,1,1)+qcdban(1,1,1) 
      qedn(2,1,1)=qcdabn(2,1,1)+qcdban(2,1,1) 

      msq1= abs(qcdabn(1,1,1))**2+abs(qcdabn(2,1,1))**2
      msq2= abs(qcdban(1,1,1))**2+abs(qcdban(2,1,1))**2 
      msqq= abs(qedn(1,1,1))**2+abs(qedn(2,1,1))**2

      mmsqv_cs(1,+1,+1)=msq1
      mmsqv_cs(2,+1,+1)=msq2
      if ((Qflag) .and. (colourchoice .eq. 1)) then
        mmsqv_cs(0,+1,+1)=zip
      else
        mmsqv_cs(0,+1,+1)=-ninth*msqq
      endif


      return
      end

