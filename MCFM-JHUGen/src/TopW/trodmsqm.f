      double precision function trodmsqm(j1,j2,j3,j4,j5,j6,j7,p,tmass,
     & mtop,manti)
      implicit none
C---  Author R.K. Ellis: April 2012
C---  matrix element squared summed over colors and spins
      include 'constants.f'
      include 'masses.f'
C-----The dot and spinor products transferred are those derived from q
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,j,k,hb,hc,hg
      double complex qcda(2,2,2),qcdb(2,2,2),qedi(2,2,2),qedf(2,2,2),
     & adk(2,2,2),bdk(2,2,2),idk(2,2,2),fdk(2,2,2),mtop(2,2),manti(2,2)
      double precision prop,tmass,p(mxpart,4),dot
c---calculate the W propagator
      prop=((2d0*dot(p,9,10)-wmass**2)**2+(wmass*wwidth)**2)
C---These two calls exploit the symmetry under 1<->2,3<->4,6<->7,za<->zb 
C---and overall sign change
      call Wbb(j1,j2,j3,j4,j5,j6,j7,tmass,za,zb,1,qedi,qedf,qcda,qcdb)
      call Wbb(j2,j1,j4,j3,j5,j7,j6,tmass,zb,za,2,qedi,qedf,qcda,qcdb)

C---Attach decays to amplitudes
      adk(:,:,:)=czip
      bdk(:,:,:)=czip
      idk(:,:,:)=czip
      fdk(:,:,:)=czip

      do hb=1,2
      do hc=1,2
      do hg=1,2
      do j=1,2
      do k=1,2
      adk(hb,hc,hg)=adk(hb,hc,hg)+mtop(hb,j)*qcda(j,k,hg)*manti(k,hc)
      bdk(hb,hc,hg)=bdk(hb,hc,hg)+mtop(hb,j)*qcdb(j,k,hg)*manti(k,hc)
      idk(hb,hc,hg)=idk(hb,hc,hg)+mtop(hb,j)*qedi(j,k,hg)*manti(k,hc)
      fdk(hb,hc,hg)=fdk(hb,hc,hg)+mtop(hb,j)*qedf(j,k,hg)*manti(k,hc)
      enddo
      enddo
      enddo
      enddo
      enddo

C  zero out matrix element
      trodmsqm=zip
      do hb=1,2
      do hc=1,2
      do hg=1,2
      trodmsqm=trodmsqm+
     & V*xn/eight*(cdabs(adk(hb,hc,hg))**2+cdabs(bdk(hb,hc,hg))**2)
     &+V/(eight*xn)*(cdabs(idk(hb,hc,hg))**2+cdabs(fdk(hb,hc,hg))**2
     &-two*(cdabs(idk(hb,hc,hg)+fdk(hb,hc,hg)))**2)
      enddo
      enddo
      enddo

      trodmsqm=trodmsqm/prop
      return
      end
