      subroutine zqqbQQba_msqij(j1,j2,j3,j4,j5,j6,j7,msquubddb)
********************************************************************
* return squared matrix element for                                *
* 0 -> q(p1)+qb(p2)+Q(p3)+Qb(p4)+gamma(p5)+lb(p6)+l(p7)            *
* q and Q have different flavour                                   *
* above is Nagy-Trocsanyi momentum assignment                      *
* NO average, NO identical particle factor included                *                                        
********************************************************************
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'ewcharge.f'
      include 'new_pspace.f'
      include 'ipsgen.f'
      integer qi,qj,j1,j2,j3,j4,j5,j6,j7
      double precision msquubddb(4),msq(2,2),t,ix,iy
      double precision msqAA(2,2),msqAB(2,2),msqBB(2,2)
      integer i,j,k
      integer hq,Qh,lh,hg,ihel
      double complex propQQ,propQL,m70hA(2,2,16),m70hB(2,2,16),Cfifj
      double complex Ai(2,2,2,2,2,2),Bi(2,2,2,2,2,2)
      double complex Aii(2,2,2,2),Bii(2,2,2,2)
c-----call helicity amplitude from both channel
      call xzqqQQa_qq(j1,j2,j3,j4,j5,j6,j7,Ai,Bi)
      call xzqqQQa_ql(j1,j2,j3,j4,j5,j6,j7,Aii,Bii)
c-----propagator
      propQQ=s(6,7)/Dcmplx(s(6,7)-zmass**2,zwidth*zmass)
      propQL=t(5,6,7)/Dcmplx(t(5,6,7)-zmass**2,zwidth*zmass)
c-----helicity amplitudes
      do qi=1,2
      do qj=1,2
      ihel=1
      do hq=1,2
      do Qh=1,2
      do hg=1,2
      do lh=1,2
c--------compensation for (3,4) or (6,7) exchange
         if (hq*Qh.eq.2) then
            ix=-1D0
         else
            ix=1D0
         endif
         if (hq*lh.ne.2) then
            iy=-1D0*dabs(q1)
         else
            iy=1D0*dabs(q1)
         endif
c--------helicity amplitude
         m70hA(qi,qj,ihel) = 
     .        Cfifj(qi,hq,lh,propQQ)*Ai(qi,qj,hq,Qh,hg,lh)
     .      + Cfifj(qj,Qh,lh,propQQ)*Bi(qi,qj,Qh,hq,hg,lh)
         m70hB(qi,qj,ihel) = 
     .      + iy*(  Cfifj(qi,hq,lh,propQL)*Aii(hq,Qh,hg,lh)
     .            + ix*Cfifj(qj,Qh,lh,propQL)*Bii(Qh,hq,hg,lh))
         ihel=ihel+1
      enddo
      enddo
      enddo
      enddo
      ihel=ihel-1
c-----square them up and sum them up
      msq(qi,qj)=zip
      msqAA(qi,qj)=zip
      msqBB(qi,qj)=zip
      msqAB(qi,qj)=zip
      do i=1,16
         msqAA(qi,qj)=msqAA(qi,qj)+cdabs(m70hA(qi,qj,i))**2
         msqBB(qi,qj)=msqBB(qi,qj)+cdabs(m70hB(qi,qj,i))**2
         msqAB(qi,qj)= msqAB(qi,qj)
     .   +2d0*dreal(dconjg(m70hA(qi,qj,i))*m70hB(qi,qj,i))
      enddo
      do i=1,16
         if (ipsgen.eq.1) then
            msq(qi,qj)=msqAA(qi,qj)
         elseif (ipsgen.eq.2) then
            msq(qi,qj)=msqBB(qi,qj)+msqAB(qi,qj)
         else
            write(6,*) 'Parameter ipsgen should be 1 or 2'
            write(6,*) 'ipsgen = ',ipsgen
            stop
         endif
c         msq(qi,qj)=msqAA(qi,qj)+msqBB(qi,qj)+msqAB(qi,qj)
c         if (new_pspace) then
c            msq(qi,qj)=msqBB(qi,qj)+msqAB(qi,qj)
c         else
c            msq(qi,qj)=msqAA(qi,qj)
c         endif
      enddo
c-----include color factor, gauge couplings
c-----no averaging
      msq(qi,qj)=8d0*8d0*esq**3*gsq**2*msq(qi,qj)
      enddo
      enddo
c-----
      msquubddb(1)=msq(2,1)
      msquubddb(2)=msq(1,2)
      msquubddb(3)=msq(2,2)
      msquubddb(4)=msq(1,1)
c-----done
      return
      end


      subroutine zqqbQQba_msqii(j1,j2,j3,j4,j5,j6,j7,msquubuub)
********************************************************************
* return squared matrix element for                                *
* 0 -> q(p1)+qb(p2)+q(p3)+qb(p4)+gamma(p5)+lb(p6)+l(p7)            *
* q and Q have the same flavour                                    *           
* above is Nagy-Trocsanyi momentum assignment                      *
* NO average, NO identical particle factor included                *                                        
********************************************************************
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'ewcharge.f'
      include 'new_pspace.f'
      include 'ipsgen.f'
      integer qi,qj,j1,j2,j3,j4,j5,j6,j7
      double precision msquubuub(2),t,ix,iy,m70hsq(64),msq(2)
      double precision m70hsqAA(64),m70hsqBB(64),m70hsqAB(64)
      integer i,j,k
      integer hq,Qh,lh,hg,ihel
      double complex propQQ,propQL,Cfifj
      double complex m70hA(2,2,2,2,2),m70hB(2,2,2,2,2)
      double complex Ai(2,2,2,2,2,2),Bi(2,2,2,2,2,2)
      double complex Ci(2,2,2,2,2,2),Di(2,2,2,2,2,2)
      double complex Aii(2,2,2,2),Bii(2,2,2,2),Cii(2,2,2,2),Dii(2,2,2,2)
c-----call helicity amplitude from both channel
      call xzqqQQa_qq(j1,j2,j3,j4,j5,j6,j7,Ai,Bi)
      call xzqqQQa_ql(j1,j2,j3,j4,j5,j6,j7,Aii,Bii)
      call xzqqQQa_qq(j1,j4,j3,j2,j5,j6,j7,Ci,Di)
      call xzqqQQa_ql(j1,j4,j3,j2,j5,j6,j7,Cii,Dii)
c-----propagator
      propQQ=s(6,7)/Dcmplx(s(6,7)-zmass**2,zwidth*zmass)
      propQL=t(5,6,7)/Dcmplx(t(5,6,7)-zmass**2,zwidth*zmass)
c-----
      do qi=1,2
      qj=qi
c-----helicity amplitudes
      do hq=1,2
      do Qh=1,2
      do hg=1,2
      do lh=1,2
c--------compensation for (3,4) or (6,7) exchange
         if (hq*Qh.eq.2) then
            ix=-1D0
         else
            ix=1D0
         endif
         if (hq*lh.ne.2) then
            iy=-1D0*dabs(q1)
         else
            iy=1D0*dabs(q1)
         endif
c--------helicity amplitude
      m70hA(1,hq,Qh,hg,lh)=
     .  Cfifj(qi,hq,lh,propQQ)*Ai(qi,qj,hq,Qh,hg,lh)
     . +Cfifj(qj,Qh,lh,propQQ)*Bi(qi,qj,Qh,hq,hg,lh)
      m70hB(1,hq,Qh,hg,lh)=
     .  iy*( Cfifj(qi,hq,lh,propQL)*Aii(hq,Qh,hg,lh)
     .      +ix*Cfifj(qj,Qh,lh,propQL)*Bii(Qh,hq,hg,lh))
      m70hA(2,hq,Qh,hg,lh)=
     . -( Cfifj(qi,hq,lh,propQQ)*Ci(qi,qj,hq,Qh,hg,lh)
     .   +Cfifj(qj,Qh,lh,propQQ)*Di(qi,qj,Qh,hq,hg,lh))
      m70hB(2,hq,Qh,hg,lh)=
     . -iy*( Cfifj(qi,hq,lh,propQL)*Cii(hq,Qh,hg,lh)
     .      +ix*Cfifj(qj,Qh,lh,propQL)*Dii(Qh,hq,hg,lh))
      enddo
      enddo
      enddo
      enddo
c-----square them up 
      ihel=1
      do hq=1,2
      do Qh=1,2
      do hg=1,2
      do lh=1,2
      if (hq.eq.Qh) then
         m70hsqAA(ihel)=
     .   +8D0*( cdabs(m70hA(1,hq,Qh,hg,lh))**2
     .         +cdabs(m70hA(2,hq,Qh,hg,lh))**2)
     .   +(-8D0/3D0)*2D0*dreal( m70hA(1,hq,Qh,hg,lh)
     .                            *Dconjg(m70hA(2,hq,Qh,hg,lh)))
         m70hsqBB(ihel)=
     .   +8D0*( cdabs(m70hB(1,hq,Qh,hg,lh))**2
     .         +cdabs(m70hB(2,hq,Qh,hg,lh))**2)
     .   +(-8D0/3D0)*2D0*dreal( m70hB(1,hq,Qh,hg,lh)
     .                            *Dconjg(m70hB(2,hq,Qh,hg,lh)))
         m70hsqAB(ihel)=
     .   +8D0*( 2D0*dreal( 
     .           Dconjg(m70hA(1,hq,Qh,hg,lh))*m70hB(1,hq,Qh,hg,lh)
     .          +Dconjg(m70hA(2,hq,Qh,hg,lh))*m70hB(2,hq,Qh,hg,lh)))
     .   +(-8D0/3D0)*2D0*dreal( 
     .           Dconjg(m70hA(1,hq,Qh,hg,lh))*m70hB(2,hq,Qh,hg,lh)
     .          +Dconjg(m70hB(1,hq,Qh,hg,lh))*m70hA(2,hq,Qh,hg,lh))
         ihel=ihel+1
      elseif (hq.ne.Qh) then
         m70hsqAA(ihel)=8D0*cdabs(m70hA(1,hq,Qh,hg,lh))**2
         m70hsqBB(ihel)=8D0*cdabs(m70hB(1,hq,Qh,hg,lh))**2
         m70hsqAB(ihel)=8D0*2D0*dreal( Dconjg(m70hA(1,hq,Qh,hg,lh))
     .                                   *m70hB(1,hq,Qh,hg,lh) )
         ihel=ihel+1
         m70hsqAA(ihel)=8D0*cdabs(m70hA(2,hq,Qh,hg,lh))**2
         m70hsqBB(ihel)=8D0*cdabs(m70hB(2,hq,Qh,hg,lh))**2
         m70hsqAB(ihel)=8D0*2D0*dreal( Dconjg(m70hA(2,hq,Qh,hg,lh))
     .                                   *m70hB(2,hq,Qh,hg,lh) )
         ihel=ihel+1
      else
         m70hsqAA(ihel)=zip
         m70hsqBB(ihel)=zip
         m70hsqAB(ihel)=zip
      endif
      enddo
      enddo
      enddo
      enddo
      ihel=ihel-1
      do i=1,24
c         m70hsq(i)=m70hsqAA(i)+m70hsqBB(i)+m70hsqAB(i)
         if     (ipsgen .eq. 1) then
            m70hsq(i)=m70hsqAA(i)
         elseif (ipsgen .eq. 2) then
            m70hsq(i)=m70hsqBB(i)+m70hsqAB(i)
         else
            write(6,*) 'Parameter ipsgen should be 1 or 2'
            write(6,*) 'ipsgen = ',ipsgen
            stop
         endif
c         if (new_pspace) then
c            m70hsq(i)=m70hsqBB(i)+m70hsqAB(i)
c         else
c            m70hsq(i)=m70hsqAA(i)
c         endif
      enddo
c-----sum them up
      msq(qi)=zip
      do i=1,24
         msq(qi)=msq(qi)+m70hsq(i)
      enddo
c-----include color factor, gauge couplings
c-----no averaging
      msq(qi)=8d0*esq**3*gsq**2*msq(qi)
c-----
      enddo
      msquubuub(2)=msq(2)
      msquubuub(1)=msq(1)
c-----done
      return
      end



      double complex function Cfifj(qi,hqi,hl,prop)
*--------------------------------------------------------*
*  function to return gamma*/Z coupling to quark/lepton  *
*--------------------------------------------------------*
      implicit none
      integer qi,hqi,hl
      double complex prop
      include 'constants.f'
      include 'zcouple.f'
      include 'ewcharge.f'
c-----fill them up  
      if (hqi.eq.1.and.hl.eq.1) then
         Cfifj=q1*Q(qi)+l(qi)*l1*prop
      elseif (hqi.eq.1.and.hl.eq.2) then
         Cfifj=q1*Q(qi)+l(qi)*r1*prop
      elseif (hqi.eq.2.and.hl.eq.1) then
         Cfifj=q1*Q(qi)+r(qi)*l1*prop
      elseif (hqi.eq.2.and.hl.eq.2) then
         Cfifj=q1*Q(qi)+r(qi)*r1*prop
      endif
c-----done
      return
      end


