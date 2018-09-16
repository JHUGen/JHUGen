      subroutine multichan(x1,x2,x3,x4,pin,pout,wtdip,*)
      implicit none
      include 'types.f'
************************************************************************
*                                                                      *
*     Given a n-parton phase space point, generate an (n+1)-parton     *
*     configuration using a multi-channel procedure based on the       *
*     dipole structure of the process at hand                          *
*                                                                      *
*     x1,x2,x3: random numbers uniform on [0,1]                        *
*     x4: random number uniform on [0,1], used for choosing channel    *
*                                                                      *
*     Author: J.M.Campbell                                             *
*       Date: 19th March 2009                                          *
*                                                                      *
*     Based on routine originally used in:                             *
*                                                                      *
*     "Generalized unitarity at work: first NLO QCD results            *
*      for hadronic W+3 jet production"                                *
*      R. K. Ellis, K. Melnikov and G. Zanderighi                      *
*      arXiv:0901.4101 (2009)                                          *
*                                                                      *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'npart.f'
      include 'ptilde.f'
      integer::maxdip,dipconfig(maxd,3),ichan,ip,jp,kp,nu,m,iskip,ipp,
     & kpp,mp
      real(dp):: x1,x2,x3,x4,pin(mxpart,4),pout(mxpart,4),wtdip,
     & tiny,ranchan,ran2,vtilde,x,u,z,y,phi,sik,pt,v1(4),v2(4),vperp(4),
     & qi(4),qj(4),qk(4),dot,k(4),kt(4),ks(4),kDk,ksDks,ktDp,ksDp,
     & dotpr,xjac,sij,sjk,sim,sjm,qm(4),ptsq
      integer diptype
      integer, parameter:: inin=1, infi=2, fiin=3, fifi=4, phot=5
      parameter(tiny=1.e-15_dp)
     
c--- for checking collinear limits
c      x1=1.e-5_dp
c      x2=0.5_dp

c--- set dipole weight and momenta to zero, in case of early return
      wtdip=zip
      pout(:,:)=zip
      
c--- Step 1: obtain number of dipoles and their configurations
      call dipoleconfig(maxdip,dipconfig)
            
c--- Step 2: choose dipole channel
c---  (note: this is uniform and not chosen adaptively) 
      ranchan=x4*real(maxdip,dp)-tiny
      ichan=int(ranchan)+1
      
      if ((ichan < 1) .or. (ichan > maxdip)) then
        write(6,*) 'Chosen channel is out of bounds: ',ichan
        write(6,*) ' (should be between 1 and',maxdip,')'
        stop
      endif
           
c--- Step 3: determine type of dipole (ii/if/fi/ff)   
      ip=dipconfig(ichan,1)
      jp=dipconfig(ichan,2)
      kp=dipconfig(ichan,3)
c      write(6,*) 'Chosen dipole: ',ip,jp,kp
c      write(6,*) 'Chosen channel ',ichan,': ',ip,jp,kp
      
      if     ((ip <= 2) .and. (kp <= 2)) then
        diptype=inin
      elseif ((ip <= 2) .and. (kp > 2)) then
        diptype=infi
      elseif ((ip > 2) .and. (kp <= 2)) then
        diptype=fiin
      elseif ((ip > 2) .and. (kp > 2)) then
        diptype=fifi
      endif

c--- special notation for identified photon dipole
      if (kp == 0) then
        diptype=phot
      endif

c--- Step 4: generate dipole phase-space according to type

c--- NOTE: the section below still implicitly assumes that the
c---       emitted parton is in the last position, since the 
c---       other momenta are assumed to exist in pin.
c---       Need to modify so that this is not the case.
c---       IN FACT OKAY FOR inin AND FIXED NOW FOR infi

      if    (diptype == inin) then
c---    initial-initial dipoles
        vtilde=half*(one-(one-two*x1)*sqrt(one+four*x1*(one-x1)))
        x=one-x2**2
        phi=two*pi*x3
        
        if (vtilde > one-x) then
          return 1    ! c.f. (5.151) of C.-S.
        endif
        
        sik=two*dot(pin,ip,kp)
c--- catch generation of exceptionally collinear points
        if (abs(sik/(two*pin(ip,4)*pin(kp,4))) .lt. tiny) return 1
        ptsq=sik*vtilde*(one-vtilde-x)/(x+tiny)
        if (ptsq < 0) return 1
        pt=sqrt(ptsq)
        
        call getperp(pin,ip,kp,v1,v2,*999)
        
        do nu=1,4
        vperp(nu)=cos(phi)*v1(nu)+sin(phi)*v2(nu)
        
        qi(nu)=pin(ip,nu)/(x+tiny)
        qj(nu)=-(one-vtilde-x)/(x+tiny)*pin(ip,nu)
     &         -vtilde*pin(kp,nu)
     &         +pt*vperp(nu)
        qk(nu)=pin(kp,nu)
        enddo
        
        wtdip=8._dp*sik/x**2/16._dp/pi**2
        
c--- note: we may end up with initial state momentum qi with
c---       energy larger than that of the beam; that should be
c---       rejected outside of this routine
        
      elseif (diptype == infi) then
c---    initial-final dipoles
        u=half*(one-(one-two*x1)*sqrt(one+four*x1*(one-x1)))
        x=one-x2**2
        phi=two*pi*x3
        
c--- account for dipoles with emitted partons not in the last position
        kpp=kp
        if (kp > jp) kpp=kp-1
        
        sik=two*dot(pin,ip,kpp)
c--- catch generation of exceptionally collinear points
        if (abs(sik/pin(ip,4)/pin(kpp,4)) .lt. tiny) return 1
        ptsq=-sik*u*(one-u)*(one-x)/(x+tiny)
        if (ptsq < 0) return 1
        pt=sqrt(ptsq)
        
        call getperp(pin,ip,kpp,v1,v2,*999)
        
        do nu=1,4
        vperp(nu)=cos(phi)*v1(nu)+sin(phi)*v2(nu)

        qi(nu)=pin(ip,nu)/(x+tiny)
        qj(nu)=-(one-u)*(one-x)/(x+tiny)*pin(ip,nu)
     &         +u*pin(kpp,nu)+pt*vperp(nu)
        qk(nu)=-u*(one-x)/(x+tiny)*pin(ip,nu)
     &         +(one-u)*pin(kpp,nu)-pt*vperp(nu)
        enddo
        
        wtdip=-8._dp*sik/x**2/16._dp/pi**2

c--- note: we may end up with initial state momentum qi with
c---       energy larger than that of the beam; that should be
c---       rejected outside of this routine
        
      elseif (diptype == fiin) then
c---    final-initial dipoles
        z=half*(one-(one-two*x2)*sqrt(one+four*x2*(one-x2)))
        x=one-x1**2
        phi=two*pi*x3
        
        sik=two*dot(pin,ip,kp)
c--- catch generation of exceptionally collinear points
        if (abs(sik/(two*pin(ip,4)*pin(kp,4))) .lt. tiny) return 1
        ptsq=-sik*z*(one-z)*(one-x)/(x+tiny)
        if (ptsq < 0) return 1
        pt=sqrt(ptsq)
        
        call getperp(pin,ip,kp,v1,v2,*999)
        
        do nu=1,4
        vperp(nu)=cos(phi)*v1(nu)+sin(phi)*v2(nu)
        
        qk(nu)=pin(kp,nu)/(x+tiny)
        qi(nu)=-(one-z)*(one-x)/(x+tiny)*pin(kp,nu)
     &         +z*pin(ip,nu)+pt*vperp(nu)
        qj(nu)=-z*(one-x)/(x+tiny)*pin(kp,nu)
     &         +(one-z)*pin(ip,nu)-pt*vperp(nu)
        enddo
        
        wtdip=-8._dp*sik/x**2/16._dp/pi**2

c--- note: we may end up with initial state momentum qi with
c---       energy larger than that of the beam; that should be
c---       rejected outside of this routine
        
      elseif (diptype == fifi) then
c---    final-final dipoles
        z=half*(one-(one-two*x1)*sqrt(one+four*x1*(one-x1)))
        y=x2**2
        phi=two*pi*x3
        
c--- account for dipoles with emitted partons not in the last position
        ipp=ip
        kpp=kp
        if (jp < ip) ipp=jp
        if (jp < kp) kpp=jp
        if (ipp == kpp) then
          write(6,*) 'unanticipated dipole in multichan.f'
          write(6,*) 'ip,jp,kp = ',ip,jp,kp
          stop
        endif

        sik=two*dot(pin,ipp,kpp)
c--- catch generation of exceptionally collinear points
        if (abs(sik/(two*pin(ipp,4)*pin(kpp,4))) .lt. tiny) return 1
        ptsq=sik*z*(one-z)*y
        if (ptsq < 0) return 1
        pt=sqrt(ptsq)
        
        call getperp(pin,ipp,kpp,v1,v2,*999)
        
        do nu=1,4
        vperp(nu)=cos(phi)*v1(nu)+sin(phi)*v2(nu)
        
        qi(nu)=z*pin(ipp,nu)+y*(one-z)*pin(kpp,nu)+pt*vperp(nu)
        qj(nu)=(one-z)*pin(ipp,nu)+y*z*pin(kpp,nu)-pt*vperp(nu)
        qk(nu)=(one-y)*pin(kpp,nu)
        enddo
        
        wtdip=8._dp*sik*(one-y)/16._dp/pi**2

      elseif (diptype == phot) then
c---    identified photon dipoles
        z=x2**2
        y=x1
        phi=two*pi*x3
        
        if (npart == 3) then

c--- work out other spectator
        mp=3+4+5-ip-jp
        
        sim=two*dot(pin,3,4)
        ptsq=sim*y*(one-y)*(1-z)
        if (ptsq < 0) return 1
        pt=sqrt(ptsq)
        
        call getperp(pin,ip,mp,v1,v2,*999)
        
        do nu=1,4
        vperp(nu)=cos(phi)*v1(nu)+sin(phi)*v2(nu)
        
        qi(nu)=z*pin(ip,nu)
        qj(nu)=(one-z)*(one-y)*pin(ip,nu)+y*pin(mp,nu)+pt*vperp(nu)
        qm(nu)=(one-z)*y*pin(ip,nu)+(one-y)*pin(mp,nu)-pt*vperp(nu)
        enddo

        wtdip=two*z*sim/16._dp/pi**2
        
        else
          write(6,*) 'Identified photon dipole not written for npart>3' 
          stop
        endif

      else
c---    error or unimplemented type
        write(6,*) 'Unrecognized dipole type in multichan.f: ',diptype
        stop
      endif

c--- check for accidental NaN
      if (qi(4) .ne. qi(4)) then
        write(6,*) 'multichan: qi NaN, diptype = ',diptype
c        write(6,*) 'z,pt,v1,v2',z,pt,v1,v2
        return 1
      endif
      if (qj(4) .ne. qj(4)) then
        write(6,*) 'multichan: qj NaN, diptype = ',diptype
        return 1
      endif
      if (qk(4) .ne. qk(4)) then
        write(6,*) 'multichan: qk NaN, diptype = ',diptype
        return 1
      endif
      
c--- Step 5: reinstate all momenta
      
cc---    note: For now we assume that emitted partons in the
cc---          dipoles should always be added at the end of the list of
cc---          n-parton momenta.
cc---          This is true for simple cases (e.g. Drell-Yan, Wbb), but
cc---          not for more complicated ones.
c        if (jp .ne. npart+2) then
c          write(6,*) 'Dipole multichannel phase space, multichan.f,'
c          write(6,*) 'assumes that emitted radiation is the last'
c          write(6,*) 'momentum, but this is not the case: jp=',jp
c          stop
c        endif

      if    (diptype == inin) then
c---    initial-initial dipoles
        do nu=1,4
          pout(ip,nu)=qi(nu)
          pout(kp,nu)=qk(nu)
          pout(jp,nu)=qj(nu)
c---       must now perform LT on remaining final state partons
c---        (c.f. relevant section of transform.f)
          k(nu)=-qi(nu)-qk(nu)-qj(nu)
          kt(nu)=-x*qi(nu)-qk(nu)
          ks(nu)=k(nu)+kt(nu)
        enddo
        kDk=dotpr(k,k)
        ksDks=dotpr(ks,ks)
        iskip=0
        do m=3,npart+2
          if (m == jp) then
            iskip=1
          else
            ktDp=kt(4)*pin(m-iskip,4)-kt(1)*pin(m-iskip,1)
     &          -kt(2)*pin(m-iskip,2)-kt(3)*pin(m-iskip,3)
            ksDp=ks(4)*pin(m-iskip,4)-ks(1)*pin(m-iskip,1)
     &          -ks(2)*pin(m-iskip,2)-ks(3)*pin(m-iskip,3)
            do nu=1,4
              pout(m,nu)=pin(m-iskip,nu)-two*ksDp*ks(nu)/ksDks
     &                                  +two*ktDp*k(nu)/kDk
            enddo
          endif
        enddo
        
      elseif (( diptype == infi)
     &    .or. (diptype == fiin)
     &    .or. (diptype == fifi) ) then
c---    all other dipoles are straightforward
        do nu=1,4
          pout(ip,nu)=qi(nu)
          pout(kp,nu)=qk(nu)
          pout(jp,nu)=qj(nu)
        enddo
        iskip=0
        do m=1,npart+2
          if     (m == jp) then
            iskip=1
          elseif ((m .ne. ip) .and. (m .ne. kp)) then
            do nu=1,4
              pout(m,nu)=pin(m-iskip,nu)
            enddo
          endif
        enddo

      elseif (diptype == phot) then
c---    identified photon dipole
        if (npart == 3) then
          do nu=1,4
            pout(1,nu)=pin(1,nu)
            pout(2,nu)=pin(2,nu)
            pout(ip,nu)=qi(nu)
            pout(jp,nu)=qj(nu)
            pout(mp,nu)=qm(nu)
          enddo
        else
          write(6,*) 'Should not get here, npart=',npart
          stop
        endif

      else
c---    error or unimplemented type
        write(6,*) 'Unrecognized dipole type in multichan.f: ',diptype
        stop
      endif


c--- Step 6: compute the sum of inverse Jacobians for all channels
      xjac=zip
      
      do m=1,maxdip

      ip=dipconfig(m,1)
      jp=dipconfig(m,2)
      kp=dipconfig(m,3)
      if     ((ip <= 2) .and. (kp <= 2)) then
        diptype=inin
      elseif ((ip <= 2) .and. (kp > 2)) then
        diptype=infi
      elseif ((ip > 2) .and. (kp <= 2)) then
        diptype=fiin
      elseif ((ip > 2) .and. (kp > 2)) then
        diptype=fifi
      endif
c--- special notation for identified photon dipole
      if (kp == 0) then
        diptype=phot
      endif

c--- compute relevant invariants
      if (diptype .ne. phot) then
      sij=two*dot(pout,ip,jp)
      sik=two*dot(pout,ip,kp)
      sjk=two*dot(pout,jp,kp)
      else
      mp=3+4+5-ip-jp
      sij=two*dot(pout,ip,jp)
      sim=two*dot(pout,ip,mp)
      sjm=two*dot(pout,jp,mp)
      endif

      if     (diptype == inin) then
c---    initial-initial dipoles
        x=one+(sij+sjk)/sik
        vtilde=-sij/sik
        xjac=xjac+(sqrt(vtilde)+sqrt(one-vtilde))
     &            /(sqrt(vtilde*(one-vtilde)*(one-x))+tiny)
      elseif (diptype == infi) then
c---    initial-final dipoles
        x=one+sjk/(sij+sik)
        u=sij/(sij+sik)
        xjac=xjac+(sqrt(u)+sqrt(one-u))
     &            /(sqrt(u*(one-u)*(one-x))+tiny)

      elseif (diptype == fiin) then
c---    final-initial dipoles
        x=one+sij/(sjk+sik)
        z=sik/(sik+sjk)
        xjac=xjac+(sqrt(z)+sqrt(one-z))
     &            /(sqrt(z*(one-z)*(one-x))+tiny)

      elseif (diptype == fifi) then
c---    final-final dipoles
        y=sij/(sij+sjk+sik)
        z=sik/(sjk+sik)
        xjac=xjac+(sqrt(z)+sqrt(one-z))
     &            /(sqrt(z*(one-z)*y)+tiny)

      elseif (diptype == phot) then
c---    identified photon dipoles
        z=(sij+sim)/(sij+sim+sjm)
        xjac=xjac+one/(sqrt(z)+tiny)

      else
c---    error or unimplemented type
        write(6,*) 'Unrecognized dipole type in multichan.f: ',diptype
        stop
      endif

      enddo
      
c--- Step 7: compute final weight
      wtdip=wtdip*real(maxdip,dp)/xjac

c--- catch NaN, typically generated by rounding errors in generation
c--- of momenta causing p^2<0
      if (wtdip /= wtdip) then
        wtdip=zip
        return 1
      endif

c--- for checking
c      if (ichan == 8) then
c      write(6,*) 'wtdip = ',wtdip
c      write(6,*) 'PIN'
c      call writeout(pin)
c      write(6,*)
c      write(6,*) 'POUT'
c      call writeout(pout)
c      pause
c      endif
            
      return

c--- alternative return to reject point
  999 return 1
  
      
      end
      
