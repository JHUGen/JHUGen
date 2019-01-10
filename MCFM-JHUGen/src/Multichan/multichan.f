      subroutine multichan(x1,x2,x3,x4,pin,pout,wtdip,*)
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
      implicit none
      include 'constants.f'
      include 'npart.f'
      include 'ptilde.f'
      integer maxdip,dipconfig(maxd,3),ichan,ip,jp,kp,nu,m,iskip,kpp,mp
      double precision x1,x2,x3,x4,pin(mxpart,4),pout(mxpart,4),wtdip,
     . tiny,ranchan,ran2,vtilde,x,u,z,y,phi,sik,pt,v1(4),v2(4),vperp(4),
     . qi(4),qj(4),qk(4),dot,k(4),kt(4),ks(4),kDk,ksDks,ktDp,ksDp,
     . dotpr,xjac,sij,sjk,sim,sjm,qm(4)
      character*2 diptype
      parameter(tiny=1d-15)
     
c--- for checking collinear limits
c      x1=1d-5
c      x2=0.5d0

c--- set dipole weight and momenta to zero, in case of early return
      wtdip=0d0
      do m=1,mxpart     
      do nu=1,4     
        pout(m,nu)=0d0
      enddo     
      enddo     
      
c--- Step 1: obtain number of dipoles and their configurations      
      call dipoleconfig(maxdip,dipconfig)
            
c--- Step 2: choose dipole channel
c---  (note: this is uniform and not chosen adaptively) 
      ranchan=x4*dfloat(maxdip)-tiny
      ichan=int(ranchan)+1
      
      if ((ichan .lt. 1) .or. (ichan .gt. maxdip)) then
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
      
      if     ((ip .le. 2) .and. (kp .le. 2)) then
        diptype='ii'
      elseif ((ip .le. 2) .and. (kp .gt. 2)) then
        diptype='if'
      elseif ((ip .gt. 2) .and. (kp .le. 2)) then
        diptype='fi'
      elseif ((ip .gt. 2) .and. (kp .gt. 2)) then
        diptype='ff'
      endif

c--- special notation for identified photon dipole
      if (kp .eq. 0) then
        diptype='ph'
      endif

c--- Step 4: generate dipole phase-space according to type

c--- NOTE: the section below still implicitly assumes that the
c---       emitted parton is in the last position, since the 
c---       other momenta are assumed to exist in pin.
c---       Need to modify so that this is not the case.
c---       IN FACT OKAY FOR 'ii' AND FIXED NOW FOR 'if'

      if    (diptype .eq. 'ii') then
c---    initial-initial dipoles
        vtilde=0.5d0*(1d0-(1d0-2d0*x1)*dsqrt(1d0+4d0*x1*(1d0-x1)))
        x=1d0-x2**2
        phi=2d0*pi*x3
        
        if (vtilde .gt. 1d0-x) then
          return 1    ! c.f. (5.151) of C.-S.
        endif
        
        sik=2d0*dot(pin,ip,kp)
        pt=dsqrt(sik*vtilde*(1d0-vtilde-x)/(x+tiny))
        
        call getperp(pin,ip,kp,v1,v2)
        
        do nu=1,4
        vperp(nu)=dcos(phi)*v1(nu)+dsin(phi)*v2(nu)
        
        qi(nu)=pin(ip,nu)/(x+tiny)
        qj(nu)=-(1d0-vtilde-x)/(x+tiny)*pin(ip,nu)
     .         -vtilde*pin(kp,nu)
     .         +pt*vperp(nu)
        qk(nu)=pin(kp,nu)
        enddo
        
        wtdip=8d0*sik/x**2/16d0/pi**2

c--- note: we may end up with initial state momentum qi with
c---       energy larger than that of the beam; that should be
c---       rejected outside of this routine
        
      elseif (diptype .eq. 'if') then
c---    initial-final dipoles
        u=0.5d0*(1d0-(1d0-2d0*x1)*dsqrt(1d0+4d0*x1*(1d0-x1)))
        x=1d0-x2**2
        phi=2d0*pi*x3
        
c--- account for dipoles with emitted partons not in the last position        
        kpp=kp
        if (kp .gt. jp) kpp=kp-1
        
        sik=2d0*dot(pin,ip,kpp)
        pt=dsqrt(-sik*u*(1d0-u)*(1d0-x)/(x+tiny))
        
        call getperp(pin,ip,kpp,v1,v2)
        
        do nu=1,4
        vperp(nu)=dcos(phi)*v1(nu)+dsin(phi)*v2(nu)

        qi(nu)=pin(ip,nu)/(x+tiny)
        qj(nu)=-(1d0-u)*(1d0-x)/(x+tiny)*pin(ip,nu)
     .         +u*pin(kpp,nu)+pt*vperp(nu)
        qk(nu)=-u*(1d0-x)/(x+tiny)*pin(ip,nu)
     .         +(1d0-u)*pin(kpp,nu)-pt*vperp(nu)
        enddo
        
        wtdip=-8d0*sik/x**2/16d0/pi**2

c--- note: we may end up with initial state momentum qi with
c---       energy larger than that of the beam; that should be
c---       rejected outside of this routine
        
      elseif (diptype .eq. 'fi') then
c---    final-initial dipoles
        z=0.5d0*(1d0-(1d0-2d0*x2)*dsqrt(1d0+4d0*x2*(1d0-x2)))
        x=1d0-x1**2
        phi=2d0*pi*x3
        
        sik=2d0*dot(pin,ip,kp)
        pt=dsqrt(-sik*z*(1d0-z)*(1d0-x)/(x+tiny))
        
        call getperp(pin,ip,kp,v1,v2)
        
        do nu=1,4
        vperp(nu)=dcos(phi)*v1(nu)+dsin(phi)*v2(nu)
        
        qk(nu)=pin(kp,nu)/(x+tiny)
        qi(nu)=-(1d0-z)*(1d0-x)/(x+tiny)*pin(kp,nu)
     .         +z*pin(ip,nu)+pt*vperp(nu)
        qj(nu)=-z*(1d0-x)/(x+tiny)*pin(kp,nu)
     .         +(1d0-z)*pin(ip,nu)-pt*vperp(nu)
        enddo
        
        wtdip=-8d0*sik/x**2/16d0/pi**2

c--- note: we may end up with initial state momentum qi with
c---       energy larger than that of the beam; that should be
c---       rejected outside of this routine
        
      elseif (diptype .eq. 'ff') then
c---    final-final dipoles
        z=0.5d0*(1d0-(1d0-2d0*x1)*dsqrt(1d0+4d0*x1*(1d0-x1)))
        y=x2**2
        phi=2d0*pi*x3
        
        sik=2d0*dot(pin,ip,kp)
        pt=dsqrt(sik*z*(1d0-z)*y)
        
        call getperp(pin,ip,kp,v1,v2)
        
        do nu=1,4
        vperp(nu)=dcos(phi)*v1(nu)+dsin(phi)*v2(nu)
        
        qi(nu)=z*pin(ip,nu)+y*(1d0-z)*pin(kp,nu)+pt*vperp(nu)
        qj(nu)=(1d0-z)*pin(ip,nu)+y*z*pin(kp,nu)-pt*vperp(nu)
        qk(nu)=(1d0-y)*pin(kp,nu)
        enddo
        
        wtdip=8d0*sik*(1d0-y)/16d0/pi**2

      elseif (diptype .eq. 'ph') then
c---    identified photon dipoles
        z=x2**2
        y=x1
        phi=2d0*pi*x3
        
        if (npart .eq. 3) then

c--- work out other spectator
        mp=3+4+5-ip-jp
        
        sim=2d0*dot(pin,3,4)
        pt=dsqrt(sim*y*(1d0-y)*(1-z))
        
        call getperp(pin,ip,mp,v1,v2)
        
        do nu=1,4
        vperp(nu)=dcos(phi)*v1(nu)+dsin(phi)*v2(nu)
        
        qi(nu)=z*pin(ip,nu)
        qj(nu)=(1d0-z)*(1d0-y)*pin(ip,nu)+y*pin(mp,nu)+pt*vperp(nu)
        qm(nu)=(1d0-z)*y*pin(ip,nu)+(1d0-y)*pin(mp,nu)-pt*vperp(nu)
        enddo

        wtdip=2d0*z*sim/16d0/pi**2
        
        else
          write(6,*) 'Identified photon dipole not written for npart>3' 
          stop
        endif

      else
c---    error or unimplemented type
        write(6,*) 'Unrecognized dipole type in multichan.f: ',diptype
        stop
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

      if    (diptype .eq. 'ii') then
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
          if (m .eq. jp) then
            iskip=1
          else
            ktDp=kt(4)*pin(m-iskip,4)-kt(1)*pin(m-iskip,1)
     .          -kt(2)*pin(m-iskip,2)-kt(3)*pin(m-iskip,3)
            ksDp=ks(4)*pin(m-iskip,4)-ks(1)*pin(m-iskip,1)
     .          -ks(2)*pin(m-iskip,2)-ks(3)*pin(m-iskip,3)
            do nu=1,4
              pout(m,nu)=pin(m-iskip,nu)-2d0*ksDp*ks(nu)/ksDks
     .                                  +2d0*ktDp*k(nu)/kDk
            enddo
          endif
        enddo
        
      elseif (( diptype .eq. 'if')
     .    .or. (diptype .eq. 'fi')
     .    .or. (diptype .eq. 'ff') ) then
c---    all other dipoles are straightforward
        do nu=1,4
          pout(ip,nu)=qi(nu)
          pout(kp,nu)=qk(nu)
          pout(jp,nu)=qj(nu)
        enddo
        iskip=0
        do m=1,npart+2
          if     (m .eq. jp) then
            iskip=1
          elseif ((m .ne. ip) .and. (m .ne. kp)) then
            do nu=1,4
              pout(m,nu)=pin(m-iskip,nu)
            enddo
          endif
        enddo

      elseif (diptype .eq. 'ph') then
c---    identified photon dipole
        if (npart .eq. 3) then
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
      xjac=0d0
      
      do m=1,maxdip

      ip=dipconfig(m,1)
      jp=dipconfig(m,2)
      kp=dipconfig(m,3)
      if     ((ip .le. 2) .and. (kp .le. 2)) then
        diptype='ii'
      elseif ((ip .le. 2) .and. (kp .gt. 2)) then
        diptype='if'
      elseif ((ip .gt. 2) .and. (kp .le. 2)) then
        diptype='fi'
      elseif ((ip .gt. 2) .and. (kp .gt. 2)) then
        diptype='ff'
      endif
c--- special notation for identified photon dipole
      if (kp .eq. 0) then
        diptype='ph'
      endif

c--- compute relevant invariants
      if (diptype .ne. 'ph') then
      sij=two*dot(pout,ip,jp)
      sik=two*dot(pout,ip,kp)
      sjk=two*dot(pout,jp,kp)
      else
      mp=3+4+5-ip-jp
      sij=two*dot(pout,ip,jp)
      sim=two*dot(pout,ip,mp)
      sjm=two*dot(pout,jp,mp)
      endif

      if     (diptype .eq. 'ii') then
c---    initial-initial dipoles
        x=1d0+(sij+sjk)/sik
        vtilde=-sij/sik
        xjac=xjac+(dsqrt(vtilde)+dsqrt(1d0-vtilde))
     .            /(dsqrt(vtilde*(1d0-vtilde)*(1d0-x))+tiny)
      elseif (diptype .eq. 'if') then
c---    initial-final dipoles
        x=1d0+sjk/(sij+sik)
        u=sij/(sij+sik)
        xjac=xjac+(dsqrt(u)+dsqrt(1d0-u))
     .            /(dsqrt(u*(1d0-u)*(1d0-x))+tiny)

      elseif (diptype .eq. 'fi') then
c---    final-initial dipoles
        x=1d0+sij/(sjk+sik)
        z=sik/(sik+sjk)
        xjac=xjac+(dsqrt(z)+dsqrt(1d0-z))
     .            /(dsqrt(z*(1d0-z)*(1d0-x))+tiny)

      elseif (diptype .eq. 'ff') then
c---    final-final dipoles
        y=sij/(sij+sjk+sik)
        z=sik/(sjk+sik)
        xjac=xjac+(dsqrt(z)+dsqrt(1d0-z))
     .            /(dsqrt(z*(1d0-z)*y)+tiny)

      elseif (diptype .eq. 'ph') then
c---    identified photon dipoles
        z=(sij+sim)/(sij+sim+sjm)
        xjac=xjac+1d0/(dsqrt(z)+tiny)

      else
c---    error or unimplemented type
        write(6,*) 'Unrecognized dipole type in multichan.f: ',diptype
        stop
      endif

      enddo
      
c--- Step 7: compute final weight
      wtdip=wtdip*dfloat(maxdip)/xjac

c--- for checking
c      if (ichan .eq. 8) then
c      write(6,*) 'wtdip = ',wtdip
c      write(6,*) 'PIN'
c      call writeout(pin)
c      write(6,*)
c      write(6,*) 'POUT'
c      call writeout(pout)
c      pause
c      endif
            
      return
      end
      
