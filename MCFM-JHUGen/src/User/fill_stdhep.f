************************************************************************
*  Routine to fill the STDHEP common block with event information      *
*  Inputs are: p  - momenta of the particles                           *
*             j,k - identity of incoming partons                       *
*          wgt_jk - the weight corresponding to these incoming partons *
************************************************************************
      subroutine fill_stdhep(p,ij,ik,wgt_jk)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'npart.f'
      include 'stdhep.f'
      include 'masses.f'
      include 'plabel.f'
      integer:: ij,ik,n,nu,eventnumber,jetlabel_to_stdhep
      real(dp):: p(mxpart,4),wgt_jk,mass
      data eventnumber/0/
      save eventnumber

c--- increase the event number by 1
      eventnumber=eventnumber+1

c--- nevhep is the STDHEP event number
      nevhep=eventnumber

c--- nhep is the number of entries in this record
c---    number of initial state particles = 2
c---    number of  final  state particles = npart
c---      additional particle for weight  = 1
cc    nhep=npart+3
      nhep=npart+2

c--- isthep is the status code for each particle:
c---    0 - null
c---    1 - final state
c---    2 - intermediate state

c--- we'll set the initial state ones to 0, as well as the
c--- fictitious weight particle

c--- For expandability, we will put the weight in as particle 1
c---  and "STDHEP particle #" = "MCFM particle #" + 1
c---  DSW, 31.07.2002. Use
      isthep(1)=21
      isthep(2)=21
cc    isthep(3)=0
cc    do n=4,nhep
cc      isthep(n)=1
cc    enddo
      do n=3,nhep
         isthep(n)=1
      enddo

c--- idhep is the particle ID number, as per the PDG standard
c--- (g,d,u,s,c,b,t) = (21,1,2,3,4,5,6)
      if (ij == 0) then
cc      idhep(2)=21
        idhep(1)=21
      else
cc      idhep(2)=ij
        idhep(1)=ij
      endif

      if (ik == 0) then
cc      idhep(3)=21
        idhep(2)=21
      else
cc      idhep(3)=ik
        idhep(2)=ik
      endif

c --- DSW, 31.07.2002
c --- Force uubar for testing purposes :
      idhep(1)=2
      idhep(2)=-2

cc    do n=4,nhep
cc      idhep(n)=jetlabel_to_stdhep(plabel(n-1))
cc    enddo
      do n=3,nhep
        idhep(n)=jetlabel_to_stdhep(plabel(n))
      enddo

c--- jmohep, jdahep are the mother and daughter particles
c--- (not relevant here)
      do n=1,nhep
        jmohep(1,n)=0
        jmohep(2,n)=0
        jdahep(1,n)=0
        jdahep(2,n)=0
      enddo

c--- phep(1..4) = (x,y,z,E) momentum 4-vector
c---   phep(5)  = mass

c--- "Particle 1" contains the weight, as noted above
cc    do nu=1,4
cc      phep(nu,1)=0._dp
cc    enddo
cc    phep(5,1)=wgt_jk

c--- note that all matrix elements are calculated for massless
c--- particles, so that there will be a mismatch between p.p
c--- and this assigned mass. The energy is thus re-scaled,
c--- which induces an error. This should be investigated.
cc    do n=2,nhep
cc      mass=0._dp
cc      if     ((plabel(n-1) == 'el').or.(plabel(n-1) == 'ea')) then
cc        mass=mel
cc      elseif ((plabel(n-1) == 'ml').or.(plabel(n-1) == 'ma')) then
cc        mass=mmu
cc      elseif ((plabel(n-1) == 'tl').or.(plabel(n-1) == 'ta')) then
cc        mass=mtau
cc      elseif ((plabel(n-1) == 'bq').or.(plabel(n-1) == 'ba')) then
cc        mass=mb
cc      endif
cc      phep(5,n)=mass
cc      do nu=1,3
cc        phep(nu,n)=p(n-1,nu)
cc      enddo
c--- here's the re-scaling
cc      phep(4,n)=sqrt(p(n-1,4)**2+mass**2)
cc    enddo

      do n=1,nhep
        mass=0._dp
        if     ((plabel(n) == 'el').or.(plabel(n) == 'ea')) then
          mass=mel
        elseif ((plabel(n) == 'ml').or.(plabel(n) == 'ma')) then
          mass=mmu
        elseif ((plabel(n) == 'tl').or.(plabel(n) == 'ta')) then
          mass=mtau
        elseif ((plabel(n) == 'bq').or.(plabel(n) == 'ba')) then
          mass=mb
        endif
        phep(5,n)=mass
        do nu=1,3
          phep(nu,n)=p(n,nu)
        enddo
c--- here's the re-scaling
        phep(4,n)=sqrt(p(n,4)**2+mass**2)
      enddo



c--- vhep(1..4) is vertex information, zero here
      do n=1,nhep
        do nu=1,4
          vhep(nu,n)=0._dp
        enddo
      enddo

c--- write-out common block (to unit 6) for checking, if necessary
c     call write_stdhep(6)
c     pause

      return
      end

