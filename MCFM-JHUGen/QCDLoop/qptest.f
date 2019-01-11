      implicit none
      include 'types.f'
      include 'qpconstants.f'
      complex(kind=qp) qpI1,qpI2,qpI3,qpI4,Ival(-2:0)
      double complex qlI1,qlI2,qlI3,qlI4,Ivald(-2:0)
      real(kind=qp)  p1sq,p2sq,p3sq,p4sq,s12,s23,
     . m1sq,m2sq,m3sq,m4sq,musq 
      integer ep
      call qpinit
      call qlinit
      musq=1.1_qp
      m1sq=3._qp
      m2sq=5._qp
      m3sq=9._qp
      m4sq=1.3_qp
      p1sq=1.5_qp
      p2sq=-1.7_qp
      p3sq=2.3_qp
      p4sq=2.9_qp
      s12=37._qp
      s23=-15.7_qp

c      write(6,*) 'Tadpole'
c      Ivald(-2:0)=complex(0d0,0d0)
c      do ep=-1,0
c      Ival(ep)=qpI1(m1sq,musq,ep)
c      Ivald(ep)=qlI1(dble(m1sq),dble(musq),ep)
c      write(6,*) ep,Ival(ep)/Ivald(ep)
c      enddo

c      write(6,*) 'Bubble'
c      do ep=-1,0
c      Ival(ep)=qpI2(p1sq,m1sq,m2sq,musq,ep)
c      Ivald(ep)=qlI2(dble(p1sq),dble(m1sq),dble(m2sq),dble(musq),ep)
c      write(6,*) ep,Ival(ep)/Ivald(ep)
c      enddo

      write(6,*) 'Triangle'
      do ep=0,0
      Ival(ep)=
     & qpI3(p1sq,p2sq,p3sq,m1sq,m2sq,m3sq,musq,ep)
      Ivald(ep)=qlI3(dble(p1sq),dble(p2sq),dble(p3sq),
     & dble(m1sq),dble(m2sq),dble(m3sq),dble(musq),ep)
      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

c      write(6,*) 'Box'
c      do ep=0,0
c      Ival(ep)=
c     . qpI4(p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,m2sq,m3sq,m4sq,musq,ep)
c      Ivald(ep)=qlI4(dble(p1sq),dble(p2sq),dble(p3sq),dble(p4sq),
c     & dble(s12),dble(s23),
c     & dble(m1sq),dble(m2sq),dble(m3sq),dble(m4sq),dble(musq),ep)

c      write(6,*) ep,Ival(ep)/Ivald(ep),Ivald(ep)
c      enddo

      write(6,*)
      write(6,*) 'Divergent integrals'
      write(6,*) 'Triangle 1'
      do ep=-2,0
      Ival(ep)=qpI3(zip,zip,p3sq,zip,zip,zip,musq,ep)      
      Ivald(ep)=qlI3(0d0,0d0,dble(p3sq),0d0,0d0,0d0,dble(musq),ep)
      write(6,*) 
     . ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Triangle 2'
      do ep=-2,0
      Ival(ep)=qpI3(zip,p2sq,p3sq,zip,zip,zip,musq,ep)      
      Ivald(ep)=qlI3(0d0,dble(p2sq),dble(p3sq),
     & 0d0,0d0,0d0,dble(musq),ep)
      write(6,*) 
     . ep,Ival(ep)/Ivald(ep)
      enddo


      write(6,*) 'Triangle 3'
      do ep=-2,0
      Ival(ep)=qpI3(zip,p2sq,p3sq,zip,zip,m3sq,musq,ep)      
      Ivald(ep)=qlI3(0d0,dble(p2sq),dble(p3sq),
     & 0d0,0d0,dble(m3sq),dble(musq),ep)
      write(6,*) 
     . ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Triangle 4'
      do ep=-2,0
      Ival(ep)=qpI3(zip,p2sq,m3sq,zip,zip,m3sq,musq,ep)      
      Ivald(ep)=qlI3(0d0,dble(p2sq),dble(m3sq),
     & 0d0,0d0,dble(m3sq),dble(musq),ep)
      write(6,*) 
     . ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Triangle 5'
      do ep=-2,0
      Ival(ep)=qpI3(zip,p2sq,m3sq,zip,zip,m3sq,musq,ep)      
      Ivald(ep)=qlI3(0d0,dble(p2sq),dble(m3sq),
     & 0d0,0d0,dble(m3sq),dble(musq),ep)
      write(6,*) 
     . ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Triangle 6'
      do ep=-2,0
      Ival(ep)=qpI3(m2sq,zip,m3sq,zip,m2sq,m3sq,musq,ep)      
      Ivald(ep)=qlI3(dble(m2sq),0d0,dble(m3sq),
     & 0d0,dble(m2sq),dble(m3sq),dble(musq),ep)
      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box1'
      do ep=-2,0
      Ival(ep)=
     .qpI4(zip,zip,zip,zip,s12,s23,zip,zip,zip,zip,musq,ep)
      Ivald(ep)=qlI4(0d0,0d0,0d0,0d0,
     & dble(s12),dble(s23),
     & 0d0,0d0,0d0,0d0,dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box2'
      do ep=-2,0
      Ival(ep)=
     .qpI4(zip,zip,zip,p4sq,s12,s23,zip,zip,zip,zip,musq,ep)
      Ivald(ep)=qlI4(0d0,0d0,0d0,dble(p4sq),
     & dble(s12),dble(s23),
     & 0d0,0d0,0d0,0d0,dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box3'
      do ep=-2,0
      Ival(ep)=
     .qpI4(zip,zip,zip,p4sq,s12,s23,zip,zip,zip,zip,musq,ep)
      Ivald(ep)=qlI4(0d0,0d0,0d0,dble(p4sq),
     & dble(s12),dble(s23),
     & 0d0,0d0,0d0,0d0,dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box4'
      do ep=-2,0
      Ival(ep)=
     .qpI4(zip,p2sq,zip,p4sq,s12,s23,zip,zip,zip,zip,musq,ep)
      Ivald(ep)=qlI4(0d0,dble(p2sq),0d0,dble(p4sq),
     & dble(s12),dble(s23),
     & 0d0,0d0,0d0,0d0,dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box5'
      do ep=-2,0
      Ival(ep)=
     .qpI4(zip,p2sq,p3sq,p4sq,s12,s23,zip,zip,zip,zip,musq,ep)
      Ivald(ep)=qlI4(0d0,dble(p2sq),dble(p3sq),dble(p4sq),
     & dble(s12),dble(s23),
     & 0d0,0d0,0d0,0d0,dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box6'
      do ep=-2,0
      Ival(ep)=
     .qpI4(zip,zip,m4sq,m4sq,s12,s23,zip,zip,zip,m4sq,musq,ep)
      Ivald(ep)=qlI4(0d0,0d0,dble(m4sq),dble(m4sq),
     & dble(s12),dble(s23),
     & 0d0,0d0,0d0,dble(m4sq),dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box7'
      do ep=-2,0
      Ival(ep)=
     .qpI4(zip,zip,m4sq,p4sq,s12,s23,zip,zip,zip,m4sq,musq,ep)
      Ivald(ep)=qlI4(0d0,0d0,dble(m4sq),dble(p4sq),
     & dble(s12),dble(s23),
     & 0d0,0d0,0d0,dble(m4sq),dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box8'
      do ep=-2,0
      Ival(ep)=
     .qpI4(zip,zip,p3sq,p4sq,s12,s23,zip,zip,zip,m4sq,musq,ep)
      Ivald(ep)=qlI4(0d0,0d0,dble(p3sq),dble(p4sq),
     & dble(s12),dble(s23),
     & 0d0,0d0,0d0,dble(m4sq),dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box9'
      do ep=-2,0
      Ival(ep)=
     .qpI4(zip,p2sq,p3sq,m4sq,s12,s23,zip,zip,zip,m4sq,musq,ep)
      Ivald(ep)=qlI4(0d0,dble(p2sq),dble(p3sq),dble(m4sq),
     & dble(s12),dble(s23),
     & 0d0,0d0,0d0,dble(m4sq),dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box10'
      do ep=-2,0
      Ival(ep)=
     .qpI4(zip,p2sq,p3sq,p4sq,s12,s23,zip,zip,zip,m4sq,musq,ep)
      Ivald(ep)=qlI4(0d0,dble(p2sq),dble(p3sq),dble(p4sq),
     & dble(s12),dble(s23),
     & 0d0,0d0,0d0,dble(m4sq),dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box11'
      do ep=-2,0
      Ival(ep)=
     .qpI4(zip,p2sq,p3sq,p4sq,s12,s23,zip,zip,zip,m4sq,musq,ep)
      Ivald(ep)=qlI4(0d0,dble(p2sq),dble(p3sq),dble(p4sq),
     & dble(s12),dble(s23),
     & 0d0,0d0,0d0,dble(m4sq),dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box12'
      do ep=-2,0
      Ival(ep)=
     .qpI4(zip,m3sq,p3sq,p4sq,s12,s23,zip,zip,m3sq,m4sq,musq,ep)
      Ivald(ep)=qlI4(0d0,dble(m3sq),dble(p3sq),dble(p4sq),
     & dble(s12),dble(s23),
     & 0d0,0d0,dble(m3sq),dble(m4sq),dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box13'
      do ep=-2,0
      Ival(ep)=
     .qpI4(zip,p2sq,p3sq,p4sq,s12,s23,zip,zip,m3sq,m4sq,musq,ep)
      Ivald(ep)=qlI4(0d0,dble(p2sq),dble(p3sq),dble(p4sq),
     & dble(s12),dble(s23),
     & 0d0,0d0,dble(m3sq),dble(m4sq),dble(musq),ep)


      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box14'
      do ep=-2,0
      Ival(ep)=
     .qpI4(zip,p2sq,p3sq,p4sq,s12,s23,zip,zip,zip,m4sq,musq,ep)
      Ivald(ep)=qlI4(0d0,dble(p2sq),dble(p3sq),dble(p4sq),
     & dble(s12),dble(s23),
     & 0d0,0d0,0d0,dble(m4sq),dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box15'
      do ep=-2,0
      Ival(ep)=
     .qpI4(m2sq,p2sq,p3sq,m4sq,s12,s23,zip,m2sq,zip,m4sq,musq,ep)
      Ivald(ep)=qlI4(dble(m2sq),dble(p2sq),dble(p3sq),dble(m4sq),
     & dble(s12),dble(s23),
     & 0d0,dble(m2sq),0d0,dble(m4sq),dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo

      write(6,*) 'Box16'
      do ep=-2,0
      Ival(ep)=
     . qpI4(m2sq,p2sq,p3sq,m4sq,s12,s23,zip,m2sq,m3sq,m4sq,musq,ep)
      Ivald(ep)=qlI4(dble(m2sq),dble(p2sq),dble(p3sq),dble(m4sq),
     & dble(s12),dble(s23),
     & 0d0,dble(m2sq),dble(m3sq),dble(m4sq),dble(musq),ep)

      write(6,*) ep,Ival(ep)/Ivald(ep)
      enddo




      stop
      end
