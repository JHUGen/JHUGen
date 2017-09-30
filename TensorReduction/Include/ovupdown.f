      double precision,parameter:: 
     & up(4,4)=reshape((
     &        /-1d0, 0d0, 0d0, 0d0,
     &         0d0,-1d0, 0d0, 0d0,
     &         0d0, 0d0,-1d0, 0d0,
     &         0d0, 0d0, 0d0,+1d0/),(/4,4/)),
     & dn(4,4)=reshape((
     &       /+1d0, 0d0, 0d0, 0d0,
     &         0d0,+1d0, 0d0, 0d0,
     &         0d0, 0d0,+1d0, 0d0,
     &         0d0, 0d0, 0d0,+1d0/),(/4,4/))

