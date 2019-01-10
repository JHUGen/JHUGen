      function isALepton(idpart)
      implicit none
      logical isALepton
      integer idpart
      isALepton =
     & idpart.eq.11 .or.
     & idpart.eq.13 .or.
     & idpart.eq.15
      return
      end

      function isANeutrino(idpart)
      implicit none
      logical isANeutrino
      integer idpart
      isANeutrino =
     & idpart.eq.12 .or.
     & idpart.eq.14 .or.
     & idpart.eq.16
      return
      end

      function isUpTypeQuark(idpart)
      implicit none
      logical isUpTypeQuark
      integer idpart
      isUpTypeQuark =
     & idpart.eq.2 .or.
     & idpart.eq.4 .or.
     & idpart.eq.6
      return
      end

      function isUpTypeLightQuark(idpart)
      implicit none
      logical isUpTypeLightQuark
      integer idpart
      isUpTypeLightQuark =
     & idpart.eq.2 .or.
     & idpart.eq.4
      return
      end

      function isDnTypeQuark(idpart)
      implicit none
      logical isDnTypeQuark
      integer idpart
      isDnTypeQuark =
     & idpart.eq.1 .or.
     & idpart.eq.3 .or.
     & idpart.eq.5
      return
      end

      function isAQuark(idpart)
      implicit none
      logical isAQuark
      integer idpart
      isAQuark =
     & idpart.gt.0 .and.
     & idpart.le.6
      return
      end

      function isALightQuark(idpart)
      implicit none
      logical isALightQuark
      integer idpart
      isALightQuark =
     & idpart.gt.0 .and.
     & idpart.lt.6
      return
      end

      function isAGluon(idpart)
      implicit none
      logical isAGluon
      integer idpart
      isAGluon = idpart.eq.21
      return
      end

      function isAnUnknownJet(idpart)
      implicit none
      logical isAnUnknownJet
      integer idpart
      isAnUnknownJet = idpart.eq.0
      return
      end

      function isAJet(idpart)
      implicit none
      logical isAJet
      logical isAGluon
      integer idpart
      isAJet = abs(idpart).le.6 .or. isAGluon(idpart)
      return
      end


