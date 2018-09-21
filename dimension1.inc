c*****************************************************************************
c*                                                        AXITRA Version 1.0 *
c*                          DIMENSION1                                       *  
c*                                                                           *
c* 	Declaration de variables generales :                                 *
c*		constantes, parametres du modele                             *
c*                                                                           *
c*****************************************************************************

      complex*16 ai,omega,omega2,a1,ur,ut,uz,xlnf
      real*8     cosr,hc,kr,kr2
      logical    tconv

      common /dim1a/ ai,xlnf,pi,pi2,omega,omega2,a1,kr,kr2,uconv,
     &               tconv
      common /dim1b/ xr(nrp),yr(nrp),zr(nrp),cosr(nrp),sinr(nrp)
      common /dim1c/ xs,ys,zs
      common /dim1d/ zrc(nrp),rr(nrp),nc,nr,nrup,isc,zsc,irc(nrp)
      common /dim1e/ hc(ncp),vp(ncp),vs(ncp),vp2(ncp),vs2(ncp),
     &               rho(ncp),qp(ncp),qs(ncp)
      common /dim1f/ ur(nrp,3),ut(nrp,3),uz(nrp,3)

