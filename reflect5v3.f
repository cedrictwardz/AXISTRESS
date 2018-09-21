c***********************************************************
c*                                      AXITRA Version 3.0 *
c*              SUBROUTINE REFLECT5V3                      *
c*                                                         *
c*        Calcul des deplacements avec diverses rotations  *
c*        et recombinaisons.                               *
c*        Multiplication par les termes frequentiel        *
c*        et angulaire :                                   *
c*                     u=u*C3(theta)*CFF(omega )           *
c***********************************************************
      
      subroutine reflect5v3

      include "parameter.inc"
      include "dimension1.inc"
      include "dimension2.inc"

      real*8 cor,cor2,co2r
      complex*16 urx,utx,uzx,ury,uty,uzy,urz,utz,uzz

      do 10 ir=1,nr

      u(ir,1)=u(ir,1)*cff(isc)
      u(ir,2)=u(ir,2)*cff(isc)
      u(ir,3)=u(ir,3)*cff(isc)
      u(ir,4)=u(ir,4)*cff(isc)
      u(ir,5)=u(ir,5)*cff(isc)
 10   continue

c+++++++++++++	
c	Deplacements dus aux forces Fx, Fy, Fz
c+++++++++++++

      do 20 ir=1,nr

      cor=cosr(ir)
      sir=sinr(ir)

c	Fx
      ur(ir,1)=  cor * u(ir,1)
      ut(ir,1)= -sir * u(ir,2)
      uz(ir,1)=  cor * u(ir,3)

c	Fy
      ur(ir,2)=  sir * u(ir,1)
      ut(ir,2)=  cor * u(ir,2)
      uz(ir,2)=  sir * u(ir,3)

c	Fz
      ur(ir,3)= u(ir,4)
      ut(ir,3)= 0.
      uz(ir,3)= u(ir,5)

 20   continue
      
      return
      end


