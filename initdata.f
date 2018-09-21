c******************************************************************************
c*                                                    AXITRA Version 3.0      *
c*                     SUBROUTINE INITDATA                                    *
c*                                                                            *
c*    Initialisation de divers parametres                                     *
c*                                                                            *
c******************************************************************************

      subroutine initdata

      include "parameter.inc"
      include "dimension1.inc"
c++++++++++++
c		 on reordonne les stations par profondeur
c 		 croissante. Les nrup premieres sont au dessus
c		 de la source.
c++++++++++++
      if (nr.gt.1) then
      !do 10 ir1=1,nr-1
      !do 10 ir2=ir1+1,nr
      !if (zr(ir1).gt.zr(ir2)) then
      !xtr=xr(ir1)
      !ytr=yr(ir1)
      !ztr=zr(ir1)
      !xr(ir1)=xr(ir2)
      !yr(ir1)=yr(ir2)
      !zr(ir1)=zr(ir2)
      !xr(ir2)=xtr
      !yr(ir2)=ytr
      !zr(ir2)=ztr
      !endif
!  10  continue
      endif
      do 12 ir=1,nr
 12   if (zr(ir).lt.zs) nrup=ir
c++++++++++++
c                 distances radiales / source
c++++++++++++
      do 1 ir=1,nr
 1    rr(ir)=sqrt((xr(ir)-xs)*(xr(ir)-xs)+(yr(ir)-ys)*(yr(ir)-ys))
      do 13 ir=1,nr
      if (rr(ir).ne.0.) then
      cosr(ir)=(xr(ir)-xs)/rr(ir)
      sinr(ir)=(yr(ir)-ys)/rr(ir)
      else
      cosr(ir)=1.
      sinr(ir)=0.
      endif
 13   continue

      do 2 ic=1,nc
      vs2(ic)=vs(ic)*vs(ic)
      vp2(ic)=vp(ic)*vp(ic)
 2    continue

c++++++++++++
c                  conversion interface -> epaisseur des couches     
c++++++++++++

      if (hc(1).eq.0.) then
      do 3 ic=1,nc-1
      hc(ic)=hc(ic+1)-hc(ic)
 3    continue
      endif

c++++++++++++
c                 isc = indice de la couche contenant la source
c         zsc = profondeur de la source par rapport au sommet de la couche
c++++++++++++

      hh=0.
      isc=1
      zsc=zs
      do 5 ic=1,nc-1
      hh=hc(ic)
      if (zsc.gt.hh) then
      zsc=zsc-hh
      isc=ic+1
      else
      goto 51
      endif
 5    continue
 51   continue

c++++++++++++
c                 irc() = indice de la couche contenant la station
c         zrc() = profondeur de la station par rapport au sommet de la couche
c++++++++++++

      do 6 ir=1,nr
      hh=0.
      irc(ir)=1
      zrc(ir)=zr(ir)
      do 61 ic=1,nc-1
      hh=hc(ic)
      if (zrc(ir).gt.hh) then
      zrc(ir)=zrc(ir)-hh
      irc(ir)=ic+1
      else
      goto 62
      endif
 61   continue
 62   continue
      
  6   continue

      return
      end
