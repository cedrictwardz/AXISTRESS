c**************************************************************************
c*                                                     AXITRA Version 3.0 *
c*                                                                        *
c*                     SUBROUTINE REFLECT4V3                              *
c*                                                                        *
c*  Calcul des potentiels et des deplacements a chaque recepteur.         *
c*  - Matrice de passage des potentiels de la couche source aux couches   *
c*    recepteur (FTUP et FTDO)                                            *
c*  - Calcul des potentiels dans toutes les couches (PU et PD)            *
c*  - Calcul de 11 deplacements (termes intermediaires) a chaque recepteur*
c*    (U)                                                                 *
c*                                                                        *
c**************************************************************************


      subroutine reflect4v3

      include "parameter.inc"
      include "dimension1.inc"
      include "dimension2.inc"

      complex*16   egam,enu,s1phiu,s1phid,pu,pd,push,pdsh,
     &             s1psiu,s1psid,s2phiu,s2phid,s2psiu,s2psid,s3phiu,
     &             s3phid,s3psiu,s3psid,s4phid,s4phiu,s4psid,s4psiu,
     &             s5,s6,ftup,ftdo,ftupsh,ftdosh
      
      dimension pu(nrp,2,4),pd(nrp,2,4),push(nrp,2),pdsh(nrp,2),
     &          ftup(ncp,2,2),ftdo(ncp,2,2),ftupsh(ncp),ftdosh(ncp)


c      
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  Matrices de passage des vecteurs potentiels dans la 
c  couche source aux vecteurs potentiel dans chaque couche
c
c                    [ftup] et [ftdo]
c       
c               ------------------------
c
c  Les vecteurs potentiels pu() et pd() sont obtenus a
c  partir des vecteurs potentiels su() et sd() dans la 
c  couche source par :
c
c  Couche (n) au dessus de la couche source :
c
c   pu(n) = [fup(n)]*[fup(n+1)]* *[fup(isc-1] . su
c
c   d'ou l'on tire pd(n) par  pd(n) = [nt(n)] . pu(n)
c
c  Couche (m) au dessous de la couche source :
c
c   pd(m) = [fdo(m)]*[fdo(m-1)]* *[fdo(isc+1)] . sd
c
c   d'ou l'on tire pu(m) par  pu(m) = [mt(m)] . pd(m)
c
c                -------------------------
c   On pose :
c
c        [ftup(n)] = [fup(n)]*...*[fup(isc-1)]*[tud]
c
c        [ftdo(m)] = [fdo(m)]*...*[fdo(isc+1)]*[tdu]
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c		Couches au dessus de la couche source
c
      ftup(isc,1,1)=1.
      ftup(isc,1,2)=0.
      ftup(isc,2,1)=0.
      ftup(isc,2,2)=1.
      ftupsh(isc)=1.

      do 10 ic=isc-1,1,-1

      ftup(ic,1,1)=fup(ic,1,1)*ftup(ic+1,1,1)+fup(ic,1,2)*ftup(ic+1,2,1)
      ftup(ic,1,2)=fup(ic,1,1)*ftup(ic+1,1,2)+fup(ic,1,2)*ftup(ic+1,2,2)
      ftup(ic,2,1)=fup(ic,2,1)*ftup(ic+1,1,1)+fup(ic,2,2)*ftup(ic+1,2,1)
      ftup(ic,2,2)=fup(ic,2,1)*ftup(ic+1,1,2)+fup(ic,2,2)*ftup(ic+1,2,2)
      ftupsh(ic)=fupsh(ic)*ftupsh(ic+1)

 10   continue

c		Couches au dessous de la couche source
c
      ftdo(isc,1,1)=1.
      ftdo(isc,1,2)=0.
      ftdo(isc,2,1)=0.
      ftdo(isc,2,2)=1.
      ftdosh(isc)=1.

      do 11 ic=isc+1,nc

      ftdo(ic,1,1)=fdo(ic,1,1)*ftdo(ic-1,1,1)+fdo(ic,1,2)*ftdo(ic-1,2,1)
      ftdo(ic,1,2)=fdo(ic,1,1)*ftdo(ic-1,1,2)+fdo(ic,1,2)*ftdo(ic-1,2,2)
      ftdo(ic,2,1)=fdo(ic,2,1)*ftdo(ic-1,1,1)+fdo(ic,2,2)*ftdo(ic-1,2,1)
      ftdo(ic,2,2)=fdo(ic,2,1)*ftdo(ic-1,1,2)+fdo(ic,2,2)*ftdo(ic-1,2,2)
      ftdosh(ic)=fdosh(ic)*ftdosh(ic-1)

 11   continue

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  Vecteurs potentiel montant (pu) et descendant (pd),
c  dans chaque couche recepteur, pour les 6 sources elementaires
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c		Recepteurs au dessus de la source

      do 20 ir=1,nrup

      ic=irc(ir)
      pu(ir,1,1)=ftup(ic,1,1)*su1(1)+ftup(ic,1,2)*su1(2)
      pu(ir,2,1)=ftup(ic,2,1)*su1(1)+ftup(ic,2,2)*su1(2)
      pd(ir,1,1)=nt(ic,1,1)*pu(ir,1,1)+nt(ic,1,2)*pu(ir,2,1)
      pd(ir,2,1)=nt(ic,2,1)*pu(ir,1,1)+nt(ic,2,2)*pu(ir,2,1)

      pu(ir,1,2)=ftup(ic,1,1)*su2(1)+ftup(ic,1,2)*su2(2)
      pu(ir,2,2)=ftup(ic,2,1)*su2(1)+ftup(ic,2,2)*su2(2)
      pd(ir,1,2)=nt(ic,1,1)*pu(ir,1,2)+nt(ic,1,2)*pu(ir,2,2) 
      pd(ir,2,2)=nt(ic,2,1)*pu(ir,1,2)+nt(ic,2,2)*pu(ir,2,2)

      pu(ir,1,3)=ftup(ic,1,1)*su3(1)+ftup(ic,1,2)*su3(2)
      pu(ir,2,3)=ftup(ic,2,1)*su3(1)+ftup(ic,2,2)*su3(2)
      pd(ir,1,3)=nt(ic,1,1)*pu(ir,1,3)+nt(ic,1,2)*pu(ir,2,3) 
      pd(ir,2,3)=nt(ic,2,1)*pu(ir,1,3)+nt(ic,2,2)*pu(ir,2,3)

      pu(ir,1,4)=ftup(ic,1,1)*su4(1)+ftup(ic,1,2)*su4(2)
      pu(ir,2,4)=ftup(ic,2,1)*su4(1)+ftup(ic,2,2)*su4(2)
      pd(ir,1,4)=nt(ic,1,1)*pu(ir,1,4)+nt(ic,1,2)*pu(ir,2,4) 
      pd(ir,2,4)=nt(ic,2,1)*pu(ir,1,4)+nt(ic,2,2)*pu(ir,2,4)

      push(ir,1)=ftupsh(ic)*su1sh
      pdsh(ir,1)=ntsh(ic)*push(ir,1)

      push(ir,2)=ftupsh(ic)*su2sh
      pdsh(ir,2)=ntsh(ic)*push(ir,2)

 20   continue

c             Recepteurs au dessous de la source

      do 21 ir=nrup+1,nr
 
      ic=irc(ir)
 
      pd(ir,1,1)=ftdo(ic,1,1)*sd1(1)+ftdo(ic,1,2)*sd1(2)
      pd(ir,2,1)=ftdo(ic,2,1)*sd1(1)+ftdo(ic,2,2)*sd1(2)
      pu(ir,1,1)=mt(ic,1,1)*pd(ir,1,1)+mt(ic,1,2)*pd(ir,2,1)
      pu(ir,2,1)=mt(ic,2,1)*pd(ir,1,1)+mt(ic,2,2)*pd(ir,2,1)
 
      pd(ir,1,2)=ftdo(ic,1,1)*sd2(1)+ftdo(ic,1,2)*sd2(2)
      pd(ir,2,2)=ftdo(ic,2,1)*sd2(1)+ftdo(ic,2,2)*sd2(2)
      pu(ir,1,2)=mt(ic,1,1)*pd(ir,1,2)+mt(ic,1,2)*pd(ir,2,2)
      pu(ir,2,2)=mt(ic,2,1)*pd(ir,1,2)+mt(ic,2,2)*pd(ir,2,2)
 
      pd(ir,1,3)=ftdo(ic,1,1)*sd3(1)+ftdo(ic,1,2)*sd3(2)
      pd(ir,2,3)=ftdo(ic,2,1)*sd3(1)+ftdo(ic,2,2)*sd3(2)
      pu(ir,1,3)=mt(ic,1,1)*pd(ir,1,3)+mt(ic,1,2)*pd(ir,2,3)
      pu(ir,2,3)=mt(ic,2,1)*pd(ir,1,3)+mt(ic,2,2)*pd(ir,2,3)
 
      pd(ir,1,4)=ftdo(ic,1,1)*sd4(1)+ftdo(ic,1,2)*sd4(2)
      pd(ir,2,4)=ftdo(ic,2,1)*sd4(1)+ftdo(ic,2,2)*sd4(2)
      pu(ir,1,4)=mt(ic,1,1)*pd(ir,1,4)+mt(ic,1,2)*pd(ir,2,4)
      pu(ir,2,4)=mt(ic,2,1)*pd(ir,1,4)+mt(ic,2,2)*pd(ir,2,4)
 
      pdsh(ir,1)=ftdosh(ic)*sd1sh
      push(ir,1)=mtsh(ic)*pdsh(ir,1)

      pdsh(ir,2)=ftdosh(ic)*sd2sh
      push(ir,2)=mtsh(ic)*pdsh(ir,2)
 
 21   continue

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  Deplacements pour chaque sources du tenseur, exprime a l'aide de 
c  de source intermediaires. Chaque source intermediaire correspond
c  aux rayonnements des trois potentiels PHI, PSI, KHI de chaque
c  moment du tenseur.
c
c  ex : Mxy -> PHI0, PSI0, KHI0
c
c              PHI0 -> PHI, PSI apres conversion sur les interfaces
c              PSI0 -> PHI, PSI       "                "
c              KHI0 -> KHI            "                "
c
c                       -------------------------
c
c		u = C2(kr)*C4(kr.r)*C5(kr,z)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c                          Termes : C2(kr)
c                          ---------------
c		Constantes dependant du nombre d'onde et
c		de la couche source

      cs2=c2(isc)
      cs4=kr2
      cs5=ai*cgam(isc)
      cs3=ckb2(isc)/cs5
      cs6=ckb2(isc)
      cs8=ai*cnu(isc)
      cs7=cs8*kr2
      cs9=(ckb2(isc)-2.*kr2)/cs5
      do 30 ir=1,nr

      zc=zrc(ir)
      ic=irc(ir)

c                          Termes : C5(kr,z)
c                          -----------------
 
c		Dephasage par rapport au sommet de la couche pour
c		les ondes PHI, PSI et KHI

      egam=exp(-ai*cgam(ic)*zc)
      enu=exp(-ai*cnu(ic)*zc)
      cr1=ai*cgam(ic)
      cr2=cgam(ic)*cgam(ic)
      cr3=ai*cnu(ic)

c		termes sources 
      s1phiu=pu(ir,1,1)/enu
      s1phid=pd(ir,1,1)*enu
      s1psiu=pu(ir,2,1)/egam
      s1psid=pd(ir,2,1)*egam
      s2phiu=pu(ir,1,2)/enu
      s2phid=pd(ir,1,2)*enu
      s2psiu=pu(ir,2,2)/egam
      s2psid=pd(ir,2,2)*egam
      s3phiu=pu(ir,1,3)/enu
      s3phid=pd(ir,1,3)*enu
      s3psiu=pu(ir,2,3)/egam
      s3psid=pd(ir,2,3)*egam
      s4phiu=pu(ir,1,4)/enu
      s4phid=pd(ir,1,4)*enu
      s4psiu=pu(ir,2,4)/egam
      s4psid=pd(ir,2,4)*egam
      s5=push(ir,1)/egam+pdsh(ir,1)*egam
      s6=push(ir,2)/egam+pdsh(ir,2)*egam
 
c				Source phi
      cz1=(s1phiu+s1phid)+cr1*(s1psiu-s1psid)
      cz1b=cr3*(s1phiu-s1phid)+cs4*(s1psiu+s1psid)
c				Source phi*sign(z-z0)
      cz2=(s2phiu+s2phid)+cr1*(s2psiu-s2psid)
      cz2b=cr3*(s2phiu-s2phid)+cs4*(s2psiu+s2psid)
c				Source psi
      cz3=(s3phiu+s3phid)+cr1*(s3psiu-s3psid)
      cz3b=cr3*(s3phiu-s3phid)+cs4*(s3psiu+s3psid)
c				Source psi*sign(z-z0)
      cz4=(s4phiu+s4phid)+cr1*(s4psiu-s4psid)
      cz4b=cr3*(s4phiu-s4phid)+cs4*(s4psiu+s4psid)


c	Fx et Fy composantes horizontales de la source
c	ur => u(1) 
c       ut => u(2)
c	uz => u(3)

      cu=cs2*cz1 + cz4
      cu2=cs3*s5
      cdu1=  k5(ir)*cu + k2(ir)*cu2
      cdu2=-(k2(ir)*cu + k5(ir)*cu2)
      cdu3=j1(ir)*(cs2*cz1b + cz4b)
      u(ir,1)=cdu1 + u(ir,1)
      u(ir,2)=cdu2 + u(ir,2)
      u(ir,3)=cdu3 + u(ir,3)
      du1=abs(u(ir,1))*uconv
      tconv=((abs(cdu1).le.du1).and.(tconv))
      du2=abs(u(ir,2))*uconv
      tconv=((abs(cdu2).le.du2).and.(tconv))
      du3=abs(u(ir,3))*uconv
      tconv=((abs(cdu3).le.du3).and.(tconv))

c	Fz composante verticale de la source
c	ur => u(4)
c	ut = 0
c	uz => u(5) 

      cdu4=-j1(ir)*cs4*(cz2  + cz3/cs5)
      cdu5= j0(ir)*kr* (cz2b + cz3b/cs5)
      u(ir,4)=cdu4 + u(ir,4)
      u(ir,5)=cdu5 + u(ir,5)
      du4=abs(u(ir,4))*uconv
      tconv=((abs(cdu4).le.du4).and.(tconv))
      du5=abs(u(ir,5))*uconv
      tconv=((abs(cdu5).le.du5).and.(tconv))

 30   continue

      return
      end
