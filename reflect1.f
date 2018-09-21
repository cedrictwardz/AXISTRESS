c******************************************************************************
c                                                          AXITRA Version 3.0 *
c                         SUBROUTINE REFLECT1                                 *
c                                                                             *
c               Calcul des coefficients de reflexion/transmission             *
c		   Matrice de Reflexion/Transmission et Dephasage             *
c       (Les coefficients de reflexion/transmission utilisent les memes       *
c            termes intermediaires que Aki-Richards p149, MAIS :              *
c            Aki utilise la convention inverse pour la TF (exp(-iwt)),        *
c        et travaille avec le parametre de rai et les angles d'incidences)    *
c                                                                             *
c      Le potentiel PSI utilise pour l'onde SV est defini par :               *
c                  u = rot ( rot (PSI) )                                      *
c      i.e. un terme de derivation supplementaire par rapport a la convention *
c      habituelle : u= rot (PSI)                                              *
c                                                                             *
c   On deduit les coefficients de REF/TRANS de ceux definis par la convention *
c      classique en divisant le potentiel PSI par 1./ai/kr = coef             *
c                                                                             *
c       Ordre de stockage :                                                   *
c                              pp=(1,1)   sp=(1,2)                            *
c                              ps=(2,1)   ss=(2,2)                            *
c******************************************************************************


      subroutine reflect1

      include "parameter.inc" 
      include "dimension1.inc"
      include "dimension2.inc"
c Coefficient pour la convention sur PSI (coef) et sur la TF (aki=-1.)
      coef=1./ai
      aki=-1.
      
c               CONDITIONS AUX LIMITES a la profondeur z=0, coefficients
c               de reflexion/transmission
c      2 possibilites : 1) surface libre (mettre en commentaire de B1 a B2)
c                       2) 1/2 espace sup. infini (commentaires de A1 a A2)

cA1                    SURFACE LIBRE 
      cf1=(ckb2(1)-2.*kr2)
      cf2=cf1*cf1
      cf3=4.*cnu(1)*kr2*cgam(1)
      cdd=cf2+cf3

      ru(1,1,1)=(-cf2+cf3)/cdd
      ru(1,2,1)=4.*cnu(1)*cf1/cdd*coef*aki
      ru(1,2,2)=(cf2-cf3)/cdd*aki
      ru(1,1,2)=4.*kr2*cgam(1)*cf1/cdd/coef
      tu(1,1,1)=0.
      tu(1,1,2)=0.
      tu(1,2,1)=0.
      tu(1,2,2)=0.
      rush(1)=1.
      tush(1)=0.
cA2

cB1                   1/2 ESPACE SUP. INFINI  
c     ru(1,1,1)=0.
c     ru(1,2,1)=0.
c     ru(1,2,2)=0.
c     ru(1,1,2)=0.
c     tu(1,1,1)=1.
c     tu(1,1,2)=0.
c     tu(1,2,1)=0.
c     tu(1,2,2)=1.
c     rush(1)=0.
c     tush(1)=1.
cB2

c               Coefficients aux interfaces entre couches

      do 24 ic=2,nc

      cb1=kr2/ckb2(ic-1)
      cb2=kr2/ckb2(ic)
      ca1d=rho(ic-1)*(1.-2.*cb1)
      ca2d=rho(ic)*(1.-2.*cb2)
      ca=ca2d-ca1d
      cb=ca2d+2.*rho(ic-1)*cb1
      cc=ca1d+2.*rho(ic)*cb2
      cd=2.*(rho(ic)/ckb2(ic)-rho(ic-1)/ckb2(ic-1))
      ce=cb*cnu(ic-1)+cc*cnu(ic)
      cf=cb*cgam(ic-1)+cc*cgam(ic)
      cg=ca-cd*cnu(ic-1)*cgam(ic)
      ch=ca-cd*cnu(ic)*cgam(ic-1)
      cdd=ce*cf+cg*ch*kr2

      rd(ic,1,1)= (cf*(cb*cnu(ic-1)-cc*cnu(ic))-
     &            ch*kr2*(ca+cd*cnu(ic-1)*cgam(ic)))/cdd
      rd(ic,1,2)=-2.*kr2*cgam(ic-1)*
     &            (ca*cb+cc*cd*cnu(ic)*cgam(ic))/cdd/coef*aki
      rd(ic,2,2)=-(ce*(cb*cgam(ic-1)-cc*cgam(ic))-
     &            cg*kr2*(ca+cd*cnu(ic)*cgam(ic-1)))/cdd*aki
      rd(ic,2,1)=-2.*cnu(ic-1)*
     &            (ca*cb+cc*cd*cnu(ic)*cgam(ic))/cdd*coef
      td(ic,1,1)= 2.*rho(ic-1)*cnu(ic-1)*cf/cdd
      td(ic,1,2)=-2.*rho(ic-1)*cgam(ic-1)*cg*kr2/cdd/coef*aki
      td(ic,2,2)= 2.*rho(ic-1)*cgam(ic-1)*ce/cdd
      td(ic,2,1)= 2.*rho(ic-1)*cnu(ic-1)*ch/cdd*coef*aki

      ru(ic,1,1)=-(cf*(cb*cnu(ic-1)-cc*cnu(ic))+
     &            cg*kr2*(ca+cd*cnu(ic)*cgam(ic-1)))/cdd
      ru(ic,1,2)= 2.*kr2*cgam(ic)*
     &            (ca*cc+cb*cd*cnu(ic-1)*cgam(ic-1))/cdd/coef
      ru(ic,2,2)= (ce*(cb*cgam(ic-1)-cc*cgam(ic))+
     &            ch*kr2*(ca+cd*cnu(ic-1)*cgam(ic)))/cdd*aki
      ru(ic,2,1)= 2.*cnu(ic)*
     &            (ca*cc+cb*cd*cnu(ic-1)*cgam(ic-1))/cdd*coef*aki
      tu(ic,1,1)= 2.*rho(ic)*cnu(ic)*cf/cdd
      tu(ic,1,2)= 2.*rho(ic)*cgam(ic)*ch*kr2/cdd/coef
      tu(ic,2,2)= 2.*rho(ic)*cgam(ic)*ce/cdd
      tu(ic,2,1)=-2.*rho(ic)*cnu(ic)*cg/cdd*coef

c   Modification pour calculateur a faible dynamique [1.e-300; 1.e+300]

      cdeph=exp(-ai*cnu(ic-1)*hc(ic-1))
      rdeph=dreal(cdeph)
      adeph=dimag(cdeph)
      if (dabs(rdeph).lt.1.e-150) rdepth=0.
      me1(ic-1)=cmplx(rdeph,adeph)
      
      cdeph=exp(-ai*cgam(ic-1)*hc(ic-1))
      rdeph=dreal(cdeph)
      adeph=dimag(cdeph)
      if (dabs(rdeph).lt.1.e-150) rdepth=0.
      me2(ic-1)=cmplx(rdeph,adeph)

c   Version normal (ex Cray)
 
c     me1(ic-1)=exp(-ai*cnu(ic-1)*hc(ic-1))
c     me2(ic-1)=exp(-ai*cgam(ic-1)*hc(ic-1))

      cs1=rho(ic-1)/ckb2(ic-1)*cgam(ic-1)
      cs2=rho(ic)/ckb2(ic)*cgam(ic)
      cdelt=cs1+cs2

      rush(ic)=(cs2-cs1)/cdelt
      rdsh(ic)=-rush(ic)
      tush(ic)=2.*cs2/cdelt
      tdsh(ic)=2.*cs1/cdelt

 24   continue

      return
      end
