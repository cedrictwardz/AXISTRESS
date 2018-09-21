c******************************************************************************
c*                                                         AXITRA Version 3.0 *
c*			PROGRAMME AXITRA                                      *
c*									      *
c*   	Calcul de sismogrammes synthetiques en milieu stratifie a symetrie    *
c*      cylindrique. Version avec force unidirectionelle uniquement.          *
c*	Propagation par la methode de la reflectivite, avec coordonnees       *
c*      cylindriques (r, theta, z)				              *
c*      Attenuation sur les ondes P et S                                      *
c*									      *
c*      auteur : Coutant O. (Vers 3.0, 1.03.91)                               *
c*	Bibliographie :                                                       *
c*                      Kennett GJRAS vol57, pp557R, 1979                     *
c*			Bouchon JGR vol71, n4, pp959, 1981                    *
c*									      *
c******************************************************************************

      program axitra

      include   "parameter.inc"
      include   "dimension1.inc"
 
      namelist/input/nc,nfreq,tl,aw,nr,xl,ikmax,uconv,xs,ys,zs,nsrc
c     data        ai,pi,pi2/(0.,1.),3.14159265359,6.28318530718/
      ai=(0.,1.)
      pi=3.14159265359
      pi2=6.28318530718

      open (in,form='formatted',file='axi.data')
      open (out,form='formatted',file='axi.head')
      open (12,form='unformatted',file='axi.res')
      rewind(out)
c++++++++++
c           LECTURE DES PARAMETRES D'ENTREE
c              
c               sismogramme : nfreq,tl,xl
c               recepteurs  : nr,xr(),yr(),zr()
c               source      : xs,ys,zs
c               modele      : nc,hc(),vp(),vs(),rho()
c                            
c               si hc(1)=0 on donne les profondeurs des interfaces, sinon
c               on donne les epaisseurs des couches
c++++++++++

c     read input parameters
      read(in,input)

c     read layer properties
      do 1 ic=1,nc
 1     read(in,*) hc(ic),vp(ic),vs(ic),rho(ic),qp(ic),qs(ic)

c     read coordinates of receivers
      do 2 ir=1,nr
 2     read(in,*) xr(ir),yr(ir),zr(ir)

c     write parameters read
      write(out,input)
      write(out,*) 'hc,vp,vs,rho,Qp,Qs'
      do 3 ic=1,nc
 3    write(out,1001) hc(ic),vp(ic),vs(ic),rho(ic),qp(ic),qs(ic)

c		Test sur les dimensions

      if ((nr.gt.nrp).or.(nc.gt.ncp)) then
      write(6,*) 'nombre de parametres superieur au dimensionnement'
      stop
      endif

c++++++++++
c           INITIALISATIONS
c++++++++++

      dfreq=1./tl
      aw=-pi*aw/tl
      freq=-dfreq
      pil=pi2/xl
      ikmin=100

      call initdata

      write(out,*) 'xr,yr,zr'
      do 4 ir=1,nr
 4    write(out,*) xr(ir),yr(ir),zr(ir)
      
c               ***************************
c               ***************************
c               **  BOUCLE EN FREQUENCE  **
c               ***************************
c               ***************************
      
      do 10 if=1,nfreq
	write(0,*) 'freq ',if,'/',nfreq
      freq=freq+dfreq
      rw=pi2*freq
      omega=cmplx(rw,aw)
      omega2=omega*omega
      a1=.5/omega2/xl
      zom=sqrt(rw*rw+aw*aw)
      if (if.eq.1) then
      phi=-pi/2
      else
      phi=atan(aw/rw)
      endif
      xlnf=(ai*phi+alog(real(zom)))/pi


c            ******************************************
c            ******************************************
c            **  RESOLUTION PAR BOUCLE EXTERNE EN Kr **
c            ******************************************
c            ******************************************

      do 20 ik=0,ikmax
      
      kr=(ik+.258)*pil
      kr2=kr*kr
      tconv=.true.


c+++++++++++++
c              Calcul de nombreux coefficients et des fonctions de Bessel
c+++++++++++++

      call reflect0

c+++++++++++++
c              Calcul des coefficients de reflexion/transmission
c	       Matrice de Reflection/Transmission et Dephasage
c+++++++++++++

      call reflect1
      
c+++++++++++++
c              Calcul des matrices de reflectivite : mt(),mb(),nt(),nb()
c              (rapport des differents potentiels montant/descendant
c                        en haut et en bas de chaque couche)
c+++++++++++++

      call reflect2

c+++++++++++++
c	       Calcul des matrices de passage des vecteurs potentiel 
c		source, aux vecteurs potentiel PHI, PSI et KHI au sommet
c		de chaque couche
c+++++++++++++

      call reflect3

c+++++++++++++
c	       Calcul des potentiels et des deplacement dus aux sources du
c		tenseur, en chaque recepteur (termes en kr, r, z)
c+++++++++++++
    
      call reflect4v3
     
      if ((tconv).and.(ik.gt.ikmin)) goto 21
 20   continue
 21   continue

c+++++++++++++
c               Calcul des deplacements aux recepteurs 
c+++++++++++++

      call reflect5v3

c+++++++++++++
c		Sortie des resultats
c+++++++++++++
      write(out,*) 'freq =',freq,'iter =',ik
c     write(out,*) 'ux'
      write(12) ((ur(ir,is),is=1,3),ir=1,nr)
c     write(out,*) 'uy'
      write(12) ((ut(ir,is),is=1,3),ir=1,nr)
c     write(out,*) 'uz'
      write(12) ((uz(ir,is),is=1,3),ir=1,nr)
c     if (ik.ge.ikmax) then
c     write(6,*) 'Depassement du nombre d iteration maximum'
c     stop
c     endif
      
 1000 format(3(2e14.8,1x))
 1001 format(8f9.3)

 10   continue

      close(in)
      close(out)
      close(12)

      call source()
      call stress_tensor(nsrc)

      stop
      end
