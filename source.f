c******************************************************************************
c*                                                         AXITRA Version 3.0 *
c*			PROGRAMME SOURCE                                      *
c*									      *
c*   	Calcul de la reponse frequentielle due a une source particuliere      *
c*	                                                                      *
c*      Entree : Fonction de transfert Deplacement/Force ??
c*      Sortie : sismogramme en vitesse (m/s)                                 *
c*	                                                                      *
c******************************************************************************

      subroutine source( )

      parameter   (ncp=15,nrp=3000,ntp=2048,nlsrc=100)

      dimension   hc(ncp),vp(ncp),vs(ncp),rho(ncp),xr(nrp),
     &            yr(nrp),zr(nrp),a(6),qp(ncp),qs(ncp)
      complex*16  urs(ntp,nrp,6),uts(ntp,nrp,6),uzs(ntp,nrp,6),ai
      complex*16  ux(ntp),uy(ntp),uz(ntp),uu,omega,us(ntp),shx
      character   filein*11,fileout*10
      character   filename*300
      integer     iwk(ntp)
      real*8      lambda
      real*8      uur,uui
      namelist /input/nc,nfreq,tl,aw,nr,xl,ikmax,uconv,xs,ys,zs,nsrc

c     added by Hugo, to compute several sub-faults
      complex*16  unx(ntp), uny(ntp), unz(ntp)
      real        pr(nrp,3),xsr(1,3),d(nrp,3),teta(nrp)
      complex*16  un(ntp), ue(ntp)    !N-S and E-W components
      real        sca   !scale factor, not sure is necessary
      integer     iforc, nforc
      integer     reclent
      real        vector_out(ntp)

      pi=3.1415926535
      pi2=2.*pi
      ai=(0.,1.)

      filein='axi.data'
      fileout='axi.sis'


      open (10,form='formatted',file=filein)
      open (9,form='unformatted',file='axi.res')
      open (11,form='unformatted',file=fileout)
      rewind(11)

c++++++++++++
c                 Lecture de la fonction source
c++++++++++++

      write(6,*) 'fonction source ?'
      write(6,*) '1 : Ricker '
      write(6,*) '2 : step '
      write(6,*) '3 : fonction stockee sur fichier <axi.sou>'
      write(6,*) '4 : triangle '
      write(6,*) '5 : Gauss'
c      read(5,*) ics
c     fixed to Gaussian
      ics = 5
      write(6,*) 'Fixed to Gaussian shape ics: ', ics


      write(6,*) 'Sismogramme ?'
      write(6,*) '1: en deplacement (m)'
      write(6,*) '2: en vitesse (m/s)'
      write(6,*) '3: en acceleration (m/s/s)'
c      read(5,*) icc
c     fixed to deplacement
      icc = 1
      write(6,*) 'Fixed to displacement icc: ', icc

c	Le parametre t0 permet d'introduire un offset dans la
c	fonction source et le sismogramme. Il n'est pas utilise ici
 	    write(6,*) 'dt0'
c 	    read(5,*) dt0
c     delay of 1. sec
      dt0=1.
      write(6,*) 'delay dt0: ', dt0
      
c++++++++++++
c                 Lecture des fonctions de transfert
c++++++++++++

      read(10,input)

c      read(10,*)
      do 3 ic=1,nc
  3    read(10,*) hc(ic),vp(ic),vs(ic),rho(ic),qp(ic),qs(ic)
      do 4 ir=1,nr
      read(10,*) xr(ir),yr(ir),zr(ir)
!     THIS PART IS MODIFIED BY HUGO TO KNOW THE X Y Z components
      pr(ir,1:3) = [xr(ir), yr(ir), zr(ir)]
      xsr(1,1:3) = [xs, ys, zs]
      d(ir,:) = pr(ir,:) - xsr(1,:)
      teta(ir) = atan(d(ir,2)/d(ir,1))
       if ( (d(ir,1) .gt. 0) .and. (d(ir,2) .gt. 0 ) ) then
        !nothing to do, angle is in first quadrant
       elseif ( (d(ir,1) .gt. 0) .and. (d(ir,2) .lt. 0 ) ) then
        !second quadrant
        teta(ir) = pi2 + teta(ir)
       elseif ( (d(ir,1) .lt. 0) .and. (d(ir,2) .lt. 0 ) ) then
        !third quadrant
        teta(ir) = teta(ir) + pi
       elseif ( (d(ir,1) .lt. 0) .and. (d(ir,2) .gt. 0 ) ) then
        !fourht quadrant
        teta(ir) = teta(ir) + pi
       else
       endif
 4    print *, ir, pr(ir,:)!, d(ir,:), teta(ir)
!---------------------------------------------------------------

      close(10)
      close(11)

      do 5 if=1,nfreq
      read(9)((urs(if,ir,is),is=1,3),ir=1,nr)
      read(9)((uts(if,ir,is),is=1,3),ir=1,nr)
      read(9)((uzs(if,ir,is),is=1,3),ir=1,nr)
 1000 format(3(2e14.8,1x))
 5    continue
      close(9)

c     not usefull for computation I guess
c      do 44 ir=1,nr
c      filename = '';
c      write(filename,1066) 'f_xr_yr_zr_xt_yt_zt_xz_yz_zz_R',ir,'r.txt'
c      open (100000+ir+100,form='formatted',file=trim(filename))
c      filename = '';
c      write(filename,1066) 'f_xr_yr_zr_xt_yt_zt_xz_yz_zz_R',ir,'t.txt'
c      open (100000+ir+200,form='formatted',file=trim(filename))
c      filename = '';
c      write(filename,1066) 'f_xr_yr_zr_xt_yt_zt_xz_yz_zz_R',ir,'z.txt'
c      open (100000+ir+300,form='formatted',file=trim(filename))
c 1066 format(A,I0,A)
c      do 55 if=1,nfreq
c      write(100000+ir+100,'(i0)') if
c      write(100000+ir+200,'(i0)') if
c      write(100000+ir+300,'(i0)') if
c      write(100000+ir+100,1055) (real(urs(if,ir,is)),is=1,3)
c      write(100000+ir+200,1055) (real(uts(if,ir,is)),is=1,3)
c      write(100000+ir+300,1055) (real(uzs(if,ir,is)),is=1,3)
c      write(100000+ir+100,1055) (imag(urs(if,ir,is)),is=1,3)
c      write(100000+ir+200,1055) (imag(uts(if,ir,is)),is=1,3)
c      write(100000+ir+300,1055) (imag(uzs(if,ir,is)),is=1,3)
c 55   continue
c 1055 format(3(f23.18,2x))
c      close(100000+ir+100)
c      close(100000+ir+200)
c      close(100000+ir+300)
c 44   continue
c     end of useless writing stuff

c++++++++++++
c                 isc = indice de la couche contenant la source
c++++++++++++

      hh=0.
      isc=1
      zsc=zs
      do 9 ic=1,nc-1
      hh=hc(ic)
      if (zsc.gt.hh) then
       zsc=zsc-hh
       isc=ic+1
      else
       goto 91
      endif
  9   continue
  91  continue
      write(6,*) 'Source in layer= ',isc

c++++++++++++
c                 Parametres du sismogramme, nbre de points = 2**mm
c++++++++++++
      write(6,*) 'nfreq= ',nfreq
      xmm=log(real(nfreq))/log(2.)
      mm=int(xmm)+2
      mm=min(mm,11)
      nt=2**mm
c      write(6,*) 'xmm=',xmm
c      write(6,*) 'mm=',mm
      write(6,*) 'nt= ',nt
      
      dfreq=1./tl
      write(6,*) 'dfreq= ',dfreq
c     write(11) nr,nt,tl/nt
      aw=-pi*aw/tl
      rmu=vs(isc)*vs(isc)*rho(isc)
      rlambda=vp(isc)*vp(isc)*rho(isc)-2.*rmu
      
c++++++++++++
c                 Choix du type de source
c++++++++++++

c     call system('clear')
c      write (6,*) 'composante de la force ?'
c      write (6,*)
c      write(6,*) 'Fx ?'
c      read(5,*) a(1)
c      write(6,*) 'Fy ?'
c      read(5,*) a(2)
c      write(6,*) 'Fz ?'
c      read(5,*) a(3)
c      write(6,*)
      write(6,*) 'Force 1:X 2:Y 3:Z '


c++++++++++++
c                 Fonction source et choix de la frequence de la source
c++++++++++++
      if ((ics.eq.1).or.(ics.eq.2).or.(ics.eq.4)) t0=tl*3.5/nfreq
      if (ics.eq.3) open (12,form='formatted',file='axi.sou')
     	if (ics.eq.1) then
     	  write(6,*) "t0?"
	      	read(5,*) t0
     	endif
     	if (ics.eq.5) then
c        write(6,*) "tp"
c        read(5,*) t0
        t0 = 0.5
        write(6,*) 'Period Gaussian tp: ', tp
      endif
      us = (0.,0.)
      do if=1,nfreq
        freq=float(if-1)/tl
        omega=cmplx(pi2*freq,aw)
        if (icc.eq.1) deriv=1.
        if (icc.eq.2) deriv=ai*omega
        if (icc.eq.3) deriv=(ai*omega)*(ai*omega)

c               Source = Ricker en deplacement
     	if (ics.eq.1) then
        uu=omega*t0
        uu=uu*uu/pi2/pi2
        uu=exp(-uu)
        uu=omega*omega*uu*tl/nt
        us(if)=exp(-ai*omega*dt0)*uu*deriv
	     endif

c		Source = triangle en deplacement
	     if (ics.eq.4) then
	       t0=.5
	       uu=exp(ai*omega*t0/4.)
	       uu=(uu-1./uu)/2./ai
	       uu=uu/(omega*t0/2.)
	       us(if)=uu*uu**exp(-ai*omega*dt0)*deriv
	     endif

c               Source = step en deplacement
	     if (ics.eq.2) then
        shx=exp(omega*pi*t0/3.)
        shx=1./(shx-1./shx)
        uu=-ai*t0*pi*shx*tl/nt
        us(if)=exp(-ai*omega*dt0)*uu*deriv
	     endif

c               Source = fichier 'axi.sou'
c		sismogramme dans l'unite choisie dans le fichier
      if (ics.eq.3) then
        read(12,*) uur,uui
        us(if)=exp(-ai*omega*dt0)*cmplx(uur,uui)
	     endif
	     
c     gausiana
	     if (ics.eq.5) then 
	       tp=t0 * pi/4
        ts=dt0
        fmax = (nt-1) * dfreq
c       f = linspace(0,fmax,nf);
        omega_p = 2*pi / tp
        omega = 2*pi*freq
        b = omega / omega_p;
c     %   factorDeEscala = (tp/pi^.5); % Factor de escala teÛrico
        us(if) =  exp(-b**2) * exp(-ai * omega * ts) 
c       write(6,*) us(if)
	     endif
      enddo
      if (ics.eq.3) close(12)
      
c  graficar funciona de amplitud
      open(5429321,form='formatted',file='funcAmp.txt')
      do if=1,nfreq
        write(5429321,1054) real(us(if)), imag(us(if))
      enddo
      close(5429321)
 1054 format(2(f23.18,2x))
c++++++++++++
c                Calcul des sismogrammes
c++++++++++++

c++++++++++++
c                Opening output files
c++++++++++++
      reclent = nt*4
      open (20,form='unformatted',file='axi.x',
     &       access='DIRECT',recl=reclent)
      open (21,form='unformatted',file='axi.y',
     &       access='DIRECT',recl=reclent)
      open (22,form='unformatted',file='axi.z',
     &       access='DIRECT',recl=reclent)
      open (23,form='unformatted',file='axi.n',
     &       access='DIRECT',recl=reclent)
      open (24,form='unformatted',file='axi.e',
     &       access='DIRECT',recl=reclent)





c     loop over three forces
      iforcir = 1
      nforc = 3
      do 300 iforc = 1,nforc
         if ( iforc .eq. 1 ) then
          a(1:3) = [1., 0., 0.]
         elseif ( iforc .eq. 2 ) then
          a(1:3) = [0., 1., 0.]
         else
          a(1:3) = [0., 0., 1.]
         endif

      do 30 ir=1,nr
      ux = (0.,0.);      uy = (0.,0.);      uz = (0.,0.)
c      do if=1,nt
c      ux(if)=(0.,0.)
c      uz(if)=(0.,0.)
c      uy(if)=(0.,0.)
c      enddo
      
        do is=1,3
c         write(6,*) a(is)
          do if=1,nfreq
c          write(6,*) urs(if,ir,is), uts(if,ir,is), uzs(if,ir,is)
          ux(if)=a(is)*urs(if,ir,is)+ux(if)
          uy(if)=a(is)*uts(if,ir,is)+uy(if)
          uz(if)=a(is)*uzs(if,ir,is)+uz(if)
          enddo
        enddo
c               multiplication par la fonction source    
      do if=1,nfreq
c        write(6,*) ir,if,us(if)
c        write(6,*) ux(if),uy(if),uz(if)
        ux(if)= us(if)*ux(if)
        uy(if)= us(if)*uy(if)
        uz(if)= us(if)*uz(if)
c        write(6,*) ux(if)
      enddo
c		on complete le spectre pour les hautes frequences
c               avec inversion du signe de la partie imaginaire pour
c               la FFT inverse qui n'existe pas avec fft2cd
      do if=nt+2-nfreq,nt
        ux(if)=conjg(ux(nt+2-if))
	       uy(if)=conjg(uy(nt+2-if))
	       uz(if)=conjg(uz(nt+2-if))
      enddo
 
      call fft2cd(ux,mm,iwk)
      call fft2cd(uy,mm,iwk)
      call fft2cd(uz,mm,iwk)


!=====================================
C     THIS IS NEEDED TO SCALE SISMOGRAMS
      sca = (tl/(nt-1))/( real(nfreq))/sqrt(real(nt*2)) ! escala
!=====================================
 
      do it=1,nt
        ck=float(it-1)/nt
        cc=exp(-aw*tl*ck)/tl
        ux(it)=ux(it)*cc*sca
        uy(it)=uy(it)*cc*sca
        uz(it)=uz(it)*cc*sca
!------ ROTATION TO CARTESIAN COORDINATES
        un(it) = sin(teta(ir))*ux(it)-cos(teta(ir))*uy(it);
        ue(it) = cos(teta(ir))*ux(it)+sin(teta(ir))*uy(it);
!----------------------------------------------------------------
      enddo


c      do it=1,nt
c       vector_out(1:nt) = ((sngl(dreal(ux(1:nt)))))
c       write(20,rec=iforcir) vector_out
c       vector_out(1:nt) = ((sngl(dreal(ux(1:nt)))))
c       write(20,rec=iforcir) ((sngl(dreal(ux(it)))))
c       vector_out(1:nt) = ((sngl(dreal(ux(1:nt)))))
c       write(21,rec=iforcir) ((sngl(dreal(uy(it)))))
       vector_out(1:nt) = ((sngl(dreal(uz(1:nt)))))
       write(22,rec=iforcir) vector_out(1:nt)
       vector_out(1:nt) = ((sngl(dreal(un(1:nt)))))
       write(23,rec=iforcir) vector_out(1:nt)
       vector_out(1:nt) = ((sngl(dreal(ue(1:nt)))))
       write(24,rec=iforcir) vector_out(1:nt)
c      enddo
      iforcir = iforcir + 1
  30  continue
 300  continue

c      close(20)
c      close(21)
      close(22)
      close(23)
      close(24)


	return
	endsubroutine source
