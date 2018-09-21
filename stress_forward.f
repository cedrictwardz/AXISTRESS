       subroutine stress_tensor()

       IMPLICIT NONE
       integer*4 i,j,k,ic,ncp,iunitx,iunity,iunitz,iunit
       integer nc
       parameter (ncp=16)

       real*4, dimension(:,:,:), allocatable :: tst
       real*4, dimension(:,:), allocatable :: sis
       real*4, dimension(:), allocatable :: dxx,dxy,dxz
       real*4, dimension(:), allocatable :: dyx,dyy,dyz
       real*4, dimension(:), allocatable :: dzx,dzy,dzz,ldil
       real*4 x_step, hc(ncp),vp(ncp),vs(ncp),rho(ncp),qp(ncp),qs(ncp)
       real*4 mu(ncp),lam(ncp),z,mu2
       integer*4 tsam, ttsam, nr
       character*9 filein
!     added by Hugo, to compute several sub-faults
      integer nfai, ifixx, jfixx
      CHARACTER*4 fai, ifor

      integer iforc, nforc

       filein='stress.in'
       open(15,file=filein,action='read',status='old')
       read(15,*) tsam
       read(15,*) nr
       read(15,*) x_step
       x_step = x_step*1000.
       read(15,*) nc
       do 3 ic=1,nc
       read(15,*) hc(ic),vp(ic),vs(ic),rho(ic),qp(ic),qs(ic)
       rho(ic) = rho(ic)/1000.
       vp(ic) = vp(ic)*1000.
       vs(ic) = vs(ic)*1000.
       mu(ic) = (vs(ic)**2)*rho(ic)
       lam(ic) = (vp(ic)**2)*rho(ic) - (2. * mu(ic))
 3     print *, mu(ic), lam(ic)
       read(15,*) z
       close(15)

       iunitx=20
       iunity=21
       iunitz=22
       iunit=23

       nforc = 3

       ttsam=tsam*nr*nforc

       open(iunitx,file='axi.e',action='read',status='old')
       open(iunity,file='axi.n',action='read',status='old')
       open(iunitz,file='axi.z',action='read',status='old')

       allocate(sis(ttsam,3))
       allocate(dxx(tsam),dxy(tsam),dxz(tsam))
       allocate(dyx(tsam),dyy(tsam),dyz(tsam))
       allocate(dzx(tsam),dzy(tsam),dzz(tsam))
       allocate(ldil(tsam))

       read(iunitx,*) sis(:,1)
       read(iunity,*) sis(:,2)
       read(iunitz,*) sis(:,3)

       close(iunitx)
       close(iunity)
       close(iunitz)

       !-------Derivatives with respect to X-------------!

       allocate(tst(1,6,tsam))


       do iforc = 1,nforc
        write(ifor,'(I4.4)') iforc


       i = 1 + (iforc-1)*(tsam*nr)
       j = tsam + (iforc-1)*(tsam*nr)
       k = tsam * 2 + (iforc-1)*(tsam*nr)
       ifixx = j+1
       jfixx = k

       print *, i, j, k, iforc

       !      dxx 2nd order FD
       dxx = (sis(i:j,1) - sis(ifixx:jfixx,1))/(x_step)
       !      dyx 2nd order FD
       dyx = (sis(i:j,2) - sis(ifixx:jfixx,2))/(x_step)
       !      dzx 2nd order FD
       dzx = (sis(i:j,3) - sis(ifixx:jfixx,3))/(x_step)

       !-------Derivatives with respect to Y-------------!

       i = tsam * 2 + 1 + (iforc-1)*(tsam*nr)
       j = tsam * 3 + (iforc-1)*(tsam*nr)

       !                dxy
       dxy = (sis(i:j,1) - sis(ifixx:jfixx,1))/(x_step)
       !                dyy
       dyy = (sis(i:j,2) - sis(ifixx:jfixx,2))/(x_step)
       !                dzy
       dzy = (sis(i:j,3) - sis(ifixx:jfixx,3))/(x_step)

       !-------Derivatives with respect to Z-------------!

       i = tsam * 3 + 1 + (iforc-1)*(tsam*nr)
       j = tsam * 4 + (iforc-1)*(tsam*nr)

       !                dxz
       dxz = (sis(i:j,1) - sis(ifixx:jfixx,1))/(x_step)
       !                dyz
       dyz = (sis(i:j,2) - sis(ifixx:jfixx,2))/(x_step)
       !                dzz
       dzz = (sis(i:j,3) - sis(ifixx:jfixx,3))/(x_step)


       ! Estimation of stress tensor 

       ! Lame parameters
       do ic=2,nc
         if ( (z .gt. hc(ic-1)) .and. (z .lt. hc(ic)) ) then
          mu2 = 2. * mu(ic-1)
          ldil = lam(ic-1) * (dxx + dyy + dzz)
         endif
       enddo
       !print *, mu2, lam(ic-1)

       !            Sxx
       tst(1,1,:) = ldil + ( mu2 * dxx )
       !            Sxy
       tst(1,2,:) = mu2 * 0.5 * ( dxy + dyx )
       !            Sxz
       tst(1,3,:) = mu2 * 0.5 * ( dxz + dzx )
       !            Syy
       tst(1,4,:) = ldil + ( mu2 * dyy )
       !            Syz
       tst(1,5,:) = mu2 * 0.5 * ( dyz + dzy )
       !            Szz
       tst(1,6,:) = ldil + ( mu2 * dzz )


       open(iunit,file='st_'//ifor//'.ascii',status='unknown')
       do i=1,tsam 
         write(iunit,*)  tst(1,1:6,i)
       enddo
       close(iunit)

       enddo
       !print *, 'Finished'
       
       deallocate(sis,dxx,dyx,dzx,dxy,dyy,dzy,dxz,dyz,dzz)
       deallocate(tst,ldil)
       end subroutine stress_tensor
