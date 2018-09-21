       subroutine stress_tensor (nsrc)

       IMPLICIT NONE
       integer*4 i,j,k,ic,ncp,iunitx,iunity,iunitz,iunit
       integer nc
       parameter (ncp=16)

       real*4, dimension(:,:), allocatable :: tst
       real*4, dimension(:,:), allocatable :: sis
       real*4, dimension(:), allocatable :: dxx,dxy,dxz
       real*4, dimension(:), allocatable :: dyx,dyy,dyz
       real*4, dimension(:), allocatable :: dzx,dzy,dzz,ldil
       real*4 x_step, hc(ncp),vp(ncp),vs(ncp),rho(ncp),qp(ncp),qs(ncp)
       real*4 mu(ncp),lam(ncp),z,mu2
       integer*4 tsam, ttsam, nr
       character*9 filein
!     added by Hugo, to compute several sub-faults
      CHARACTER*4 source

      integer iforc, nforc, ifai, nfai, isrc, nsrc
      integer reclent, irec, recwlen
      real, dimension(:), allocatable :: vector_out


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

       nfai = nr / 6

       iunitx=20
       iunity=21
       iunitz=22
       iunit=23

       write(source,'(I4.4)') nsrc
       recwlen = tsam*6*4
       allocate(vector_out(tsam*6))
       open(iunitx,file='st_'//source//'x.bin',status='unknown',
     &      form='unformatted',access='direct',recl=recwlen)
       open(iunity,file='st_'//source//'y.bin',status='unknown',
     &      form='unformatted',access='direct',recl=recwlen)
       open(iunitz,file='st_'//source//'z.bin',status='unknown',
     &      form='unformatted',access='direct',recl=recwlen)
       close(iunitx)
       close(iunity)
       close(iunitz)



       nforc = 3

       ttsam=tsam*nr*nforc

       reclent = ttsam*4
       open(iunitx,file='axi.e',action='read',status='old',
     &      access='direct',form='unformatted',recl=reclent)
       open(iunity,file='axi.n',action='read',status='old',
     &      access='direct',form='unformatted',recl=reclent)
       open(iunitz,file='axi.z',action='read',status='old',
     &      access='direct',form='unformatted',recl=reclent)

       allocate(sis(ttsam,3))
       allocate(dxx(tsam),dxy(tsam),dxz(tsam))
       allocate(dyx(tsam),dyy(tsam),dyz(tsam))
       allocate(dzx(tsam),dzy(tsam),dzz(tsam))
       allocate(ldil(tsam))

       read(iunitx,rec=1) sis(:,1)
       read(iunity,rec=1) sis(:,2)
       read(iunitz,rec=1) sis(:,3)

       close(iunitx)
       close(iunity)
       close(iunitz)

       !-------Derivatives with respect to X-------------!

       allocate(tst(6,tsam))


       do ifai = 1,nfai
       do iforc = 1,nforc

       i = 1 + (iforc-1)*(tsam*nr) + (ifai-1)*(tsam*6)
       j = tsam * 1 + (iforc-1)*(tsam*nr) + (ifai-1)*(tsam*6)
       k = tsam * 2 + (iforc-1)*(tsam*nr) + (ifai-1)*(tsam*6)


       print *, i, j, k, iforc

       !      dxx 2nd order FD
       dxx = (sis(j+1:k,1) - sis(i:j,1))/(2.*x_step)
       !      dyx 2nd order FD
       dyx = (sis(j+1:k,2) - sis(i:j,2))/(2.*x_step)
       !      dzx 2nd order FD
       dzx = (sis(j+1:k,3) - sis(i:j,3))/(2.*x_step)

       !-------Derivatives with respect to Y-------------!

       i = tsam * 2 + 1 + (iforc-1)*(tsam*nr) + (ifai-1)*(tsam*6)
       j = tsam * 3 + (iforc-1)*(tsam*nr) + (ifai-1)*(tsam*6)
       k = tsam * 4 + (iforc-1)*(tsam*nr) + (ifai-1)*(tsam*6)
       print *, i, j, k, iforc

       !                dxy
       dxy = (sis(j+1:k,1) - sis(i:j,1))/(2.*x_step)
       !                dyy
       dyy = (sis(j+1:k,2) - sis(i:j,2))/(2.*x_step)
       !                dzy
       dzy = (sis(j+1:k,3) - sis(i:j,3))/(2.*x_step)

       !-------Derivatives with respect to Z-------------!

       i = tsam * 4 + 1 + (iforc-1)*(tsam*nr) + (ifai-1)*(tsam*6)
       j = tsam * 5 + (iforc-1)*(tsam*nr) + (ifai-1)*(tsam*6)
       k = tsam * 6 + (iforc-1)*(tsam*nr) + (ifai-1)*(tsam*6)
       print *, i, j, k, iforc

       !                dxz
       dxz = (sis(j+1:k,1) - sis(i:j,1))/(2.*x_step)
       !                dyz
       dyz = (sis(j+1:k,2) - sis(i:j,2))/(2.*x_step)
       !                dzz
       dzz = (sis(j+1:k,3) - sis(i:j,3))/(2.*x_step)


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
       tst(1,:) = ldil + ( mu2 * dxx )
       !            Syy
       tst(2,:) = ldil + ( mu2 * dyy )
       !            Szz
       tst(3,:) = ldil + ( mu2 * dzz )
       !            Sxy
       tst(4,:) = mu2 * 0.5 * ( dxy + dyx )
       !            Sxz
       tst(5,:) = mu2 * 0.5 * ( dxz + dzx )
       !            Syz
       tst(6,:) = mu2 * 0.5 * ( dyz + dzy )

      vector_out=[tst(1,:),tst(2,:),tst(3,:),tst(4,:),tst(5,:),tst(6,:)]

      irec = ifai 
      print *, irec, 'irec'
       if ( iforc .eq. 1 ) then
       open(iunit,file='st_'//source//'x.bin',status='old',
     &      form='unformatted',access='direct',recl=recwlen)
        write(iunit,rec=irec) vector_out
       elseif ( iforc .eq. 2 ) then
       open(iunit,file='st_'//source//'y.bin',status='old',
     &      form='unformatted',access='direct',recl=recwlen)
        write(iunit,rec=irec) vector_out
       elseif ( iforc .eq. 3 ) then
       open(iunit,file='st_'//source//'z.bin',status='old',
     &      form='unformatted',access='direct',recl=recwlen)
        write(iunit,rec=irec) vector_out
       endif
       close(iunit)

       enddo !nforc
       enddo !nfai
       !print *, 'Finished'
       
       deallocate(sis,dxx,dyx,dzx,dxy,dyy,dzy,dxz,dyz,dzz)
       deallocate(tst,ldil)
       deallocate(vector_out)
       end subroutine stress_tensor
