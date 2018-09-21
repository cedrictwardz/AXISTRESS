
      program cut

      implicit none
      integer iunit, i, j, k, ii, jj, irec
      integer tsam, ncomp, nfai, nsta
      real, dimension(:), allocatable :: st1, st2, vector
      integer reclenr, reclenw
      character*4 sta

      nsta = 2

      tsam = 2048
      ncomp = 6
      nfai = 3
      reclenr = tsam*ncomp*nfai*4
      reclenw = tsam*ncomp*4

      allocate(st1(reclenr/4))
      allocate(st2(reclenr/4))
      allocate(vector(tsam*ncomp))

      do i = 1,nsta
       write(sta,'(I4.4)') i
       iunit=22
       open(iunit,file='st_'//sta//'x.bin',form='unformatted',&
  &         access='direct',recl=reclenr)
       read(iunit,rec=1) st1
       close(iunit)
       iunit=22
       open(iunit,file='st_x.bin',form='unformatted',&
  &         access='direct',recl=reclenw)
       do j = 1,nfai
        irec = i+(j-1)*nsta
        ii = tsam*ncomp*(j-1)+1
        jj = tsam*ncomp*j
        print *, ii,jj, i, j, irec
        vector = st1(ii:jj)
        write(iunit,rec=irec) vector
       enddo
      enddo
      close(iunit)


      do i = 1,nsta
       write(sta,'(I4.4)') i
       iunit=22
       open(iunit,file='st_'//sta//'y.bin',form='unformatted',&
  &         access='direct',recl=reclenr)
       read(iunit,rec=1) st1
       close(iunit)
       iunit=22
       open(iunit,file='st_y.bin',form='unformatted',&
  &         access='direct',recl=reclenw)
       do j = 1,nfai
        irec = i+(j-1)*nsta
        ii = tsam*ncomp*(j-1)+1
        jj = tsam*ncomp*j
        print *, ii,jj, i, j, irec
        vector = st1(ii:jj)
        write(iunit,rec=irec) vector
       enddo
      enddo
      close(iunit)


      do i = 1,nsta
       write(sta,'(I4.4)') i
       iunit=22
       open(iunit,file='st_'//sta//'z.bin',form='unformatted',&
  &         access='direct',recl=reclenr)
       read(iunit,rec=1) st1
       close(iunit)
       iunit=22
       open(iunit,file='st_z.bin',form='unformatted',&
  &         access='direct',recl=reclenw)
       do j = 1,nfai
        irec = i+(j-1)*nsta
        ii = tsam*ncomp*(j-1)+1
        jj = tsam*ncomp*j
        print *, ii,jj, i, j, irec
        vector = st1(ii:jj)
        write(iunit,rec=irec) vector
       enddo
      enddo
      close(iunit)



      deallocate(st1,st2,vector)
      end program cut
