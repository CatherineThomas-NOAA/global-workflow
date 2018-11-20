!23456
      module fv3wam_grids_module
      IMPLICIT none
!
      integer, parameter    :: nzt62 = 82
      integer, parameter    :: nyt62 = 94
      integer, parameter    :: nxt62 = 192 
!
      integer, parameter    :: nzfv3 = 127
      integer, parameter    :: nxc096 = 96
      integer, parameter    :: nyc096 = 96
!
      integer, parameter    :: nxc192 = 192
      integer, parameter    :: nyc192 = 192
!
      integer, parameter    :: nxc384 = 384
      integer, parameter    :: nyc384 = 384

      integer, parameter    :: nxc768 = 768
      integer, parameter    :: nyc768 = 768
!
      integer               :: nlat, nlon, nlev, ntime
      integer               :: nlats, nlons, nlevs, nilevs
!
      real, allocatable     :: newhyam(:), newhybm(:)
      real, allocatable     :: newhyai(:), newhybi(:)
      real, allocatable     :: newlev(:), newlevi(:)
      real, allocatable     :: newlat(:), newlon(:)
      contains
!
      subroutine READ_FV3WAMGRIDS(CRES_FV3)
      use netcdf
!
      character(len=4) :: CRES_FV3              ! C096, C192 , C384, C786 
!
      integer :: ncid, status, RecordDimID
      integer :: iLevDimID, LevDimID, LatDimID, LonDimID
      integer :: latid, lonid, levid, ilevid
      character(len=132) filename/'wamt62_grids.nc'/
!
      if (CRES_FV3.eq.'C096') then
          filename='C096_grids.nc'
         filename='wamt62_grids.nc'
      endif
!
      if (CRES_FV3.eq.'C192') then
         filename='C192_grids.nc'
      endif
!
      if (CRES_FV3.eq.'C384') then
         filename='C384_grids.nc'
      endif
!
      if (CRES_FV3.eq.'C768') then
         filename='C768_grids.nc'
      endif 
      filename='wamt62_grids.nc'
      print*, filename  
!
        status = nf90_open(trim(FileName), nf90_nowrite, ncid)

        status = nf90_inquire(ncid, unlimitedDimId = RecordDimID)
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_inq_dimid(ncid, "ilev", iLevDimID)
        if (status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_dimid(ncid, "lev", LevDimID)
        if (status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_dimid(ncid, "lat", LatDimID)
        if (status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_dimid(ncid, "lon", LonDimID)
        if (status /= nf90_noerr) call handle_err(status)


       ! How many values of "ilev" are there?

       status = nf90_inquire_dimension(ncid, iLevDimID, len = niLevs)
            if (status /= nf90_noerr) call handle_err(status)

         status = nf90_inquire_dimension(ncid, LevDimID, len = nLevs)
            if (status /= nf90_noerr) call handle_err(status)

         status = nf90_inquire_dimension(ncid, LatDimID, len = nLats)
             if (status /= nf90_noerr) call handle_err(status)

         status = nf90_inquire_dimension(ncid, LonDimID, len = nLons)
             if (status /= nf90_noerr) call handle_err(status)

         if (nlevs .ne. nzt62 .or. nlats .ne. nyt62 .or. nlons .ne. nxt62) then
             print *, nlevs,  nzt62
             print *, nlats,  nyt62
             print *, nlons,  nxt62 
          stop ' dimensions in Filename '
         endif

         allocate(newlev(nLevs),newhyam(nLevs),newhybm(nLevs))
         allocate(newlevI(niLevs),newhyaI(niLevs),newhybI(niLevs))
         allocate(newlat(nLats),newlon(nLons))

         status = nf90_inq_varid(ncid, 'lev', levid)
         if (status /= nf90_noerr) call handle_err(status)
         status = nf90_get_var(ncid, levid, newlev)
         if (status /= nf90_noerr) call handle_err(status)

         status = nf90_inq_varid(ncid, 'ilev', ilevid)
         if (status /= nf90_noerr) call handle_err(status)
         status = nf90_get_var(ncid, ilevid, newlevi)
         if (status /= nf90_noerr) call handle_err(status)

         status = nf90_inq_varid(ncid, 'hyam', levid)
         status = nf90_get_var(ncid, levid, newhyam)
         status = nf90_inq_varid(ncid, 'hybm', levid)
         status = nf90_get_var(ncid, levid, newhybm)

         status = nf90_inq_varid(ncid, 'hyai', levid)
         status = nf90_get_var(ncid, levid, newhyai)
         status = nf90_inq_varid(ncid, 'hybi', levid)
         status = nf90_get_var(ncid, levid, newhybi)

          status = nf90_inq_varid(ncid, 'lon', lonid)
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_get_var(ncid, lonid, newlon)
          status = nf90_inq_varid(ncid, 'lat', latid)
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_get_var(ncid, latid, newlat)
!
          print *,  ' lats' , maxval(newlat), minval(newlat)
          print *,  ' lons' , maxval(newlon), minval(newlon)
          print *,  ' lev' ,   newlev(1), newlev(nlevs)
          print *,  ' ilev' ,   newlevI(1), newlevI(nlevs+1)
!

         nlat = nlats
         nlon = nlons
         nlev = nlevs
          print *, nlat, nlon, nlev, ntime, ' nlat, nlon, nlev, FV3WAM '
!
      end subroutine READ_FV3WAMGRIDS
!
      end module fv3wam_grids_module
!
      SUBROUTINE MAP_ZWAMFV3(nT, nx,ny, nz, A4in,P4in, nlev, A4, P4)
      implicit none
      integer :: nx, ny, nz, nt
      integer :: nlev
      real    :: A4in(nx, ny, nz, nt)
      real    :: A4(nx, ny, nlev, nt)
      real    :: P4in(nx, ny, nz, nt)
      real    :: P4(nx, ny, nlev, nt)
!
      integer :: i,j, k, it
      real    :: A(nz),  B(nlev)
      real    :: zA(nz), zB(nlev)  
      do it=1, nt
       do j=1, ny
       do i=1, nx
         A(1:nz) = A4in(i,j,1:nz,it)
         zA(1:nz) = P4in(i,j,1:nz,it)
         zb(1:nlev) = P4(i,j,1:nlev,it)
         call zint_1d_g5(nz, Za, A, nlev, Zb, B)
         A4( i,j,1:nlev,it) = B(1:nlev)
       enddo
       enddo
      enddo            
      end SUBROUTINE MAP_ZWAMFV3
!
      subroutine zint_1d_g5(na, Za, A, nb, Zb, B)
      implicit none 
      integer :: na, nb
      real, dimension(na) :: Za, A
      real, dimension(nb) :: Zb, B
      real :: w1, w2, X, plim1, plim2
      integer :: i,j, k, ka, ks
         plim1 = maxval(Za)
         plim2 = minval(Za)
          ks =1
        do k=nb, 1, -1
           X = Zb(k)
           if (x <= plim2 ) B(k) = A(1)  ! 1   top lid
           if (x >= plim1 ) B(k) = A(na) ! na - bottom
           if (X > plim2 .and. X < plim1) then
             do ka=ks, na-1
               if (X >= Za(ka) .and. X <= Za(ka+1) ) exit
             enddo
               ks = ka
               B(k) = A(ks)+(A(ks)-A(ks+1))/(Za(ks)-Za(ks+1))*(X-Za(ks))    
           endif
        enddo
      end subroutine zint_1d_g5 
!
      subroutine read_vgrid_geos5(filegrid, nz, hyam, hybm, hyai, hybi)
      implicit NONE
      character(len=*) :: filegrid
      integer :: nz, nzf
      real, dimension(nz)   :: hyam, hybm, zkm
      real, dimension(nz+1) :: hyai, hybi
      integer :: k
      filegrid='/scratch3/NCEPDEV/swpc/save/Valery.Yudin/SVN_VAY/FV3_TL80/geos5_72L.ff'
      open(unit=11, file=trim(filegrid), status='old')
      read(11,*) nzf
      prinT *, nzf, ' dim-from-file'
       
      read(11,*) hyam
      read(11,*) hybm
      read(11,*) zkm
      do k=1, nz
       write(6, 111) k, hyam(k), hybm(k), 1000.*(hyam(k)+hybm(k)), -7*alog(hyam(k)+hybm(k))
      enddo
111   format(I3, 4(2x, F15.9))
      close(11)
      end subroutine read_vgrid_geos5
!
     subroutine read_vgrid_gfs(filegrid, nz, hyam, hybm, hyai, hybi)
      character(len=*) :: filegrid
      integer :: nz, nzf, nzp
      real, dimension(nz)   :: hyam, hybm
      real, dimension(nz+1) :: hyai, hybi
      real, dimension(nz+1) ::     zkm, press
      real, dimension(nz+1,2) :: vc2 
!
      integer :: k
       filegrid='/scratch3/NCEPDEV/swpc/save/Valery.Yudin/SVN_VAY/FV3_TL80/fv3gfs_128L.ff'
      open(unit=11, file=trim(filegrid), status='old')
      read(11,*) nzf
      read(11,*) nzp
      prinT *, nzf, ' dim-from-file'
       
      read(11,*) hybi   ! 1 to 0 *PS
      read(11,*) hyai   ! logp
      read(11,*) press
      read(11,*) zkm
      read(11,*) vc2
      close(11)
      hyai = 1.e-5*hyai
      do k=1, nz
        hyam(k) = .5*(hyai(k)+hyai(k+1))
        hybm(k) = .5*(hybi(k)+hybi(k+1))
        write(6, 111) k, hyam(k), hybm(k), 1000.*(hyam(k)+hybm(k)), -7*alog(hyam(k)+hybm(k))
111   format(I3, 4(2x, F15.9))
      enddo
      end subroutine read_vgrid_gfs

      SUBROUTINE get_p3d(Ps, im, jm, nz,   hyam, hybm,  P3d)
       integer :: im, jm, nz
       real  :: ps(im, jm), p3d(im, jm, nz)
       real  ::  hyam(nz), hybm(nz)
       integer :: i, j,k
      do k=1, nz
          p3d(:,:,k) = 1.e5*hyam(k)+ ps(:,:)*hybm(k)
      enddo

      end SUBROUTINE get_p3d
      SUBROUTINE get_delp(Ps, im, jm, nz,   hyai, hybi,  delp)
       integer :: im, jm, nz
       real  :: ps(im, jm), delp(im, jm, nz)
       real  ::  hyai(nz+1), hybi(nz+1)
       integer :: i, j,k
      do k=1, nz
          delp(:,:,k) = 1.e5*(hyai(k+1)-hyai(k))+ ps(:,:)*(hybi(k+1)-hybi(k))
      enddo

      end SUBROUTINE get_delp
