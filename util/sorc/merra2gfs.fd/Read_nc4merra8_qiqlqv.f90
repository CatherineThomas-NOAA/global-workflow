
  subroutine read_nc4_omega(inFileName, nx, ny, nz, ps, phis, U, T,V, Omega, ZH, Delp, O3, QI, QL, QV, lon2, lat2)

  use netcdf
  use fv3wam_grids_module, only :  newLat, newLon
  use fv3wam_grids_module, only :  nlat, nlon, nlev, ntime
  use fv3wam_grids_module, only :  newlev,newhyam,newhybm
  use fv3wam_grids_module, only :  newlevi,newhyai,newhybi
 
  implicit none 
 
  character(len=256) :: inFileName
  integer         ::  nx, ny, nz
  real, dimension(nx, ny, nz) :: U, V, T, Omega, ZH, Delp, O3, QI, QL, QV
  real, dimension(nx, ny    ) :: PS, PHIS
  real             :: lat2(ny), lon2(nx)

  integer :: out_ncid
  integer :: ncid, status, iLevDimID, LevDimID, LonDimID, LatDimID, RecordDimID, DateDimID
  integer :: niLevs, nLevs, nLats, nLons, nRecords,nDate
  character(len = nf90_max_name) :: RecordDimName

  integer :: ilev, dpid
  integer :: psid, tsid, tid, qid,dateid,datesecid,latid,lonid,levid,hyamid,hybmid,hyaiid,hybiid,timeid
  integer :: phiid, landid, qflxid, hflxid, snowhid, fsdsid, soilwid, tauxid, tauyid, uid, vid

  real, allocatable :: oldLat(:), oldLon(:)

  real, allocatable :: lev(:)
  real, allocatable :: hyam(:)
  real, allocatable :: hybm(:)
  real, allocatable :: hyai(:)
  real, allocatable :: hybi(:)
  real, parameter   :: radfac = 3.1415926535897931/180.
  integer, parameter :: IREC=1

  integer :: i, j, k
  real    :: rx, ry, gx, scal
  print*,' inFileName = '//trim(inFileName)
  
   
!  call date_and_time( values=time_vals)
!  write(*,fmt='("time: ",I2.2,":",I2.2,":",I2.2)') time_vals(5:7)

  status = nf90_open(inFileName, nf90_nowrite, ncid)
  if (status /= nf90_noerr) call handle_err(status)

  ! Get ID of unlimited dimension
  status = nf90_inquire(ncid, unlimitedDimId = RecordDimID)

  if (status /= nf90_noerr) call handle_err(status)



  status = nf90_inq_dimid(ncid, "lev", LevDimID)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_dimid(ncid, "lat", LatDimID)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_dimid(ncid, "lon", LonDimID)
  if (status /= nf90_noerr) call handle_err(status)




  ! How many values of "lev" are there?
  status = nf90_inquire_dimension(ncid, LevDimID, len = nLevs)
  ! 3D scalar fields
  if (status /= nf90_noerr) call handle_err(status)

  ! How many values of "lat" are there?
  status = nf90_inquire_dimension(ncid, LatDimID, len = nLats)
  if (status /= nf90_noerr) call handle_err(status)

  ! How many values of "lon" are there?
  status = nf90_inquire_dimension(ncid, LonDimID, len = nLons)
  if (status /= nf90_noerr) call handle_err(status)

  ! What is the name of the unlimited dimension, how many records are there?
  status = nf90_inquire_dimension(ncid, RecordDimID, name = RecordDimName, len = nRecords)
  if (status /= nf90_noerr) call handle_err(status)
!
! Records in file
!  
  allocate(lev(nLevs),hyam(nLevs),hybm(nLevs))
  allocate(hyai(niLevs),hybi(niLevs))
!  allocate(time(nRecords))
  allocate(oldLon(nLons))
  allocate(oldLat(nLats))

!lon-lat-lev-time


  status = nf90_inq_varid(ncid, 'time', timeid)
  if (status /= nf90_noerr) call handle_err(status)
!  status = nf90_get_var(ncid, timeid, time)
!  if (status /= nf90_noerr) call handle_err(status)
  !print*,time


  status = nf90_inq_varid(ncid, 'lev', levid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, levid, lev)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, 'lon', lonid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_inq_varid(ncid, 'lat', latid)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_get_var(ncid, lonid, oldLon)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, latid, oldLat)
  if (status /= nf90_noerr) call handle_err(status)


  do i=1, nlons
     gx = oldLon(i)
      if (gx .lt. 0) gx = gx + 360.0
     lon2(i) = radfac*gx

  enddo
  do j=1, nlats
     lat2(j) = radfac*oldLat(j)
  enddo
!
! do we need FLIP of Lons ??
!
  status = nf90_inq_varid(ncid, 'T', tid)
  if (status /= nf90_noerr) call handle_err(status)
!  status = nf90_get_var(ncid, tid, T)
  status= nf90_get_var( ncid, tid, T, start=(/1,1,1,IREC/), count=(/nx,ny,nz,1/))

  status = nf90_inq_varid(ncid, 'QV', qid)
  if (status /= nf90_noerr) call handle_err(status)
!  status = nf90_get_var(ncid, qid, QV)
  status= nf90_get_var( ncid, qid, QV, start=(/1,1,1,IREC/), count=(/nx,ny,nz,1/))  

  status = nf90_inq_varid(ncid, 'QL', qid)
  if (status /= nf90_noerr) call handle_err(status)
!  status = nf90_get_var(ncid, qid, QL)
  status= nf90_get_var( ncid, qid, QL, start=(/1,1,1,IREC/), count=(/nx,ny,nz,1/))  
  
  status = nf90_inq_varid(ncid, 'QI', qid)
  if (status /= nf90_noerr) call handle_err(status)
!  status = nf90_get_var(ncid, qid, QI)
  status= nf90_get_var( ncid, qid, QI, start=(/1,1,1,IREC/), count=(/nx,ny,nz,1/)) 
   
  status = nf90_inq_varid(ncid, 'O3', qid)
  if (status /= nf90_noerr) call handle_err(status)
!  status = nf90_get_var(ncid, qid, O3)
   status= nf90_get_var( ncid, qid, O3, start=(/1,1,1,IREC/), count=(/nx,ny,nz,1/))  

  status = nf90_inq_varid(ncid, 'PHIS', phiid)
  if (status /= nf90_noerr) call handle_err(status)
!  status = nf90_get_var(ncid, phiid, PHIS)
  status= nf90_get_var( ncid, phiid, PHIS, start=(/1,1,IREC/), count=(/nx,ny,1/))    

  status = nf90_inq_varid(ncid, 'PS', psid)
  if (status /= nf90_noerr) call handle_err(status)
!  status = nf90_get_var(ncid, psid, PS)
  status= nf90_get_var( ncid, psid, PS, start=(/1,1,IREC/), count=(/nx,ny,1/))  

  status = nf90_inq_varid(ncid, 'U', uid)
  if (status /= nf90_noerr) call handle_err(status)
!  status = nf90_get_var(ncid, uid, U)
  status= nf90_get_var( ncid, Uid, U, start=(/1,1,1,IREC/), count=(/nx,ny,nz,1/))

  status = nf90_inq_varid(ncid, 'V', vid)
  if (status /= nf90_noerr) call handle_err(status)
!  status = nf90_get_var(ncid, vid, V)
  status= nf90_get_var( ncid, vid, V, start=(/1,1,1,IREC/), count=(/nx,ny,nz,1/))

  status = nf90_inq_varid(ncid, 'OMEGA', vid)
  if (status /= nf90_noerr) call handle_err(status)
!  status = nf90_get_var(ncid, vid, OMEGA)
  status= nf90_get_var( ncid, vid, OMEGA, start=(/1,1,1,IREC/), count=(/nx,ny,nz,1/))
  
 status = nf90_inq_varid(ncid, 'H', vid)
  if (status /= nf90_noerr) call handle_err(status)
!  status = nf90_get_var(ncid, vid, ZH)
  status= nf90_get_var( ncid, vid, ZH, start=(/1,1,1,IREC/), count=(/nx,ny,nz,1/))  

  status = nf90_inq_varid(ncid, 'DELP',dpid)
  if (status /= nf90_noerr) call handle_err(status)
!  status = nf90_get_var(ncid, dpid, delp)
   status= nf90_get_var( ncid, dpid, delp, start=(/1,1,1,IREC/), count=(/nx,ny,nz,1/))   
!
! read hyam/hybm hyai/hybi GEOS5-72L
!
  call Flip1d_g5_to360(lon2, nx)
  call Flip2d_g5_to360(ps, nx, ny)
  call Flip2d_g5_to360(phis, nx, ny)

  call Flip3d_g5_to360(T, nx, ny, nz)
  call Flip3d_g5_to360(U, nx, ny, nz)
  call Flip3d_g5_to360(V, nx, ny, nz)
  call Flip3d_g5_to360(OMEGA, nx, ny, nz)
  call Flip3d_g5_to360(ZH, nx, ny, nz)
  call Flip3d_g5_to360(DELP, nx, ny, nz)

  call Flip3d_g5_to360(O3, nx, ny, nz)
  call Flip3d_g5_to360(QV, nx, ny, nz)
  call Flip3d_g5_to360(QL, nx, ny, nz)
  call Flip3d_g5_to360(QI, nx, ny, nz)
  print *, maxval(Phis)/9.81, minval(Phis)/9.81,  ' PHIS-minmax '
  print *, maxval(ZH), minval(Zh),  ' ZH-minmax '
  print *, maxval(OMEGA), minval(OMEGA),  ' OMEGA-minmax '
  print *, maxval(T), minval(T),  ' TEMP-minmax '
  print *, maxval(U), minval(U),  ' UWIND-minmax '
  print *, maxval(QV), minval(QV),  ' QV-minmax '
   print *,  lon2(1), lon2(nx)
   lon2(1) =0.0
!     print * 
!   do i=1, nx
!    print *, i, lon2(i)
!   enddo
!    print *
  scal  = 29./48.*1.e6
  print *, maxval(O3)*scal, minval(O3)*scal,  ' O3-minmax '
!  stop
  end subroutine read_nc4_omega
!
 SUBROUTINE READ_TILE_GRIDS(TILEFILE, IM, JM, GEOLON, GEOLON_W, GEOLON_S, &
    GEOLAT, GEOLAT_W, GEOLAT_S)
  use netcdf
  implicit NONE
  real, parameter   :: radfac = 3.1415926535897931/180.
 CHARACTER(len=*) :: TILEFILE
!
 INTEGER :: IM, JM
 REAL :: GEOLON(IM,JM)
 REAL :: GEOLON_W(IM+1,JM)
 REAL :: GEOLON_S(IM,JM+1)
!
 REAL :: GEOLAT(IM,JM)
 REAL :: GEOLAT_W(IM+1,JM)
 REAL :: GEOLAT_S(IM,JM+1)
!
 real, ALLOCATABLE     :: TMPVAR(:,:)
 integer :: IMT, JMT, NX, NY
 integer :: NCID, ERROR, ID_VAR, ID_DIM
!
 ERROR=NF90_OPEN(TRIM(TILEFILE),NF90_NOWRITE,NCID)
! CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(TILEFILE) )

 ERROR=NF90_INQ_DIMID(NCID, 'nx', ID_DIM)
! CALL NETCDF_ERR(ERROR, 'ERROR READING NX ID' )

 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NX)
! CALL NETCDF_ERR(ERROR, 'ERROR READING NX' )

 ERROR=NF90_INQ_DIMID(NCID, 'ny', ID_DIM)
! CALL NETCDF_ERR(ERROR, 'ERROR READING NY ID' )

 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NY)
! CALL NETCDF_ERR(ERROR, 'ERROR READING NY' )

 IF (MOD(NX,2) /= 0) THEN
   PRINT*,'FATAL ERROR: NX IS NOT EVEN'
 ENDIF

 IF (MOD(NY,2) /= 0) THEN
   PRINT*,'FATAL ERROR: NY IS NOT EVEN'
 ENDIF

 IMT = NX/2
 JMT = NY/2

 if (IMT .ne. IM .or. JMT .ne. JM) then
    print *, 'error in GRID_TILE FILE DIMENSIONS '
    print *, IM, IMT, TRIM(TILEFILE)
    print *, JM, JMT, TRIM(TILEFILE)
   stop  'error in GRID_TILE of READ_TILE_GRIDS'
 endif

 PRINT*, "READ FV3 GRID INFO FROM: "//TRIM(TILEFILE)

 ALLOCATE(TMPVAR(NX+1,NY+1))


 ERROR=NF90_INQ_VARID(NCID, 'x', ID_VAR) 
! CALL NETCDF_ERR(ERROR, 'ERROR READING X ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, TMPVAR)
! CALL NETCDF_ERR(ERROR, 'ERROR READING X RECORD' )

 GEOLON(1:IM,1:JM)     = TMPVAR(2:NX:2,2:NY:2)*radfac
 GEOLON_W(1:IM+1,1:JM) = TMPVAR(1:NX+1:2,2:NY:2)*radfac
 GEOLON_S(1:IM,1:JM+1) = TMPVAR(2:NX:2,1:NY+1:2)*radfac

 ERROR=NF90_INQ_VARID(NCID, 'y', ID_VAR) 
! CALL NETCDF_ERR(ERROR, 'ERROR READING Y ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, TMPVAR)
! CALL NETCDF_ERR(ERROR, 'ERROR READING Y RECORD' )

 ERROR = NF90_CLOSE(NCID)


 GEOLAT(1:IM,1:JM)     = TMPVAR(2:NX:2,2:NY:2)*radfac
 GEOLAT_W(1:IM+1,1:JM) = TMPVAR(1:NX+1:2,2:NY:2)*radfac
 GEOLAT_S(1:IM,1:JM+1) = TMPVAR(2:NX:2,1:NY+1:2)*radfac
 
 DEALLOCATE(TMPVAR)
 END SUBROUTINE READ_TILE_GRIDS
