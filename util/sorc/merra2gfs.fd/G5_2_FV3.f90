program GEOS72_2_FV3128
!
! ifort -check bounds -traceback wam_grids_mod.f90 regrid_wamt62_orig2013.f90 wam_mapping.f90 Wam_geos5_das_GEOS_FV3.f90 -132 \ 
!-L/apps/netcdf/4.3.0-intel/lib -lnetcdf -lnetcdff -I/apps/netcdf/4.3.0-intel/include
! ifort  wam_grids_mod.f90 regrid_wamt62.f90 wam_mapping.f90 Wam_geos5_das.f90 -132
!
!GEOS72_2_FV3128 compile: g5gv3_ifort
! run ./a.out

  use fv3wam_grids_module, only : read_fv3wamgrids
  implicit none
!========================================================================================
! /glade/u/home/karol/GEOS_REGRID/Vmain_merra.f90
!  MERRA-orig:  /glade/p/cesm/chwg_dev/met_data/MERRA/0.5x0.6/2016
!  ls /scratch4/NCEPDEV/nems/noscrub/emc.nemspara/FV3GFS_V0_RELEASE/ICs/C96_2016081200
!
! Steps to make ICs from  ...."GEOS.fp.asm.tavg3_3d_asm_Nv.20160125_2230.V01.nc4"
!                           ftp://ftp.nccs.nasa.gov/fp/forecast/Y2018/M01/D04/H18/
! mv gfs_data.tile*nc /scratch4/NCEPDEV/stmp4/Valery.Yudin/Z90L128_FV30125/C96_2016081200
! cp gfs_ctrl.nc /scratch4/NCEPDEV/stmp4/Valery.Yudin/Z90L128_FV30125/C96_2016081200
!========================================================================================
  integer :: yearvar
  integer :: monthvar
  integer :: dayvar
  character(len=4)  ::   CRES_FV3='C096'
  character(12)     ::   finp ='file_res'
  integer           ::   levgfs
  integer           ::   fdate
!  character(len=4), parameter ::   CRES_FV3='C768'
!  character(len=4), parameter ::   CRES_FV3='C384'
  integer,          parameter ::   nzgfs =  128

  integer,          parameter ::   nzi =  72
  integer,          parameter ::   nyi = 361
  integer,          parameter ::   nxi = 576
  !integer,          parameter ::   nyi = 721
  !integer,          parameter ::   nxi = 1152
  integer          :: IMT, JMT

  real, dimension(nzi)   :: hyam, hybm
  real, dimension(nzi+1) :: hyai, hybi

  real, dimension(nzgfs)   :: ghyam, ghybm
  real, dimension(nzgfs+1) :: ghyai, ghybi

  character(len=128) :: vg_g5, vg_gfs

  real, dimension (:,:), allocatable :: GEOLON, GEOLON_W, GEOLON_S, &
                                        GEOLAT, GEOLAT_W, GEOLAT_S
!
  real, dimension(nxi, nyi, nzi)   :: U, V, T, Omega, ZH, Delp, O3, QI, QL, QV
  real, dimension(nxi, nyi     )   :: PS, PHIS
  real, dimension(nxi, nyi, nzi+1) :: ZHP, ZH2
!
!
  real, allocatable, dimension(:,:,:) :: Tic, Uic, Vic, Omic, Zhic
  real, allocatable, dimension(:,:,:) :: Usic, Vsic, Uwic, Vwic
  real, allocatable, dimension(:,:,:) :: O3ic, Qiic, Qlic, Qvic
  real, allocatable, dimension(:,:)   ::   Psic, Phisic, PsicW, PsicS

  real, allocatable, dimension(:,:,:) :: gTic, gUic, gVic, gOmic, gZhic	
  real, allocatable, dimension(:,:,:) :: gUsic, gVsic, gUwic, gVwic
  real, allocatable, dimension(:,:,:) :: gO3ic, gQiic, gQlic, gQvic

  real, allocatable, dimension(:,:,:) ::  P3d, P3dS,  P3dW
  real, allocatable, dimension(:,:,:) :: gP3d, gP3dS, gP3dW, gdelp

  real    ::  rlon(nxi), rlat(nyi)
  integer :: isd, ied, jsd, jed
  integer :: iso, ieo, jso, jeo
  integer, allocatable :: id1(:,:), id2(:,:), jdc(:,:)
  real,    allocatable ::  s2c(:,:, :)

  integer, allocatable :: id1_w(:,:), id2_w(:,:), jdc_w(:,:)
  real,    allocatable ::  s2c_w(:,:, :)

  integer, allocatable :: id1_s(:,:), id2_s(:,:), jdc_s(:,:)
  real,    allocatable ::  s2c_s(:,:, :)



  real,    allocatable :: lon_tile(:,:), lat_tile(:,:) 
  integer  :: i, j, k
!
  character(len=16) :: year, month, day
!  character(len=128), parameter :: path_in = '/glade/p/cesm/chwg_dev/met_data/GEOS5/0.9x1.25/2017'
!  character(len=128), parameter :: path_in = '/glade/p/cesm/chwg_dev/met_data/GEOS5/2017'   !!!!!!!!!!!!from march2017
!  character(len=128), parameter :: path_in = '/glade/p/cesm/chwg_dev/met_data/GEOS5/0.5x0.6/2016'   !!!!!for 02/29/2016
!  character(len=128), parameter :: path_in = '/glade/p/cesm/chwg_dev/met_data/GEOS5/0.5x0.6/2016'   !!!!!for 02/29/2016

  character(len=128), parameter :: path_in = '/scratch3/NCEPDEV/swpc/noscrub/Valery.Yudin/GEOS5_2018/GEOS_FV3'
  character(len=128), parameter :: path_out = '/scratch3/NCEPDEV/swpc/noscrub/Valery.Yudin/GEOS5_2018/GEOS_FV3'
  character(len=128) :: DIRres, Dirout, Ftile_out
  character(len=1  ) :: XNC
  integer, parameter :: newLats = 94 
  integer, parameter :: newLons = 192 
  integer :: mdays(12)
  character(len=256) :: file_tile_grids
  character(len=256) :: file_vcoord

  character(len=256) :: inFileName, outFileName
  character(len=256) :: format_year,format_day,format_month
!  openw, 1, 'file_res'
!  printf, 1, cres =>CRES_FV3
!  printf, 1, levgfs
!  printf, 1, fdate
!  printf, 1, inFileName
!  inFileName='/scratch3/NCEPDEV/swpc/noscrub/Valery.Yudin/GEOS5_2018/GEOS.fp.asm.inst3_3d_asm_Nv.20180201_0000.V01.nc4'
!
! 
  inFileName='/scratch3/NCEPDEV/swpc/noscrub/Valery.Yudin/GEOS5_2016/MERRA2_400.inst3_3d_asm_Nv.20160812.nc4'
   open(unit=66, file=trim(finp), form='formatted', status='unknown')
   read(66, *) CRES_FV3
   read(66, *) levgfs
   read(66, *) fdate
   read(66, 661) inFileName
661 format(A256)   
   close(66)
   print *, CRES_FV3,  ' CRES_FV3 '
   print *, levgfs, fdate, ' levgfs/fdate ' 
   print *, inFileName , ' G5_FileName '
   print *
 mdays(1:12)=(/ 31,29,31,30,31,30,31,31,30,31,30,31 /)
! loop over years,months, days
  print*,' before read_wamgrids' 
!
!C1152  C192  C3072  C384  C48  C768  C768_gfdl	C96/scratch4/NCEPDEV/global/save/glopara/svn/fv3gfs/trunk/global_shared.v15.0.0/fix/fix_fv3
!
  file_vcoord=trim('gfs_ctrl.nc')
 if (CRES_FV3.eq.'C096') then
   IMT =96
   JMT =96
  file_tile_grids=trim('/scratch4/NCEPDEV/global/save/glopara/svn/fv3gfs/trunk/global_shared.v15.0.0/fix/fix_fv3/C96/C96_grid.tile1.nc')
  DirRes ='/scratch4/NCEPDEV/global/save/glopara/svn/fv3gfs/trunk/global_shared.v15.0.0/fix/fix_fv3/C96/C96_grid.tile'
 endif

 if (CRES_FV3.eq.'C192') then
   IMT =192
   JMT =192
  file_tile_grids=trim('/scratch4/NCEPDEV/global/save/glopara/svn/fv3gfs/trunk/global_shared.v15.0.0/fix/fix_fv3/C192/C192_grid.tile1.nc')
  DirRes ='/scratch4/NCEPDEV/global/save/glopara/svn/fv3gfs/trunk/global_shared.v15.0.0/fix/fix_fv3/C192/C192_grid.tile'
 endif

 if (CRES_FV3.eq.'C384') then
   IMT =384	
   JMT =384
  file_tile_grids=trim('/scratch4/NCEPDEV/global/save/glopara/svn/fv3gfs/trunk/global_shared.v15.0.0/fix/fix_fv3/C384/C384_grid.tile1.nc')
  DirRes ='/scratch4/NCEPDEV/global/save/glopara/svn/fv3gfs/trunk/global_shared.v15.0.0/fix/fix_fv3/C384/C384_grid.tile'
 endif

 if (CRES_FV3.eq.'C768') then
   IMT =768
   JMT =768
  file_tile_grids=trim('/scratch4/NCEPDEV/global/save/glopara/svn/fv3gfs/trunk/global_shared.v15.0.0/fix/fix_fv3/C768/C768_grid.tile1.nc')
  DirRes ='/scratch4/NCEPDEV/global/save/glopara/svn/fv3gfs/trunk/global_shared.v15.0.0/fix/fix_fv3/C768/C768_grid.tile'
 endif


   ALLOCATE(GEOLON(IMT,JMT), GEOLAT(IMT,JMT))
   ALLOCATE(GEOLON_W(IMT+1,JMT), GEOLAT_W(IMT+1,JMT))
   ALLOCATE(GEOLON_S(IMT,JMT+1), GEOLAT_S(IMT,JMT+1))

  CALL READ_TILE_GRIDS(FILE_TILE_GRIDS, IMT, JMT, GEOLON, GEOLON_W, GEOLON_S, &
                       GEOLAT, GEOLAT_W, GEOLAT_S)

  CALL read_vgrid_gfs(vg_gfs, nzgfs, ghyam, ghybm, ghyai, ghybi)
  CALL read_vgrid_geos5(vg_g5, nzi,  hyam,  hybm,    hyai, hybi)


 
!  print *,  'GEOLAT: ', GEOLAT *180./3.14 
!  print *,  'GEOLON: ', GEOLON *180./3.14 
!  stop
!
! read the grid information  "fv_climate_nudge_init" & register diagnostics
!
!     allocate ( lon_obs(nlon_obs), lat_obs(nlat_obs), ak_obs(nlev_obs+1), bk_obs(nlev_obs+1) )
!     call read_grid ( lon_obs, lat_obs, ak_obs, bk_obs )

!  call read_fv3wamgrids(CRES_FV3)



  call Read_nc4_omega(inFileName, nxI, nyI, nzI, ps, phis, U, T,V, Omega, ZH, Delp, O3, QI, QL, QV, rlon, rlat) 
  
  
  do k=2, nzi+1
    zhp(:,:, k) = ZH(:,:, k-1)
  enddo
    zhp(:,:,1) = Phis(:,:)/9.81
  print *
  do  i=1, nxi, 16

  print *, i, rlon(i)*180./3.14
  enddo
  print * 

!=========================================
!  Do Vertical Interpolation "U, T, V, QV "
!           recompute ZH
!=========================================

   CALL  compute_zh_geos5(nxi, nyi, nzi, delp, ps, phis, T, QV, zh2)
   print *, maxval(Zhp),   minval(Zhp),   ' Zhp '
   print *, maxval(Zh),   minval(Zh),      ' Zh '
   print *, maxval(Zh2),   minval(Zh2),      ' Zh2 '
   print *, maxval(Zh2(1:nxi, 1:nyi, 2)),   minval(Zh2(1:nxi, 1:nyi, 2)),      ' Zh2-2 '
   print *, maxval(Zh2(1:nxi, 1:nyi, 1)),   minval(Zh2(1:nxi, 1:nyi, 1)),      ' Zh2-1 '
   
  
  isd =1
  ied =nxi
  jsd =1
  jed =nyi
!
  iso=1
  ieo=imt
  jso=1
  jeo=jmt
! 
  allocate(id1(imt, jmt))
  allocate(id2(imt, jmt))
  allocate(jdc(imt, jmt))
  allocate(s2c(imt, jmt,4))


  allocate(Psic(imt, jmt))
  allocate(PsicW(imt+1, jmt))
  allocate(PsicS(imt, jmt+1))
  allocate(Phisic(imt, jmt))

  allocate(id1_w(imt+1, jmt))
  allocate(id2_w(imt+1, jmt))
  allocate(jdc_w(imt+1, jmt))
  allocate(s2c_w(imt+1, jmt,4))

  allocate(id1_s(imt, jmt+1))
  allocate(id2_s(imt, jmt+1))
  allocate(jdc_s(imt, jmt+1))
  allocate(s2c_s(imt, jmt+1,4))




  allocate(Tic(imt, jmt, nzi))
  allocate(O3ic(imt, jmt, nzi))
  allocate(Qvic(imt, jmt, nzi))
  allocate(QLic(imt, jmt, nzi))
  allocate(QIic(imt, jmt, nzi))
  allocate(OMic(imt, jmt, nzi))
!
  allocate(ZHic(imt, jmt, nzi+1))
  allocate(Uic(imt, jmt, nzi))
  allocate(Vic(imt, jmt, nzi))

  allocate(Usic(imt, jmt+1, nzi))
  allocate(Vsic(imt, jmt+1, nzi))
  allocate(Uwic(imt+1, jmt, nzi))
  allocate(Vwic(imt+1, jmt, nzi))


  allocate(gTic(imt, jmt, nzgfs))
  allocate(gO3ic(imt, jmt, nzgfs))
  allocate(gQvic(imt, jmt, nzgfs))
  allocate(gQLic(imt, jmt, nzgfs))
  allocate(gQIic(imt, jmt, nzgfs))
  allocate(gOMic(imt, jmt, nzgfs))
!
  allocate(gDelp(imt, jmt, nzgfs))
  allocate(gZHic(imt, jmt, nzgfs+1))
  allocate(gUic(imt, jmt, nzgfs))
  allocate(gVic(imt, jmt, nzgfs))

  allocate(gUsic(imt, jmt+1, nzgfs))
  allocate(gVsic(imt, jmt+1, nzgfs))
  allocate(gUwic(imt+1, jmt, nzgfs))
  allocate(gVwic(imt+1, jmt, nzgfs))

  allocate(gP3dS(imt, jmt+1, nzgfs))
  allocate(gP3dW(imt+1, jmt, nzgfs))
  allocate(gP3d(imt, jmt, nzgfs))

  allocate(P3dS(imt, jmt+1, nzi))
  allocate(P3dW(imt+1, jmt, nzi))
  allocate(P3d(imt, jmt, nzi))



!uw_ic, us_ic, vw_ic, vs_ic)





  DirOUT ='gfs_data.tile'



TILE3d: DO i=1, 6

  write(Xnc, "(I1)") I
  print *, XNC
  FILE_TILE_GRIDS= trim(DirRes)//trim(XNC)//trim('.nc')
  print *, FILE_TILE_GRIDS

  CALL READ_TILE_GRIDS(FILE_TILE_GRIDS, IMT, JMT, GEOLON, GEOLON_W, GEOLON_S, &
                       GEOLAT, GEOLAT_W, GEOLAT_S)

  CALL remap_coef( isd, ied, jsd, jed, rlon, rlat,           & 
                   iso, ieo, jso, jeo, geolon, geolat,   &
                   id1, id2, jdc, s2c )  

  CALL  remap_xy_2d( isd, ied, jsd, jed, Ps,       & 
                          iso, ieo, jso, jeo, id1, id2, jdc, s2c, Psic)
  CALL  remap_xy_2d( isd, ied, jsd, jed, Phis,       & 
                          iso, ieo, jso, jeo, id1, id2, jdc, s2c, Phisic)

  CALL  remap_xy_3d( isd, ied, jsd, jed, NzI, T,       & 
                          iso, ieo, jso, jeo, id1, id2, jdc, s2c, Tic)
   print *, maxval(T),   minval(T),   ' T-glob '
   print *, maxval(Tic), minval(Tic), ' T-tile '


  CALL  remap_xy_3d( isd, ied, jsd, jed, NzI, O3,       & 
                          iso, ieo, jso, jeo, id1, id2, jdc, s2c, O3ic)
  CALL  remap_xy_3d( isd, ied, jsd, jed, NzI, Qv,       & 
                          iso, ieo, jso, jeo, id1, id2, jdc, s2c, Qvic)
  CALL  remap_xy_3d( isd, ied, jsd, jed, NzI, QL,       & 
                          iso, ieo, jso, jeo, id1, id2, jdc, s2c, QLic)
  CALL  remap_xy_3d( isd, ied, jsd, jed, NzI, Qi,       & 
                          iso, ieo, jso, jeo, id1, id2, jdc, s2c, QIic)
  CALL  remap_xy_3d( isd, ied, jsd, jed, NzI, Omega,       & 
                          iso, ieo, jso, jeo, id1, id2, jdc, s2c, OMic)

!Uic, Vic, Zhic

  CALL  remap_xy_3d( isd, ied, jsd, jed, NzI+1, ZHP,       & 
                          iso, ieo, jso, jeo, id1, id2, jdc, s2c, ZHic)
!
!                                      U_S     V_S
!  CALL GL2ANYV(0,LEVSO,U,V,LONB,LATB,CUBE_3D,CUBE_3D2,IM,(JM+1),GEOLON_S, GEOLAT_S)
!  CALL GL2ANYV(0,LEVSO,U,V,LONB,LATB,CUBE_3D,CUBE_3D2,(IM+1),JM,GEOLON_W, GEOLAT_W)
!                                      U_W     V_W
!
  CALL remap_coef( isd, ied, jsd, jed, rlon, rlat,           & 
                   iso, ieo, jso, jeo+1, geolon_s, geolat_s,   &
                   id1_s, id2_s, jdc_s, s2c_s )

  CALL remap_coef( isd, ied, jsd, jed, rlon, rlat,           & 
                   iso, ieo+1, jso, jeo, geolon_w, geolat_w,   &
                   id1_w, id2_w, jdc_w, s2c_w )
 


  CALL  remap_xy_3d( isd, ied, jsd, jed, NzI, U,       & 
                          iso, ieo, jso, jeo, id1, id2, jdc, s2c, Uic)
  CALL  remap_xy_3d( isd, ied, jsd, jed, NzI, V,       & 
                          iso, ieo, jso, jeo, id1, id2, jdc, s2c, Vic)

  CALL  remap_xy_3d( isd, ied, jsd, jed, NzI, U,       & 
                          iso, ieo+1, jso, jeo, id1_w, id2_w, jdc_w, s2c_w, UWic)

  CALL  remap_xy_3d( isd, ied, jsd, jed, NzI, V,       & 
                          iso, ieo+1, jso, jeo, id1_w, id2_w, jdc_w, s2c_w, VWic)



  CALL  remap_xy_3d( isd, ied, jsd, jed, NzI, U,       & 
                          iso, ieo, jso, jeo+1, id1_s, id2_s, jdc_s, s2c_s, USic)

  CALL  remap_xy_3d( isd, ied, jsd, jed, NzI, V,       & 
                          iso, ieo, jso, jeo+1, id1_s, id2_s, jdc_s, s2c_s, VSic)
!
! PS-shifts
!
  CALL  remap_xy_2d( isd, ied, jsd, jed, Ps,       & 
                          iso, ieo, jso, jeo+1, id1_s, id2_s, jdc_s, s2c_s, PsicS)
  CALL  remap_xy_2d( isd, ied, jsd, jed, Ps,       & 
                          iso, ieo+1, jso, jeo, id1_w, id2_w, jdc_w, s2c_w, PsicW)

  call get_p3d(Psic, imt, jmt, nzi,   hyam, hybm,    P3d)
  call get_p3d(Psic, imt, jmt, nzgfs, ghyam, ghybm, gP3d)
  call get_delp(Psic, imt, jmt, nzgfs,  ghyai, ghybi, gDelp)

  call get_p3d(PsicW, imt+1, jmt, nzi,   hyam, hybm,    P3dW)
  call get_p3d(PsicW, imt+1, jmt, nzgfs, ghyam, ghybm, gP3dW)

  call get_p3d(PsicS, imt, jmt+1, nzi,   hyam, hybm,    P3dS)
  call get_p3d(PsicS, imt, jmt+1, nzgfs, ghyam, ghybm, gP3dS)


  print *, maxval(U), minval(U), ' U-geos '
  print *, maxval(Usic), minval(Usic), ' USic-geos '
  print *, maxval(Uwic), minval(Uwic), ' UWic-geos '
  print *, maxval(Uic), minval(Uic),    ' Uic-geos '
!
! write 3D-tile output for IC-128L
!
  Ftile_out=trim(DirOut)//trim(XNC)//trim('.nc')
!                 (im, jm, nzi, A, Pa, nzg, gPa, gA)

!
! Vertical Remap
!
  call remap_zdim(imt, jmt, nzi, Tic,  P3d, nzgfs, gp3d,   gTic)
  call remap_zdim(imt, jmt, nzi, QVic, P3d, nzgfs, gp3d,  gQVic)
!
! Recompute gZH
!
  CALL  compute_zh_geos5(imt, jmt, nzgfs, Gdelp, psic, phisic, gTic, gQVic, gzhic)	
  print *, maxval(gzhic), minval(gzhic), ' gzhic '

! OMEGA + Tracers
  call remap_zdim(imt, jmt, nzi, Omic, P3d, nzgfs, gp3d,  gOmic)
  
  call remap_zdim(imt, jmt, nzi, QLic, P3d, nzgfs, gp3d,  gQLic)
  call remap_zdim(imt, jmt, nzi, QIic, P3d, nzgfs, gp3d,  gQIic)
  call remap_zdim(imt, jmt, nzi, O3ic, P3d, nzgfs, gp3d,  gO3ic)
!
! West-shit
!
  call remap_zdim(imt+1, jmt, nzi, UWic, P3dW, nzgfs, gp3dW, gUWic)
  call remap_zdim(imt+1, jmt, nzi, VWic, P3dW, nzgfs, gp3dW, gVWic)
! South-shift
  call remap_zdim(imt, jmt+1, nzi, USic, P3dS, nzgfs, gp3dS, gUSic)
  call remap_zdim(imt, jmt+1, nzi, VSic, P3dS, nzgfs, gp3dS, gVSic)


!  CALL WRITE_3DTILE_GFS(Ftile_out, imt, jmt, nzi, geolat, geolon, psIC, omIC, zhIC, QvIC, QLIC, O3IC, uwic, usic, vwic, vsic) 

  CALL WRITE_3DTILE_GFS(Ftile_out, imt, jmt, nzGFS, geolat, geolon, psIC, GomIC, GzhIC, GQvIC, GQLIC, GQIic, GO3IC, &
       Guwic, Gusic, Gvwic, Gvsic)   
   print *, Ftile_out
   ENDDO TILE3D
   
   
  deallocate(psic, phisic, psicS,psicW)
  deallocate(p3d, gp3d,  p3dW, gp3dW,p3dS, gp3dS)

  deallocate(geolat, geolon, geolat_w, geolon_w,geolat_s, geolon_s)
  deallocate(id1, id2, jdc, s2c)
  deallocate(id1_w, id2_w, jdc_w, s2c_w)
  deallocate(id1_s, id2_s, jdc_s, s2c_s)

  deallocate(Tic, Uic, Vic, Omic, Zhic)
  deallocate(Usic, Vsic, Uwic, Vwic)
  deallocate(O3ic, Qvic, QLic, QIic)

  deallocate(gTic, gUic, gVic, gOmic, gZhic)
  deallocate(gUsic, gVsic, gUwic, gVwic)
  deallocate(gO3ic, gQvic, gQLic, gQIic)


  STOP 'Done-FV3GFS_IC with GEOS5 721x1152 72L '


end program GEOS72_2_FV3128



subroutine handle_err( status )
  implicit none

  integer, intent(in) :: status

  print*,'***** handle_err status = ',status
  stop

end subroutine handle_err
