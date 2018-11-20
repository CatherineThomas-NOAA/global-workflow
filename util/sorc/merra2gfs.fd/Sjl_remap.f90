
subroutine remap_coef( isd, ied, jsd, jed, lon_in, lat_in,   & 
                         is, ie, js, je, lon_out, lat_out,   &
                         id1, id2, jdc, s2c )
!--------
! Input:
!--------
! Data Input: data must be global
  integer, intent(in):: isd, ied                  ! Data x-dimension (W->E; must be periodic)
  integer, intent(in):: jsd, jed                  ! Data y-dimension (S->N)
  real,    intent(in):: lon_in(isd:ied)           ! Data longitude (Radian; periodic)
  real,    intent(in):: lat_in(jsd:jed)           ! Data latitude (increases from SP to NP)

! Model input:
  integer, intent(in):: is, ie, js, je        ! model horizontal dimensions (un-ghosted sub-domian)
  real,    intent(in):: lon_out(is:ie,js:je)  ! model longitude (Radian)
  real,    intent(in):: lat_out(is:ie,js:je)  ! model latitude (Radian)

!--------
! Output:
!--------
!
  integer, intent(out), dimension(is:ie,js:je  ):: id1, id2, jdc
  real,    intent(out), dimension(is:ie,js:je,4):: s2c
  real, parameter :: RADIUS = 6.3712e+6
  real, parameter ::         PI     = 3.1415926535897931
!real, public, parameter :: OMEGA  = 7.2921e-5 
!real, public, parameter :: GRAV   = 9.8066
!real, public, parameter :: RDGAS  = 287.05
!real, public, parameter :: RVGAS  = 461.50
! Extra:
!real, public, parameter :: HLV = 2.5e6_r8_kind   
!real, public, parameter :: HLF = 3.3358e5_r8_kind   
!real, public, parameter :: con_cliq   =4.1855e+3_r8_kind      ! spec heat H2O liq   (J/kg/K)
!real, public, parameter :: con_csol   =2.1060e+3_r8_kind      ! spec heat H2O ice   (J/kg/K) 
!===============================================================================================

! local:
  real:: rdlon(isd:ied)
  real:: rdlat(jsd:jed)
  real:: a1, b1
  integer i, j, i1, i2, jc, i0, j0

 !pk0(1) = ak_in(1)**KAPPA 
 !pn_top = log(ak_in(1))

  do i=isd,ied-1
     rdlon(i) = 1. / (lon_in(i+1) - lon_in(i))
  enddo
     rdlon(ied) = 1. / (lon_in(isd) + 2.*PI - lon_in(ied))  ! periodic assumption

  do j=jsd,jed-1
     rdlat(j) = 1. / (lat_in(j+1) - lat_in(j))
  enddo

! * Interpolate to cubed sphere cell center
  do 5000 j=js,je

     do i=is,ie

       if ( lon_out(i,j) .gt. lon_in(ied) ) then
            i1 = ied;     i2 = isd
            a1 = (lon_out(i,j)-lon_in(ied)) * rdlon(ied)
       elseif ( lon_out(i,j) .lt. lon_in(1) ) then
            i1 = ied;     i2 = isd
            a1 = (lon_out(i,j)+2.*PI-lon_in(ied)) * rdlon(ied)
       else
            do i0=isd,ied-1
            if ( lon_out(i,j) .ge. lon_in(i0) .and. lon_out(i,j) .le. lon_in(i0+1) ) then
               i1 = i0;  i2 = i0+1
               a1 = (lon_out(i,j)-lon_in(i1)) * rdlon(i0)
               go to 111
            endif
            enddo
       endif

111    continue

       if ( lat_out(i,j) .lt. lat_in(jsd) ) then
            jc = jsd
            b1 = 0.
       elseif ( lat_out(i,j) .gt. lat_in(jed) ) then
            jc = jed-1
            b1 = 1.
       else
          do j0=jsd,jed-1
          if ( lat_out(i,j) .ge. lat_in(j0) .and. lat_out(i,j) .le. lat_in(j0+1) ) then
               jc = j0
               b1 = (lat_out(i,j)-lat_in(jc)) * rdlat(jc)
               go to 222
          endif
          enddo
       endif
222    continue

! Debug codes:
!      if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
!           write(*,*) i,j,a1, b1
!      endif

       s2c(i,j,1) = (1.-a1) * (1.-b1)
       s2c(i,j,2) =     a1  * (1.-b1)
       s2c(i,j,3) =     a1  *     b1
       s2c(i,j,4) = (1.-a1) *     b1
       id1(i,j) = i1
       id2(i,j) = i2
       jdc(i,j) = jc

     enddo
5000 continue

  end subroutine remap_coef
!
  subroutine remap_xy_3d( isd, ied, jsd, jed, km, q_in,       & 
                          is, ie, js, je, id1, id2, jdc, s2c, &
                          q_out )
!--------
! Input:
!--------
! Data Input: data must be global
  integer, intent(in):: isd, ied                  ! Data x-dimension (W->E; must be periodic)
  integer, intent(in):: jsd, jed                  ! Data y-dimension (S->N)
  integer, intent(in):: km                        ! Data vertical-dimension
  real,    intent(in), dimension(isd:ied,jsd:jed,km) :: q_in    ! Data on input grid

! Model input:
  integer, intent(in):: is, ie, js, je        ! model horizontal dimensions (un-ghosted sub-domian)

! Indices and coefficients for bilinear interpolation
  integer, intent(in), dimension(is:ie,js:je  ):: id1, id2, jdc
  real,    intent(in), dimension(is:ie,js:je,4):: s2c

!--------
! Output:
!--------
! Data Output: on local model horizontal grid and input data vertical grid
  real,    intent(out), dimension(is:ie,js:je,km):: q_out

  integer :: i, j, k

    do k=1,km
    do j=js,je
    do i=is,ie
        q_out(i,j,k) = s2c(i,j,1)*q_in(id1(i,j),jdc(i,j),  k) + s2c(i,j,2)*q_in(id2(i,j),jdc(i,j),  k) +  &
                       s2c(i,j,3)*q_in(id2(i,j),jdc(i,j)+1,k) + s2c(i,j,4)*q_in(id1(i,j),jdc(i,j)+1,k)
    enddo
    enddo
    enddo

  end subroutine remap_xy_3d
!
!
 subroutine remap_xy_2d( isd, ied, jsd, jed, q_in,            & 
                          is, ie, js, je, id1, id2, jdc, s2c, &
                          q_out )
!--------
! Input:
!--------
  integer, intent(in):: isd, ied                            ! Data x-dimension (W->E; must be periodic)
  integer, intent(in):: jsd, jed                            ! Data y-dimension (S->N)
  real,    intent(in), dimension(isd:ied,jsd:jed) :: q_in   ! Data on input grid
  integer, intent(in):: is, ie, js, je                      ! model horizontal dimensions (un-ghosted sub-domian)
  integer, intent(in), dimension(is:ie,js:je  ):: id1, id2, jdc
  real,    intent(in), dimension(is:ie,js:je,4):: s2c
!--------
! Output:
!--------
  real,    intent(out), dimension(is:ie,js:je):: q_out

! Local:
  real, dimension(isd:ied,jsd:jed,1) :: q3d_in
  real, dimension(is:ie,js:je,1)     :: q3d_out

     q3d_in(:,:,1) = q_in
     call remap_xy_3d( isd, ied, jsd, jed, 1, q3d_in, &
                          is, ie, js, je, id1, id2, jdc, s2c, &
                          q3d_out )
     q_out = q3d_out(:,:,1)

  end subroutine remap_xy_2d
  subroutine     Flip1d_g5_to360(A, nx)
  integer nx
  integer i
  real A(nx), B(nx)
  B=A
  do i=1,nx/2
     A(i) = B(i+nx/2)
     A(i+nx/2) = B(i)
  enddo
  end subroutine     Flip1d_g5_to360
  subroutine     Flip2d_g5_to360(A, nx, ny)
  integer nx, ny
  integer i, j
  real A(nx, ny)
  do j=1, ny
     call  Flip1d_g5_to360(A(:,j), nx)
  enddo
  end subroutine Flip2d_g5_to360

  subroutine     Flip3d_g5_to360(A, nx, ny, nz)
  integer nx, ny, nz
  integer i, j, k
  real A(nx, ny, nz)
  do k=1, nz
     call  Flip2d_g5_to360(A(:,:,K), nx, ny)
  enddo
  end subroutine Flip3d_g5_to360
!
!
!
  SUBROUTINE WRITE_3DTILE_GFS(FNC, im, jm, nz, geolat, geolon, ps, w, zh, Qv, QL, QI, O3, uw, us, vw, vs)
  use netcdf
  IMPLICIT NONE
  include "netcdf.inc"
  character(len=*)  :: FNC
  integer :: im, jm, nz
  real    ::  geolat(im, jm), geolon(im,jm)
  real, dimension(im, jm, nz)   :: O3, QV, QL, QI, W
  real, dimension(im, jm, nz+1) :: ZH 
  real, dimension(im, jm+1, nz)   :: US, VS
  real, dimension(im+1, jm, nz)   :: Uw, Vw
  real, dimension(im, jm)         :: Ps 

   INTEGER               :: HEADER_BUFFER_VAL = 16384
!dim real, dimensions:
!        lon = 96 ;
!        lat = 96 ;
!        lonp = 97 ;
!        latp = 97 ;
!        lev = 128 ;
!        levp = 129 ;
!        ntracer = 3 ;
!variables:
!        float lon(lon) ;
!                lon:cartesian_axis = "X" ;
!        float lat(lat) ;
!                lat:cartesian_axis = "Y" ;
!        float ps(lat, lon) ;
!        float w(lev, lat, lon) ;
!        float zh(levp, lat, lon) ;
!        float sphum(lev, lat, lon) ;
!        float o3mr(lev, lat, lon) ;
!        float liq_wat(lev, lat, lon) ;
!        float u_w(lev, lat, lonp) ;
!        float v_w(lev, lat, lonp) ;
!        float u_s(lev, latp, lon) ;
!        float v_s(lev, latp, lon) ;
!
     integer :: error, ncid
     INTEGER               :: ID_DIM, ID_VAR
     INTEGER               :: N, NX, NY
     INTEGER               :: INITAL=0, FSIZE=65536
     INTEGER               :: DIM_LON, DIM_LAT, DIM_LONP, DIM_LATP
     INTEGER               :: DIM_LEV, DIM_LEVP, DIM_TRACER
     INTEGER               :: ID_LON, ID_LAT, ID_PS
     INTEGER               :: ID_W, ID_ZH, ID_SPHUM, ID_O3MR
     INTEGER               :: ID_CLWMR, ID_ICMR, ID_U_W, ID_V_W
     INTEGER               :: ID_U_S, ID_V_S
   
!     ERROR=NF90_OPEN(TRIM(FNC),NF_NOWRITE,NCID)
     INTEGER, parameter    ::   NTRACM=3


      REAL(KIND=4), ALLOCATABLE  :: CUBE_2D_4BYTE(:,:)
      REAL(KIND=4), ALLOCATABLE  :: CUBE_3D_4BYTE(:,:,:)
      ALLOCATE( CUBE_2D_4BYTE(IM,JM))
!
     ERROR = NF__CREATE(trim(FNC), IOR(NF_NETCDF4,NF_CLASSIC_MODEL),INITAL, FSIZE, NCID)
     ERROR = NF90_DEF_DIM(NCID, 'lon', IM, DIM_LON)
     ERROR = NF90_DEF_DIM(NCID, 'lat', JM, DIM_LAT)
     ERROR = NF90_DEF_DIM(NCID, 'lonp', (IM+1), DIM_LONP)
     ERROR = NF90_DEF_DIM(NCID, 'latp', (JM+1), DIM_LATP)
     ERROR = NF90_DEF_DIM(NCID, 'lev',   NZ, DIM_LEV)
     ERROR = NF90_DEF_DIM(NCID, 'levp', NZ+1, DIM_LEVP)
     ERROR = NF90_DEF_DIM(NCID, 'ntracer', NTRACM, DIM_TRACER)
     ERROR = NF90_DEF_VAR(NCID, 'lon', NF90_FLOAT, DIM_LON, ID_LON)
     ERROR = NF90_DEF_VAR(NCID, 'lat', NF90_FLOAT, DIM_LAT, ID_LAT)
     ERROR = NF_PUT_ATT_TEXT(NCID, ID_LON, "cartesian_axis", 1, "X")
     ERROR = NF_PUT_ATT_TEXT(NCID, ID_LAT, "cartesian_axis", 1, "Y")
     ERROR = NF90_DEF_VAR(NCID, 'ps', NF90_FLOAT, (/DIM_LON, DIM_LAT/), ID_PS)
     ERROR = NF90_DEF_VAR(NCID, 'w', NF90_FLOAT, (/DIM_LON, DIM_LAT, DIM_LEV/), ID_W)
     ERROR = NF90_DEF_VAR(NCID, 'zh', NF90_FLOAT, (/DIM_LON, DIM_LAT, DIM_LEVP/), ID_ZH)
     ERROR = NF90_DEF_VAR(NCID, 'sphum', NF90_FLOAT, (/DIM_LON, DIM_LAT, DIM_LEV/), ID_sphum)
     ERROR = NF90_DEF_VAR(NCID, 'o3mr', NF90_FLOAT, (/DIM_LON, DIM_LAT, DIM_LEV/), ID_o3mr)
     ERROR = NF90_DEF_VAR(NCID, 'liq_wat', NF90_FLOAT, (/DIM_LON, DIM_LAT, DIM_LEV/), ID_clwmr)
     ERROR = NF90_DEF_VAR(NCID, 'ice_wat', NF90_FLOAT, (/DIM_LON, DIM_LAT, DIM_LEV/), ID_icmr)
     ERROR = NF90_DEF_VAR(NCID, 'u_w', NF90_FLOAT, (/DIM_LONP, DIM_LAT, DIM_LEV/), ID_U_W)
     ERROR = NF90_DEF_VAR(NCID, 'v_w', NF90_FLOAT, (/DIM_LONP, DIM_LAT, DIM_LEV/), ID_V_W)
     ERROR = NF90_DEF_VAR(NCID, 'u_s', NF90_FLOAT, (/DIM_LON, DIM_LATP, DIM_LEV/), ID_U_S)
     ERROR = NF90_DEF_VAR(NCID, 'v_s', NF90_FLOAT, (/DIM_LON, DIM_LATP, DIM_LEV/), ID_V_S)
     ERROR = NF90_ENDDEF(NCID, HEADER_BUFFER_VAL, 4, 0, 4)

     CUBE_2D_4BYTE = REAL(GEOLON,4)
     ERROR = NF90_PUT_VAR(NCID, ID_LON, CUBE_2D_4BYTE(:,1))
     CUBE_2D_4BYTE = REAL(GEOLAT,4)
     ERROR = NF90_PUT_VAR(NCID, ID_LAT, CUBE_2D_4BYTE(1,:))
     CUBE_2D_4BYTE = REAL(PS,4)
     ERROR = NF90_PUT_VAR(NCID, ID_PS, CUBE_2D_4BYTE)
     DEALLOCATE(CUBE_2D_4BYTE)

      ALLOCATE(CUBE_3D_4BYTE(IM,JM,NZ))

      CUBE_3D_4BYTE = REAL(W,4)
      ERROR = NF90_PUT_VAR(NCID, ID_W, CUBE_3D_4BYTE)

      CUBE_3D_4BYTE = REAL(QV,4)
      ERROR = NF90_PUT_VAR(NCID, ID_sphum, CUBE_3D_4BYTE)
      CUBE_3D_4BYTE = REAL(QL,4)
      ERROR = NF90_PUT_VAR(NCID, ID_clwmr, CUBE_3D_4BYTE)
      CUBE_3D_4BYTE = REAL(QI,4)
      ERROR = NF90_PUT_VAR(NCID, ID_icmr, CUBE_3D_4BYTE)
      CUBE_3D_4BYTE = REAL(O3,4)
      ERROR = NF90_PUT_VAR(NCID, ID_o3mr, CUBE_3D_4BYTE)

      DEALLOCATE(CUBE_3D_4BYTE)

      ALLOCATE(CUBE_3D_4BYTE(IM,JM,nz+1))
      CUBE_3D_4BYTE = REAL(ZH,4)
      ERROR = NF90_PUT_VAR(NCID, ID_zh, CUBE_3D_4BYTE)

      DEALLOCATE(CUBE_3D_4BYTE)

      ALLOCATE(CUBE_3D_4BYTE(IM+1,JM, nz))
      CUBE_3D_4BYTE = REAL(Uw,4)
      ERROR = NF90_PUT_VAR(NCID, ID_u_w, CUBE_3D_4BYTE)
      CUBE_3D_4BYTE = REAL(Vw,4)
      ERROR = NF90_PUT_VAR(NCID, ID_v_w, CUBE_3D_4BYTE)
      DEALLOCATE(CUBE_3D_4BYTE)

      ALLOCATE(CUBE_3D_4BYTE(IM,JM+1, nz))
      CUBE_3D_4BYTE = REAL(Us,4)
      ERROR = NF90_PUT_VAR(NCID, ID_u_s, CUBE_3D_4BYTE)
      CUBE_3D_4BYTE = REAL(Vs,4)
      ERROR = NF90_PUT_VAR(NCID, ID_v_s, CUBE_3D_4BYTE)
      DEALLOCATE(CUBE_3D_4BYTE)
      ERROR = NF90_CLOSE(NCID)
   END SUBROUTINE WRITE_3DTILE_GFS
!
      subroutine compute_zh_geos5(im, jm, levp, delp, ps, zs, t, sphum, zh) 
!xxxxxxxxxxx
! VAY-2017 WAM-related computation of geopotential.... constant Grav and variable Cp
!xxxxxxxxxxx
       implicit none 
!       integer, intent(in):: ntrac
       integer, intent(in):: levp, im,jm 
!       real,    intent(in), dimension(levp+1):: ak_in, bk_in 
       real,    intent(in), dimension(im,jm):: ps, zs 
       real,    intent(in), dimension(im,jm,levp)   :: t 
       real,    intent(in), dimension(im,jm,levp)   :: sphum 
       real,    intent(in), dimension(im,jm,levp)   :: delp 
       real,    intent(out), dimension(im,jm,levp+1):: zh 
       real, dimension(im,jm):: delpk
       ! Local:                                                         
!       real, dimension(im,levp+1):: pe0, pn0 
!       real, dimension(levp+1) :: ak, bk 
       integer i,j,k 
       real, parameter :: GRAV   = 9.80665 
       real, parameter :: RDGAS  = 287.05 
       real, parameter :: RVGAS  = 461.50 
       real, parameter :: e0 = 610.71 
       real, parameter :: hlv = 2.501e6 
       real, parameter :: tfreeze = 273.15 
       real  :: zvir 
       real:: grd, rdg
       real :: pikj, dplog, psold, dpp, psdelp

       grd = grav/rdgas
       rdg = rdgas/grav
!
! tv = T* (1.+zvir*sphum(i,j,k-1))
!
       zvir = rvgas/rdgas - 1. 
!       ak = ak_in 
!       bk = bk_in 
!
! NON-sense 1.e-9   should be Zero
!
!       ak(levp+1) = max(1.e-9, ak(levp+1)) 
                                                                        
       do j = 1, jm 
       
       
!         do i=1, im 
!           pe0(i,levp+1) = ak(levp+1) 
!           pn0(i,levp+1) = log(pe0(i,levp+1)) 
!         enddo 
!                                                                        
!         do k=levp,1, -1 
!            do i=1,im 
!              pe0(i,k) = ak(k) + bk(k)*ps(i,j) 
!              pn0(i,k) = log(pe0(i,k)) 
!            enddo 
!         enddo 
!--------------------------
! Solving dZ/dp = [R/g]* T/p from  surface(1) => TL(levp+1)
! -------------------------                                                                      
         zh(1:im,j,levp+1) = zs(1:im,j)/GRAV
!
! interfaces 
!
   
         do i = 1, im 
            pikj = ps(i,j)
	    psdelp = sum(delp(i,j,1:levp))
!            print *, pikj/psdelp, ps(i,j), psdelp
          pikj = psdelp 	        
         do k = levp, 1, -1 
	   psold = pikj	  
           dpp = delp(i,j,k) 
!           if (k==1)  dpp = dpp*.99
	   pikj = pikj-dpp
           
	   if (pikj .le.0.) pikj=0.0 
	    dplog = 2.*delp(i,j,k)/(pikj+psold) 
            zh(i,j,k) = zh(i,j,k+1)+t(i,j,k)*(1.+zvir*sphum(i,j,k))*rdg*dplog                                 
           enddo 
!            zh(i,j,1)=zh(i,j,2)
         enddo 
                                                                        
       enddo 
                                                                        
      end subroutine compute_zh_geos5 
!
      subroutine remap_zdim(im, jm, nzi, A, Pa, nzg, gPa, gA)
      implicit NONE
      integer :: im, jm, nzi, nzg
      real    ::  A(im, jm, nzi), Pa(im, jm, nzi)
      real    :: gA(im, jm, nzg), gPa(im, jm, nzg)
      integer :: i, j, k
      real    :: B(nzi),  gB(nzg)
      real    :: zB(nzi), zgB(nzg)

!      print *, im, jm, nzi, nzg, ' remap_zdim '
      do i=1, im
       do j=1, jm
       do k=1, nzi
       B(k) = A(i,j,k)
         zB(k) = -7.*alog(pA (i,j,k)*1.e-5)
       enddo
       do k=1, nzg
        zgB(k) = -7.*alog(gpA(i,j,k)*1.e-5)
      enddo
          call zmap_1d(nzi, Zb, B, nzg, Zgb, gB)
           do k=1, nzg
            gA(i,j,k) = gB(k)
             enddo
       enddo
      enddo
      end subroutine remap_zdim

      subroutine zmap_1d(na, Za, A, nb, Zb, B)
      implicit none 
      integer :: na, nb
      real, dimension(na) :: Za, A
      real, dimension(nb) :: Zb, B
      real :: w1, w2, X, plim1, plim2
      integer :: i,j, k, ka, ks
         plim1 = maxval(Za)
         plim2 = minval(Za)
          ks =na-1
 
        do k=nb, 1, -1
           X = Zb(k)
           if (x >= plim1 ) B(k) = A(1)  ! 1   top lid
           if (x <= plim2 ) B(k) = A(na) ! na - bottom
           if (X < plim1 .and. X > plim2) then
             do ka=ks, 1, -1
               if (X .ge. Za(ka+1) .and. X .le. Za(ka) ) exit
             enddo
               ks = ka
!       
               B(k) = A(ks)+(A(ks)-A(ks+1))/(Za(ks)-Za(ks+1))*(X-Za(ks))    
           endif
        enddo

        RETURN
!          print *, na, nb , 'zmap_1d  na, nb '
          print *, plim2, plim1, ' zlim2, zlim1 ' 
          print *, minval(zb), maxval(zb), ' zBlim2, zBlim1 '
115     format(i4, 4(2x, F12.5))
        do k=1,72
        write(6,115) k, B(k), zb(k), za(k), A(k)
        enddo
        print *
        do k=1,15
        write(6,115) nb+1-k, B(nb+1-k), zb(nb+1-k), za(na+1-k), A(na+1-k)
        enddo
!          print *, na, nb , ' return zmap_1d  na, nb '
        stop
      end subroutine zmap_1d

