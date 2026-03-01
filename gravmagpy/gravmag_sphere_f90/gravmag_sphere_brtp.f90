program gravmag_sphere_brtp
  use, intrinsic :: iso_fortran_env, only: real32, int32
  use gravmag_sphere_state
  use gravmag_sphere_subs
  use gravmag_sphere_physics
  implicit none

  real(real32) :: rsphere_km
  character(len=256) :: infile, outfile
  character(len=256) :: title
  integer(int32) :: ios, knk

  ! Card 2
  real(real32) :: lat0_deg, lon0_deg, dlat_deg, dlon_deg, elvo_km
  integer(int32) :: nlat, nlon

  ! Card 3
  integer(int32) :: nr, ntheta, nphi, nblim

  ! Card 4 (kept for compatibility)
  real(real32) :: htheta, hphi, phi1, phi2, rho_dummy

  ! Card 5
  integer(int32) :: ifield, nrem
  real(real32) :: M_amp_Apm, inc_deg, dec_deg

  ! Card 6
  integer(int32) :: iprint, nfile, nopt
  real(real32) :: first, conint

  ! Card 7 (NEW format for both nblim=0 and nblim=1)
  integer(int32) :: npts
  real(real32) :: depth_top_km, depth_bot_km
  real(real32), allocatable :: poly_lat(:), poly_lon(:)
  real(real32) :: lat_lo, lat_hi, lon_lo, lon_hi, lon_ref

  ! constants
  real(real32), parameter :: pi = 3.14159265358979323846_real32
  real(real32), parameter :: deg2rad = pi/180.0_real32
  real(real32), parameter :: T_to_nT = 1.0e9_real32

  ! obs grid xyz (m)
  real(real32), allocatable :: xo(:,:), yo(:,:), zo(:,:)

  ! magnetization vector (A/m), planet-fixed Cartesian
  real(real32) :: cosI, sinI, cosD, sinD
  real(real32) :: Mx, My, Mz

  ! quadrature arrays (nodes/weights)
  real(real32), allocatable :: xr(:), wr(:), xlat(:), wlat(:), xlon(:), wlon(:)

  ! radii (m)
  real(real32) :: r_top_m, r_bot_m, rsphere_m, ro_m

  ! loops
  integer(int32) :: i, j, ir, it, ip
  real(real32) :: lat_q_deg, lon_q_deg, lat_q_rad, lon_q_rad
  real(real32) :: r_q_m

  ! scalar weights (avoid name collisions)
  real(real32) :: w_r_s, w_lat_s, w_lon_s, dV_m3

  real(real32) :: xs, ys, zs
  real(real32) :: dBx, dBy, dBz
  real(real32) :: lon_out, lat_out, btot_nt

  if (command_argument_count() < 3) then
    write(*,*) "Usage: gravmag_sphere_brtp <R_sphere_km> <input.in> <output.txt>"
    stop 2
  end if

  call get_command_argument(1, title)
  read(title,*,iostat=ios) rsphere_km
  if (ios /= 0) stop "ERROR: could not parse R_sphere_km"

  call get_command_argument(2, infile)
  call get_command_argument(3, outfile)

  open(4, file=trim(infile), status='old', iostat=ios)
  if (ios /= 0) stop "ERROR: cannot open input file"

  open(2, file=trim(outfile), status='replace', iostat=ios)
  if (ios /= 0) stop "ERROR: cannot open output file"

  open(6, file='sphout', status='replace', iostat=ios)
  if (ios /= 0) stop "ERROR: cannot open sphout"

  rsphere_m = rsphere_km * 1000.0_real32
  knk = 0

  do
    ! ---- read title (skip blanks/comments) ----
    title = ''
    do
      read(4,'(A)', iostat=ios) title
      if (ios /= 0) exit
      if (len_trim(title) == 0) cycle
      if (title(1:1) == '!' .or. title(1:1) == '#') cycle
      exit
    end do
    if (ios /= 0) exit

    knk = knk + 1
    write(6,'(A,I0)') '--- BODY ', knk
    write(6,'(A)') trim(title)

    ! ---- Cards 2-6 ----
    read(4,*,iostat=ios) lat0_deg, lon0_deg, dlat_deg, dlon_deg, elvo_km, nlat, nlon
    if (ios /= 0) stop "ERROR reading Card 2"

    read(4,*,iostat=ios) nr, ntheta, nphi, nblim
    if (ios /= 0) stop "ERROR reading Card 3"

    read(4,*,iostat=ios) htheta, hphi, phi1, phi2, rho_dummy
    if (ios /= 0) stop "ERROR reading Card 4"

    read(4,*,iostat=ios) ifield, nrem, M_amp_Apm, inc_deg, dec_deg
    if (ios /= 0) stop "ERROR reading Card 5"
    if (nrem /= 1) stop "ERROR: magnetics require nrem=1"

    read(4,*,iostat=ios) iprint, nfile, nopt, first, conint
    if (ios /= 0) stop "ERROR reading Card 6"

    ! ---- magnetization vector (planet-fixed Cartesian) ----
    cosI = cos(inc_deg*deg2rad)
    sinI = sin(inc_deg*deg2rad)
    cosD = cos(dec_deg*deg2rad)
    sinD = sin(dec_deg*deg2rad)

    Mx = M_amp_Apm * (cosI*cosD)
    My = M_amp_Apm * (cosI*sinD)
    Mz = M_amp_Apm * (sinI)

    ! ---- Card 7: NEW polygon header + vertices (Option B) ----
    if (allocated(poly_lat)) deallocate(poly_lat, poly_lon)

    read(4,*,iostat=ios) npts, depth_top_km, depth_bot_km
    if (ios /= 0) stop "ERROR reading Card 7 polygon header: npts depth_top depth_bot"
    if (npts < 3) stop "ERROR: polygon needs >=3 vertices"

    allocate(poly_lat(npts), poly_lon(npts))
    do i=1,npts
      read(4,*,iostat=ios) poly_lat(i), poly_lon(i)   ! lat lon
      if (ios /= 0) stop "ERROR reading polygon vertex (lat lon)"
      call wrap180_scalar(poly_lon(i))
    end do

    ! ---- unwrap polygon lons to continuous branch ----
    lon_ref = sum(poly_lon) / real(npts,real32)
    call wrap180_scalar(lon_ref)
    call unwrap_lons_inplace(poly_lon, lon_ref)

    lat_lo = minval(poly_lat); lat_hi = maxval(poly_lat)
    lon_lo = minval(poly_lon); lon_hi = maxval(poly_lon)

    ! ---- depths -> radii ----
    if (depth_top_km < 0.0_real32 .or. depth_bot_km < 0.0_real32) stop "ERROR: depth must be >=0"
    if (depth_bot_km <= depth_top_km) stop "ERROR: need depth_bot > depth_top"

    r_top_m = rsphere_m - depth_top_km*1000.0_real32
    r_bot_m = rsphere_m - depth_bot_km*1000.0_real32
    if (r_top_m <= 0.0_real32 .or. r_bot_m <= 0.0_real32) stop "ERROR: depth too large"

    ! enforce r_top_m > r_bot_m
    if (r_bot_m > r_top_m) then
      call swap_real(r_bot_m, r_top_m)
    end if

    ro_m = (rsphere_km + elvo_km) * 1000.0_real32

    ! ---- alloc / build obs grid ----
    call ensure_alloc(nlat, nlon)
    call clear_fields()

    if (.not. allocated(xo)) then
      allocate(xo(nlat,nlon), yo(nlat,nlon), zo(nlat,nlon))
    else
      if (size(xo,1) /= nlat .or. size(xo,2) /= nlon) then
        deallocate(xo,yo,zo)
        allocate(xo(nlat,nlon), yo(nlat,nlon), zo(nlat,nlon))
      end if
    end if

    do i=1,nlat
      do j=1,nlon
        lat_out = lat0_deg + real(i-1,real32)*dlat_deg
        lon_out = lon0_deg + real(j-1,real32)*dlon_deg
        call wrap180_scalar(lon_out)

        xo(i,j) = ro_m * cos(lat_out*deg2rad) * cos(lon_out*deg2rad)
        yo(i,j) = ro_m * cos(lat_out*deg2rad) * sin(lon_out*deg2rad)
        zo(i,j) = ro_m * sin(lat_out*deg2rad)
      end do
    end do

    ! ---- quadrature on [-1,1] ----
    if (allocated(xr)) deallocate(xr,wr,xlat,wlat,xlon,wlon)
    allocate(xr(nr), wr(nr), xlat(ntheta), wlat(ntheta), xlon(nphi), wlon(nphi))

    call gauss_legendre(nr,     xr,   wr)
    call gauss_legendre(ntheta, xlat, wlat)
    call gauss_legendre(nphi,   xlon, wlon)

    ! ---- integrate over bbox, masked by polygon ----
    do ip=1,nphi
      lon_q_deg = 0.5_real32*((lon_hi-lon_lo)*xlon(ip) + (lon_hi+lon_lo))
      w_lon_s   = 0.5_real32*(lon_hi-lon_lo) * wlon(ip)   ! degrees

      do it=1,ntheta
        lat_q_deg = 0.5_real32*((lat_hi-lat_lo)*xlat(it) + (lat_hi+lat_lo))
        w_lat_s   = 0.5_real32*(lat_hi-lat_lo) * wlat(it) ! degrees

        if (.not. point_in_poly_lonlat(lon_q_deg, lat_q_deg, poly_lon, poly_lat)) cycle

        lat_q_rad = lat_q_deg*deg2rad
        lon_q_rad = lon_q_deg*deg2rad

        do ir=1,nr
          r_q_m = 0.5_real32*((r_top_m-r_bot_m)*xr(ir) + (r_top_m+r_bot_m))
          w_r_s = 0.5_real32*(r_top_m-r_bot_m) * wr(ir)  ! meters

          xs = r_q_m * cos(lat_q_rad) * cos(lon_q_rad)
          ys = r_q_m * cos(lat_q_rad) * sin(lon_q_rad)
          zs = r_q_m * sin(lat_q_rad)

          dV_m3 = (r_q_m*r_q_m) * cos(lat_q_rad) * w_r_s * (w_lat_s*deg2rad) * (w_lon_s*deg2rad)
          if (dV_m3 <= 0.0_real32) cycle

          do i=1,nlat
            do j=1,nlon
              call add_dipole_field_cart(xo(i,j), yo(i,j), zo(i,j), xs, ys, zs, &
                                         Mx*dV_m3, My*dV_m3, Mz*dV_m3, dBx, dBy, dBz)
              bx(i,j) = bx(i,j) + dBx
              by(i,j) = by(i,j) + dBy
              bz(i,j) = bz(i,j) + dBz
            end do
          end do

        end do
      end do
    end do

    ! ---- output ----
    if (nfile /= 0) then
      write(2,'(A)') '# body_id lon_deg[-180,180) lat_deg Bx_nT By_nT Bz_nT Btot_nT'
      do i=1,nlat
        do j=1,nlon
          lat_out = lat0_deg + real(i-1,real32)*dlat_deg
          lon_out = lon0_deg + real(j-1,real32)*dlon_deg
          call wrap180_scalar(lon_out)

          btot_nt = sqrt( (bx(i,j)*T_to_nT)**2 + (by(i,j)*T_to_nT)**2 + (bz(i,j)*T_to_nT)**2 )

          write(2,'(I6,1X,F12.6,1X,F12.6,1X,ES15.7,1X,ES15.7,1X,ES15.7,1X,ES15.7)') &
               knk, lon_out, lat_out, bx(i,j)*T_to_nT, by(i,j)*T_to_nT, bz(i,j)*T_to_nT, btot_nt
        end do
      end do
    end if

  end do

  close(4); close(2); close(6)

contains

  subroutine swap_real(a,b)
    real(real32), intent(inout) :: a,b
    real(real32) :: t
    t=a; a=b; b=t
  end subroutine swap_real

end program gravmag_sphere_brtp