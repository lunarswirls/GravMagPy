program gravmag_sphere_dipole_grid_fit
  use, intrinsic :: iso_fortran_env, only: int32, int64, real64
  implicit none

  integer, parameter :: wp = real64
  real(wp), parameter :: pi = 3.1415926535897932384626433832795_wp
  real(wp), parameter :: deg2rad = pi / 180.0_wp
  real(wp), parameter :: mu0_over_4pi = 1.0e-7_wp
  real(wp), parameter :: T_to_nT = 1.0e9_wp

  character(len=1024) :: csv_list_arg, pred_csv, dipole_csv, tmp
  integer(int32) :: argc, ios

  real(wp) :: rsphere_km, rsphere_m
  real(wp) :: source_depth_km, source_dlat_deg, source_dlon_deg
  real(wp) :: reg_lambda, source_dr_km
  real(wp) :: lat_pad_deg, lon_pad_deg, max_memory_mib
  real(wp) :: source_depth_km_use, source_dlat_deg_use, source_dlon_deg_use
  real(wp) :: source_dr_km_use, lat_pad_deg_use, lon_pad_deg_use
  integer(int32) :: source_nr

  integer(int32) :: nobs, obs_cap
  real(wp), allocatable :: obs_lon_deg(:), obs_lat_deg(:), obs_radius_m(:)
  real(wp), allocatable :: obs_bx_t(:), obs_by_t(:), obs_bz_t(:)
  real(wp), allocatable :: obs_x(:), obs_y(:), obs_z(:)
  real(wp), allocatable :: obs_lon_unwrapped(:)

  real(wp) :: lon_ref_deg, lat_min, lat_max, lon_min, lon_max
  real(wp) :: min_obs_radius_m, max_obs_radius_m
  real(wp) :: auto_spacing_deg
  integer(int32) :: nlat_src, nlon_src, nsrc
  integer(int32) :: ilat, ilon, ir, isrc
  real(wp), allocatable :: src_lon_deg(:), src_lat_deg(:), src_radius_m(:)
  real(wp), allocatable :: src_x(:), src_y(:), src_z(:)

  integer(int64) :: nrow64, ncol64
  integer(int32) :: nrow, ncol
  integer(int64) :: mem_g_bytes, mem_ata_bytes, mem_total_bytes
  real(wp), allocatable :: G(:,:), ATA(:,:), ATb(:), model(:), y(:), ypred(:)
  real(wp) :: diag_scale, diag_sum
  integer(int32) :: iobs, jsrc, r0, c0
  real(wp) :: K(3,3)
  integer(int32) :: info

  real(wp), allocatable :: pred_bx_t(:), pred_by_t(:), pred_bz_t(:)
  real(wp), allocatable :: res_bx_t(:), res_by_t(:), res_bz_t(:)
  real(wp) :: rmse_bx_nt, rmse_by_nt, rmse_bz_nt, rmse_btot_nt

  argc = command_argument_count()
  if (argc < 4) then
    write(*,'(A)') 'Usage: gravmag_sphere_dipole_grid_fit <R_sphere_km> <obs_csvs> <predictions.csv> <dipoles.csv> '// &
                   '[source_depth_km] [source_dlat_deg] [source_dlon_deg] [reg_lambda] [source_nr] '// &
                   '[source_dr_km] [lat_pad_deg] [lon_pad_deg] [max_memory_mib]'
    write(*,'(A)') '  obs_csvs: comma-separated CSV paths. Supported columns include lon/lat + '// &
                   'altitude_m|altitude_km|radius_m|radius_km + Bx/By/Bz in nT or Tesla.'
    stop 2
  end if

  call get_command_argument(1, tmp)
  read(tmp,*,iostat=ios) rsphere_km
  if (ios /= 0 .or. rsphere_km <= 0.0_wp) stop 'ERROR: invalid R_sphere_km'
  rsphere_m = rsphere_km * 1000.0_wp

  call get_command_argument(2, csv_list_arg)
  call get_command_argument(3, pred_csv)
  call get_command_argument(4, dipole_csv)

  source_depth_km = 0.0_wp
  source_dlat_deg = 0.0_wp
  source_dlon_deg = 0.0_wp
  reg_lambda = 1.0e-4_wp
  source_nr = 1_int32
  source_dr_km = 0.0_wp
  lat_pad_deg = -1.0_wp
  lon_pad_deg = -1.0_wp
  max_memory_mib = 768.0_wp

  if (argc >= 5) then
    call get_command_argument(5, tmp)
    read(tmp,*,iostat=ios) source_depth_km
    if (ios /= 0) stop 'ERROR: invalid source_depth_km'
  end if
  if (argc >= 6) then
    call get_command_argument(6, tmp)
    read(tmp,*,iostat=ios) source_dlat_deg
    if (ios /= 0) stop 'ERROR: invalid source_dlat_deg'
  end if
  if (argc >= 7) then
    call get_command_argument(7, tmp)
    read(tmp,*,iostat=ios) source_dlon_deg
    if (ios /= 0) stop 'ERROR: invalid source_dlon_deg'
  end if
  if (argc >= 8) then
    call get_command_argument(8, tmp)
    read(tmp,*,iostat=ios) reg_lambda
    if (ios /= 0 .or. reg_lambda < 0.0_wp) stop 'ERROR: invalid reg_lambda'
  end if
  if (argc >= 9) then
    call get_command_argument(9, tmp)
    read(tmp,*,iostat=ios) source_nr
    if (ios /= 0 .or. source_nr < 1) stop 'ERROR: invalid source_nr'
  end if
  if (argc >= 10) then
    call get_command_argument(10, tmp)
    read(tmp,*,iostat=ios) source_dr_km
    if (ios /= 0) stop 'ERROR: invalid source_dr_km'
  end if
  if (argc >= 11) then
    call get_command_argument(11, tmp)
    read(tmp,*,iostat=ios) lat_pad_deg
    if (ios /= 0) stop 'ERROR: invalid lat_pad_deg'
  end if
  if (argc >= 12) then
    call get_command_argument(12, tmp)
    read(tmp,*,iostat=ios) lon_pad_deg
    if (ios /= 0) stop 'ERROR: invalid lon_pad_deg'
  end if
  if (argc >= 13) then
    call get_command_argument(13, tmp)
    read(tmp,*,iostat=ios) max_memory_mib
    if (ios /= 0 .or. max_memory_mib <= 0.0_wp) stop 'ERROR: invalid max_memory_mib'
  end if

  obs_cap = 1024
  nobs = 0
  allocate(obs_lon_deg(obs_cap), obs_lat_deg(obs_cap), obs_radius_m(obs_cap))
  allocate(obs_bx_t(obs_cap), obs_by_t(obs_cap), obs_bz_t(obs_cap))

  call read_observation_csv_list(rsphere_m, trim(csv_list_arg), obs_cap, nobs, &
                                 obs_lon_deg, obs_lat_deg, obs_radius_m, &
                                 obs_bx_t, obs_by_t, obs_bz_t)
  if (nobs <= 0) stop 'ERROR: no observations were read from input CSVs'

  if (nobs < obs_cap) then
    call shrink_observation_arrays(nobs, obs_lon_deg, obs_lat_deg, obs_radius_m, obs_bx_t, obs_by_t, obs_bz_t)
  end if

  allocate(obs_x(nobs), obs_y(nobs), obs_z(nobs), obs_lon_unwrapped(nobs))
  lon_ref_deg = circular_mean_deg(obs_lon_deg)
  do iobs = 1, nobs
    obs_lon_unwrapped(iobs) = unwrap_lon_to_ref(obs_lon_deg(iobs), lon_ref_deg)
    call sph_to_cart(obs_radius_m(iobs), obs_lat_deg(iobs), obs_lon_deg(iobs), obs_x(iobs), obs_y(iobs), obs_z(iobs))
  end do

  min_obs_radius_m = minval(obs_radius_m)
  max_obs_radius_m = maxval(obs_radius_m)
  lat_min = minval(obs_lat_deg)
  lat_max = maxval(obs_lat_deg)
  lon_min = minval(obs_lon_unwrapped)
  lon_max = maxval(obs_lon_unwrapped)

  source_dlat_deg_use = source_dlat_deg
  if (source_dlat_deg_use <= 0.0_wp) source_dlat_deg_use = median_positive_spacing(obs_lat_deg)
  if (source_dlat_deg_use <= 0.0_wp) source_dlat_deg_use = 1.0_wp

  source_dlon_deg_use = source_dlon_deg
  if (source_dlon_deg_use <= 0.0_wp) source_dlon_deg_use = median_positive_spacing(obs_lon_unwrapped)
  if (source_dlon_deg_use <= 0.0_wp) source_dlon_deg_use = source_dlat_deg_use

  auto_spacing_deg = max(source_dlat_deg_use, source_dlon_deg_use)

  source_depth_km_use = source_depth_km
  if (source_depth_km_use <= 0.0_wp) then
    source_depth_km_use = max(5.0_wp, 2.0_wp * auto_spacing_deg * deg2rad * rsphere_km)
  end if

  source_dr_km_use = source_dr_km
  if (source_nr > 1 .and. source_dr_km_use <= 0.0_wp) then
    source_dr_km_use = max(5.0_wp, 0.5_wp * source_depth_km_use)
  end if
  if (source_nr == 1) source_dr_km_use = 0.0_wp

  lat_pad_deg_use = lat_pad_deg
  if (lat_pad_deg_use < 0.0_wp) lat_pad_deg_use = source_dlat_deg_use
  lon_pad_deg_use = lon_pad_deg
  if (lon_pad_deg_use < 0.0_wp) lon_pad_deg_use = source_dlon_deg_use

  lat_min = max(-89.75_wp, lat_min - lat_pad_deg_use)
  lat_max = min(89.75_wp, lat_max + lat_pad_deg_use)
  lon_min = lon_min - lon_pad_deg_use
  lon_max = lon_max + lon_pad_deg_use

  nlat_src = max(1_int32, ceiling((lat_max - lat_min) / source_dlat_deg_use) + 1_int32)
  nlon_src = max(1_int32, ceiling((lon_max - lon_min) / source_dlon_deg_use) + 1_int32)
  nsrc = nlat_src * nlon_src * source_nr
  if (nsrc <= 0) stop 'ERROR: invalid source grid size'

  allocate(src_lon_deg(nsrc), src_lat_deg(nsrc), src_radius_m(nsrc))
  allocate(src_x(nsrc), src_y(nsrc), src_z(nsrc))

  isrc = 0
  do ir = 1, source_nr
    do ilat = 1, nlat_src
      do ilon = 1, nlon_src
        isrc = isrc + 1
        src_lat_deg(isrc) = lat_min + real(ilat - 1, wp) * source_dlat_deg_use
        src_lon_deg(isrc) = lon_min + real(ilon - 1, wp) * source_dlon_deg_use
        src_radius_m(isrc) = min_obs_radius_m - source_depth_km_use * 1000.0_wp - &
                             real(ir - 1, wp) * source_dr_km_use * 1000.0_wp
        if (src_radius_m(isrc) <= 0.0_wp) stop 'ERROR: source radius became non-positive; reduce source depth'
        if (src_radius_m(isrc) >= min_obs_radius_m) stop 'ERROR: source grid must be below minimum observation radius'
        call sph_to_cart(src_radius_m(isrc), src_lat_deg(isrc), wrap180(src_lon_deg(isrc)), &
                         src_x(isrc), src_y(isrc), src_z(isrc))
      end do
    end do
  end do

  nrow64 = 3_int64 * int(nobs, int64)
  ncol64 = 3_int64 * int(nsrc, int64)
  if (nrow64 > huge(1_int32) .or. ncol64 > huge(1_int32)) then
    stop 'ERROR: system too large for current integer indexing'
  end if
  nrow = int(nrow64, int32)
  ncol = int(ncol64, int32)

  mem_g_bytes = nrow64 * ncol64 * 8_int64
  mem_ata_bytes = ncol64 * ncol64 * 8_int64
  mem_total_bytes = mem_g_bytes + mem_ata_bytes + (4_int64 * ncol64 + 2_int64 * nrow64) * 8_int64
  if (bytes_to_mib(mem_total_bytes) > max_memory_mib) then
    write(*,'(A,F12.3,A)') 'ERROR: estimated system memory ', bytes_to_mib(mem_total_bytes), &
                           ' MiB exceeds max_memory_mib. Increase spacing or raise the limit.'
    stop 2
  end if

  allocate(G(nrow, ncol), ATA(ncol, ncol), ATb(ncol), model(ncol), y(nrow), ypred(nrow))
  G(:,:) = 0.0_wp
  y(:) = 0.0_wp

  do iobs = 1, nobs
    r0 = 3 * (iobs - 1)
    y(r0 + 1) = obs_bx_t(iobs)
    y(r0 + 2) = obs_by_t(iobs)
    y(r0 + 3) = obs_bz_t(iobs)
    do jsrc = 1, nsrc
      call dipole_kernel_block(obs_x(iobs), obs_y(iobs), obs_z(iobs), &
                               src_x(jsrc), src_y(jsrc), src_z(jsrc), K)
      c0 = 3 * (jsrc - 1)
      G(r0 + 1:r0 + 3, c0 + 1:c0 + 3) = K(:,:)
    end do
  end do

  ATA(:,:) = matmul(transpose(G), G)
  ATb(:) = matmul(transpose(G), y)
  diag_sum = 0.0_wp
  do isrc = 1, ncol
    diag_sum = diag_sum + ATA(isrc, isrc)
  end do
  diag_scale = max(diag_sum / real(ncol, wp), 1.0e-30_wp)
  do isrc = 1, ncol
    ATA(isrc, isrc) = ATA(isrc, isrc) + reg_lambda * diag_scale + 1.0e-12_wp * diag_scale
  end do

  call solve_linear_system(ATA, ATb, model, ncol, info)
  if (info /= 0) stop 'ERROR: dipole-grid normal equations are singular or ill-conditioned'

  ypred(:) = matmul(G, model)
  allocate(pred_bx_t(nobs), pred_by_t(nobs), pred_bz_t(nobs))
  allocate(res_bx_t(nobs), res_by_t(nobs), res_bz_t(nobs))
  do iobs = 1, nobs
    r0 = 3 * (iobs - 1)
    pred_bx_t(iobs) = ypred(r0 + 1)
    pred_by_t(iobs) = ypred(r0 + 2)
    pred_bz_t(iobs) = ypred(r0 + 3)
    res_bx_t(iobs) = pred_bx_t(iobs) - obs_bx_t(iobs)
    res_by_t(iobs) = pred_by_t(iobs) - obs_by_t(iobs)
    res_bz_t(iobs) = pred_bz_t(iobs) - obs_bz_t(iobs)
  end do

  rmse_bx_nt = T_to_nT * sqrt(sum(res_bx_t * res_bx_t) / real(max(1, nobs), wp))
  rmse_by_nt = T_to_nT * sqrt(sum(res_by_t * res_by_t) / real(max(1, nobs), wp))
  rmse_bz_nt = T_to_nT * sqrt(sum(res_bz_t * res_bz_t) / real(max(1, nobs), wp))
  rmse_btot_nt = sqrt(sum((btot_vector(pred_bx_t, pred_by_t, pred_bz_t) - &
                            btot_vector(obs_bx_t, obs_by_t, obs_bz_t))**2) / real(max(1, nobs), wp)) * T_to_nT

  call write_prediction_csv(trim(pred_csv), rsphere_m, nobs, obs_lon_deg, obs_lat_deg, obs_radius_m, &
                            obs_bx_t, obs_by_t, obs_bz_t, pred_bx_t, pred_by_t, pred_bz_t)
  call write_dipole_csv(trim(dipole_csv), rsphere_m, min_obs_radius_m, nsrc, nlat_src, nlon_src, source_nr, &
                        src_lon_deg, src_lat_deg, src_radius_m, model)

  write(*,'(A,I0)') 'Observation count: ', nobs
  write(*,'(A,I0,A,I0,A,I0,A,I0)') 'Source grid: nlat=', nlat_src, ' nlon=', nlon_src, ' nr=', source_nr, ' total=', nsrc
  write(*,'(A,F8.3,A,F8.3,A,F8.3)') 'Source spacing: dlat=', source_dlat_deg_use, ' deg dlon=', source_dlon_deg_use, &
                                     ' deg depth=', source_depth_km_use
  if (source_nr > 1) write(*,'(A,F8.3,A)') 'Radial layer spacing: ', source_dr_km_use, ' km'
  write(*,'(A,F12.3)') 'Estimated memory (MiB): ', bytes_to_mib(mem_total_bytes)
  write(*,'(A,ES12.4)') 'Regularization lambda: ', reg_lambda
  write(*,'(A,F12.6)') 'RMSE Bx (nT): ', rmse_bx_nt
  write(*,'(A,F12.6)') 'RMSE By (nT): ', rmse_by_nt
  write(*,'(A,F12.6)') 'RMSE Bz (nT): ', rmse_bz_nt
  write(*,'(A,F12.6)') 'RMSE Btot (nT): ', rmse_btot_nt
  write(*,'(A)') 'Done.'

contains

  subroutine read_observation_csv_list(rsphere_m_in, csv_list, cap, nobs_out, lon_deg, lat_deg, radius_m, bx_t, by_t, bz_t)
    real(wp), intent(in) :: rsphere_m_in
    character(len=*), intent(in) :: csv_list
    integer(int32), intent(inout) :: cap, nobs_out
    real(wp), allocatable, intent(inout) :: lon_deg(:), lat_deg(:), radius_m(:), bx_t(:), by_t(:), bz_t(:)

    character(len=1024), allocatable :: paths(:)
    integer(int32) :: npaths, ip

    call split_path_list(csv_list, paths, npaths)
    if (npaths <= 0) stop 'ERROR: obs_csvs list is empty'
    do ip = 1, npaths
      call read_observation_csv_file(rsphere_m_in, trim(paths(ip)), cap, nobs_out, lon_deg, lat_deg, radius_m, bx_t, by_t, bz_t)
    end do
  end subroutine read_observation_csv_list

  subroutine read_observation_csv_file(rsphere_m_in, path, cap, nobs_io, lon_deg, lat_deg, radius_m, bx_t, by_t, bz_t)
    real(wp), intent(in) :: rsphere_m_in
    character(len=*), intent(in) :: path
    integer(int32), intent(inout) :: cap, nobs_io
    real(wp), allocatable, intent(inout) :: lon_deg(:), lat_deg(:), radius_m(:), bx_t(:), by_t(:), bz_t(:)

    integer(int32) :: unit, ios_local, nfields
    integer(int32) :: idx_lon, idx_lat, idx_alt_m, idx_alt_km, idx_rad_m, idx_rad_km
    integer(int32) :: idx_bx_nt, idx_by_nt, idx_bz_nt, idx_bx_t, idx_by_t, idx_bz_t
    character(len=4096) :: line
    character(len=256), allocatable :: fields(:)
    logical :: have_header
    real(wp) :: lonv, latv, radv, bxv_t, byv_t, bzv_t

    open(newunit=unit, file=trim(path), status='old', action='read', iostat=ios_local)
    if (ios_local /= 0) stop 'ERROR: cannot open observation CSV file'

    have_header = .false.
    do
      read(unit,'(A)',iostat=ios_local) line
      if (ios_local /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      call split_csv_line(trim(line), fields, nfields)
      call locate_csv_columns(fields, nfields, idx_lon, idx_lat, idx_alt_m, idx_alt_km, idx_rad_m, idx_rad_km, &
                              idx_bx_nt, idx_by_nt, idx_bz_nt, idx_bx_t, idx_by_t, idx_bz_t)
      have_header = .true.
      exit
    end do
    if (.not. have_header) stop 'ERROR: CSV file does not contain a readable header'

    if (idx_lon <= 0 .or. idx_lat <= 0) stop 'ERROR: CSV requires lon/lat columns'
    if (max(idx_alt_m, idx_alt_km, idx_rad_m, idx_rad_km) <= 0) then
      stop 'ERROR: CSV requires altitude_m|altitude_km|radius_m|radius_km column'
    end if
    if ((idx_bx_nt <= 0 .and. idx_bx_t <= 0) .or. (idx_by_nt <= 0 .and. idx_by_t <= 0) .or. &
        (idx_bz_nt <= 0 .and. idx_bz_t <= 0)) then
      stop 'ERROR: CSV requires Bx/By/Bz columns'
    end if

    do
      read(unit,'(A)',iostat=ios_local) line
      if (ios_local /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle

      call split_csv_line(trim(line), fields, nfields)
      call parse_observation_fields( &
        rsphere_m_in, path, fields, nfields, &
        idx_lon, idx_lat, idx_alt_m, idx_alt_km, idx_rad_m, idx_rad_km, &
        idx_bx_nt, idx_by_nt, idx_bz_nt, idx_bx_t, idx_by_t, idx_bz_t, &
        lonv, latv, radv, bxv_t, byv_t, bzv_t &
      )

      nobs_io = nobs_io + 1
      if (nobs_io > cap) call grow_observation_arrays(cap, lon_deg, lat_deg, radius_m, bx_t, by_t, bz_t)
      lon_deg(nobs_io) = wrap180(lonv)
      lat_deg(nobs_io) = latv
      radius_m(nobs_io) = radv
      bx_t(nobs_io) = bxv_t
      by_t(nobs_io) = byv_t
      bz_t(nobs_io) = bzv_t
    end do

    close(unit)
  end subroutine read_observation_csv_file

  subroutine locate_csv_columns(fields, nfields, idx_lon, idx_lat, idx_alt_m, idx_alt_km, idx_rad_m, idx_rad_km, &
                                idx_bx_nt, idx_by_nt, idx_bz_nt, idx_bx_t, idx_by_t, idx_bz_t)
    character(len=*), intent(in) :: fields(:)
    integer(int32), intent(in) :: nfields
    integer(int32), intent(out) :: idx_lon, idx_lat, idx_alt_m, idx_alt_km, idx_rad_m, idx_rad_km
    integer(int32), intent(out) :: idx_bx_nt, idx_by_nt, idx_bz_nt, idx_bx_t, idx_by_t, idx_bz_t

    integer(int32) :: i
    character(len=256) :: key
    integer(int32) :: idx_bx_plain_nt, idx_by_plain_nt, idx_bz_plain_nt
    integer(int32) :: idx_bx_sel_nt, idx_by_sel_nt, idx_bz_sel_nt
    integer(int32) :: idx_bx_sse_nt, idx_by_sse_nt, idx_bz_sse_nt

    idx_lon = 0; idx_lat = 0; idx_alt_m = 0; idx_alt_km = 0; idx_rad_m = 0; idx_rad_km = 0
    idx_bx_nt = 0; idx_by_nt = 0; idx_bz_nt = 0; idx_bx_t = 0; idx_by_t = 0; idx_bz_t = 0
    idx_bx_plain_nt = 0; idx_by_plain_nt = 0; idx_bz_plain_nt = 0
    idx_bx_sel_nt = 0; idx_by_sel_nt = 0; idx_bz_sel_nt = 0
    idx_bx_sse_nt = 0; idx_by_sse_nt = 0; idx_bz_sse_nt = 0

    do i = 1, nfields
      key = compact_key(fields(i))
      select case (trim(key))
      case ('longitude', 'longitudedeg', 'lon', 'londeg')
        idx_lon = i
      case ('latitude', 'latitudedeg', 'lat', 'latdeg')
        idx_lat = i
      case ('altitude', 'altitudem', 'altm', 'heightm', 'elevationm')
        idx_alt_m = i
      case ('altitudekm', 'altkm', 'heightkm', 'elevationkm')
        idx_alt_km = i
      case ('radius', 'radiusm')
        idx_rad_m = i
      case ('radiuskm')
        idx_rad_km = i
      case ('bx', 'bxnt')
        idx_bx_plain_nt = i
      case ('by', 'bynt')
        idx_by_plain_nt = i
      case ('bz', 'bznt')
        idx_bz_plain_nt = i
      case ('bxsel')
        idx_bx_sel_nt = i
      case ('bysel')
        idx_by_sel_nt = i
      case ('bzsel')
        idx_bz_sel_nt = i
      case ('bxsse')
        idx_bx_sse_nt = i
      case ('bysse')
        idx_by_sse_nt = i
      case ('bzsse')
        idx_bz_sse_nt = i
      case ('bxt')
        idx_bx_t = i
      case ('byt')
        idx_by_t = i
      case ('bzt')
        idx_bz_t = i
      end select
    end do

    idx_bx_nt = first_nonzero_int(idx_bx_plain_nt, idx_bx_sel_nt, idx_bx_sse_nt)
    idx_by_nt = first_nonzero_int(idx_by_plain_nt, idx_by_sel_nt, idx_by_sse_nt)
    idx_bz_nt = first_nonzero_int(idx_bz_plain_nt, idx_bz_sel_nt, idx_bz_sse_nt)
  end subroutine locate_csv_columns

  subroutine parse_observation_fields( &
    rsphere_m_in, path, fields, nfields, &
    idx_lon, idx_lat, idx_alt_m, idx_alt_km, idx_rad_m, idx_rad_km, &
    idx_bx_nt, idx_by_nt, idx_bz_nt, idx_bx_t, idx_by_t, idx_bz_t, &
    lonv, latv, radv, bxv_t, byv_t, bzv_t &
  )
    real(wp), intent(in) :: rsphere_m_in
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: fields(:)
    integer(int32), intent(in) :: nfields
    integer(int32), intent(in) :: idx_lon, idx_lat, idx_alt_m, idx_alt_km, idx_rad_m, idx_rad_km
    integer(int32), intent(in) :: idx_bx_nt, idx_by_nt, idx_bz_nt, idx_bx_t, idx_by_t, idx_bz_t
    real(wp), intent(out) :: lonv, latv, radv, bxv_t, byv_t, bzv_t

    lonv = read_csv_real(fields, nfields, idx_lon, path)
    latv = read_csv_real(fields, nfields, idx_lat, path)
    if (abs(latv) > 90.0_wp) stop 'ERROR: latitude outside [-90,90] in observation CSV'

    if (idx_rad_m > 0) then
      radv = read_csv_real(fields, nfields, idx_rad_m, path)
    else if (idx_rad_km > 0) then
      radv = 1000.0_wp * read_csv_real(fields, nfields, idx_rad_km, path)
    else if (idx_alt_m > 0) then
      radv = rsphere_m_in + read_csv_real(fields, nfields, idx_alt_m, path)
    else
      radv = rsphere_m_in + 1000.0_wp * read_csv_real(fields, nfields, idx_alt_km, path)
    end if
    if (radv <= 0.0_wp) stop 'ERROR: non-positive radius in observation CSV'

    if (idx_bx_t > 0) then
      bxv_t = read_csv_real(fields, nfields, idx_bx_t, path)
    else
      bxv_t = read_csv_real(fields, nfields, idx_bx_nt, path) / T_to_nT
    end if
    if (idx_by_t > 0) then
      byv_t = read_csv_real(fields, nfields, idx_by_t, path)
    else
      byv_t = read_csv_real(fields, nfields, idx_by_nt, path) / T_to_nT
    end if
    if (idx_bz_t > 0) then
      bzv_t = read_csv_real(fields, nfields, idx_bz_t, path)
    else
      bzv_t = read_csv_real(fields, nfields, idx_bz_nt, path) / T_to_nT
    end if
  end subroutine parse_observation_fields

  real(wp) function read_csv_real(fields, nfields, idx, path)
    character(len=*), intent(in) :: fields(:)
    integer(int32), intent(in) :: nfields, idx
    character(len=*), intent(in) :: path
    integer(int32) :: ios_local
    if (idx <= 0 .or. idx > nfields) stop 'ERROR: requested CSV column index is invalid'
    read(fields(idx),*,iostat=ios_local) read_csv_real
    if (ios_local /= 0) then
      write(*,'(A,A,A)') 'ERROR: could not parse numeric CSV value in ', trim(path), '.'
      stop 2
    end if
  end function read_csv_real

  subroutine split_path_list(list_in, paths, npaths)
    character(len=*), intent(in) :: list_in
    character(len=1024), allocatable, intent(out) :: paths(:)
    integer(int32), intent(out) :: npaths
    character(len=256), allocatable :: fields(:)
    integer(int32) :: nfields, i, nkeep

    call split_csv_line(trim(list_in), fields, nfields)
    nkeep = 0
    do i = 1, nfields
      if (len_trim(fields(i)) > 0) nkeep = nkeep + 1
    end do
    if (nkeep <= 0) then
      npaths = 0
      allocate(paths(1))
      paths(1) = ''
      return
    end if

    allocate(paths(nkeep))
    nkeep = 0
    do i = 1, nfields
      if (len_trim(fields(i)) <= 0) cycle
      nkeep = nkeep + 1
      paths(nkeep) = adjustl(trim(fields(i)))
    end do
    npaths = nkeep
  end subroutine split_path_list

  subroutine split_csv_line(line, fields, nfields)
    character(len=*), intent(in) :: line
    character(len=256), allocatable, intent(out) :: fields(:)
    integer(int32), intent(out) :: nfields

    integer(int32) :: i, start, nf, lenline

    lenline = len_trim(line)
    if (lenline <= 0) then
      nfields = 0
      allocate(fields(1))
      fields(1) = ''
      return
    end if

    nf = 1
    do i = 1, lenline
      if (line(i:i) == ',') nf = nf + 1
    end do

    allocate(fields(nf))
    fields(:) = ''
    start = 1
    nfields = 0
    do i = 1, lenline + 1
      if (i > lenline .or. line(i:i) == ',') then
        nfields = nfields + 1
        if (i > start) then
          fields(nfields) = strip_csv_quotes(adjustl(trim(line(start:i-1))))
        else
          fields(nfields) = ''
        end if
        start = i + 1
      end if
    end do
  end subroutine split_csv_line

  character(len=256) function strip_csv_quotes(text)
    character(len=*), intent(in) :: text
    integer(int32) :: lt
    strip_csv_quotes = adjustl(trim(text))
    lt = len_trim(strip_csv_quotes)
    if (lt >= 2) then
      if ((strip_csv_quotes(1:1) == '"' .and. strip_csv_quotes(lt:lt) == '"') .or. &
          (strip_csv_quotes(1:1) == '''' .and. strip_csv_quotes(lt:lt) == '''')) then
        strip_csv_quotes = strip_csv_quotes(2:lt-1)
      end if
    end if
    strip_csv_quotes = adjustl(trim(strip_csv_quotes))
  end function strip_csv_quotes

  character(len=256) function compact_key(text)
    character(len=*), intent(in) :: text
    integer(int32) :: i, k, code
    character(len=1) :: ch

    compact_key = ''
    k = 0
    do i = 1, len_trim(text)
      ch = text(i:i)
      code = iachar(ch)
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        k = k + 1
        compact_key(k:k) = achar(code + 32)
      else if ((code >= iachar('a') .and. code <= iachar('z')) .or. (code >= iachar('0') .and. code <= iachar('9'))) then
        k = k + 1
        compact_key(k:k) = ch
      end if
      if (k >= len(compact_key)) exit
    end do
  end function compact_key

  subroutine grow_observation_arrays(cap, lon_deg, lat_deg, radius_m, bx_t, by_t, bz_t)
    integer(int32), intent(inout) :: cap
    real(wp), allocatable, intent(inout) :: lon_deg(:), lat_deg(:), radius_m(:), bx_t(:), by_t(:), bz_t(:)

    integer(int32) :: newcap
    real(wp), allocatable :: a1(:), a2(:), a3(:), a4(:), a5(:), a6(:)

    newcap = max(2_int32 * cap, cap + 1024_int32)
    allocate(a1(newcap), a2(newcap), a3(newcap), a4(newcap), a5(newcap), a6(newcap))
    a1(:) = 0.0_wp; a2(:) = 0.0_wp; a3(:) = 0.0_wp
    a4(:) = 0.0_wp; a5(:) = 0.0_wp; a6(:) = 0.0_wp
    a1(1:cap) = lon_deg(1:cap)
    a2(1:cap) = lat_deg(1:cap)
    a3(1:cap) = radius_m(1:cap)
    a4(1:cap) = bx_t(1:cap)
    a5(1:cap) = by_t(1:cap)
    a6(1:cap) = bz_t(1:cap)
    call move_alloc(a1, lon_deg)
    call move_alloc(a2, lat_deg)
    call move_alloc(a3, radius_m)
    call move_alloc(a4, bx_t)
    call move_alloc(a5, by_t)
    call move_alloc(a6, bz_t)
    cap = newcap
  end subroutine grow_observation_arrays

  subroutine shrink_observation_arrays(nkeep, lon_deg, lat_deg, radius_m, bx_t, by_t, bz_t)
    integer(int32), intent(in) :: nkeep
    real(wp), allocatable, intent(inout) :: lon_deg(:), lat_deg(:), radius_m(:), bx_t(:), by_t(:), bz_t(:)

    real(wp), allocatable :: a1(:), a2(:), a3(:), a4(:), a5(:), a6(:)
    allocate(a1(nkeep), a2(nkeep), a3(nkeep), a4(nkeep), a5(nkeep), a6(nkeep))
    a1(:) = lon_deg(1:nkeep)
    a2(:) = lat_deg(1:nkeep)
    a3(:) = radius_m(1:nkeep)
    a4(:) = bx_t(1:nkeep)
    a5(:) = by_t(1:nkeep)
    a6(:) = bz_t(1:nkeep)
    call move_alloc(a1, lon_deg)
    call move_alloc(a2, lat_deg)
    call move_alloc(a3, radius_m)
    call move_alloc(a4, bx_t)
    call move_alloc(a5, by_t)
    call move_alloc(a6, bz_t)
  end subroutine shrink_observation_arrays

  real(wp) function circular_mean_deg(lon_deg)
    real(wp), intent(in) :: lon_deg(:)
    circular_mean_deg = atan2(sum(sin(lon_deg * deg2rad)), sum(cos(lon_deg * deg2rad))) / deg2rad
  end function circular_mean_deg

  real(wp) function median_positive_spacing(values)
    real(wp), intent(in) :: values(:)
    real(wp), allocatable :: sorted(:), diffs(:)
    integer(int32) :: n, i, nd

    n = size(values)
    if (n <= 1) then
      median_positive_spacing = 0.0_wp
      return
    end if

    allocate(sorted(n))
    sorted(:) = values(:)
    call quicksort_real(sorted, 1_int32, n)

    allocate(diffs(n-1))
    nd = 0
    do i = 1, n - 1
      if (abs(sorted(i+1) - sorted(i)) <= 1.0e-8_wp) cycle
      nd = nd + 1
      diffs(nd) = sorted(i+1) - sorted(i)
    end do
    if (nd <= 0) then
      median_positive_spacing = 0.0_wp
      return
    end if
    call quicksort_real(diffs, 1_int32, nd)
    if (mod(nd, 2) == 1) then
      median_positive_spacing = diffs((nd + 1) / 2)
    else
      median_positive_spacing = 0.5_wp * (diffs(nd / 2) + diffs(nd / 2 + 1))
    end if
  end function median_positive_spacing

  recursive subroutine quicksort_real(arr, left, right)
    real(wp), intent(inout) :: arr(:)
    integer(int32), intent(in) :: left, right
    integer(int32) :: i, j
    real(wp) :: pivot, temp

    if (left >= right) return
    pivot = arr((left + right) / 2)
    i = left
    j = right
    do
      do while (arr(i) < pivot)
        i = i + 1
      end do
      do while (arr(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        temp = arr(i)
        arr(i) = arr(j)
        arr(j) = temp
        i = i + 1
        j = j - 1
      end if
      if (i > j) exit
    end do
    if (left < j) call quicksort_real(arr, left, j)
    if (i < right) call quicksort_real(arr, i, right)
  end subroutine quicksort_real

  subroutine sph_to_cart(radius_m, lat_deg, lon_deg, x, y, z)
    real(wp), intent(in) :: radius_m, lat_deg, lon_deg
    real(wp), intent(out) :: x, y, z
    real(wp) :: lat_rad, lon_rad
    lat_rad = lat_deg * deg2rad
    lon_rad = wrap180(lon_deg) * deg2rad
    x = radius_m * cos(lat_rad) * cos(lon_rad)
    y = radius_m * cos(lat_rad) * sin(lon_rad)
    z = radius_m * sin(lat_rad)
  end subroutine sph_to_cart

  subroutine dipole_kernel_block(xo, yo, zo, xs, ys, zs, K)
    real(wp), intent(in) :: xo, yo, zo, xs, ys, zs
    real(wp), intent(out) :: K(3,3)

    real(wp) :: rx, ry, rz, r2, invr, invr3, invr5

    rx = xo - xs
    ry = yo - ys
    rz = zo - zs
    r2 = rx*rx + ry*ry + rz*rz
    if (r2 <= 0.0_wp) stop 'ERROR: observation coincides with source dipole location'

    invr = 1.0_wp / sqrt(r2)
    invr3 = invr * invr * invr
    invr5 = invr3 * invr * invr

    K(1,1) = mu0_over_4pi * (3.0_wp * rx * rx * invr5 - invr3)
    K(1,2) = mu0_over_4pi * (3.0_wp * rx * ry * invr5)
    K(1,3) = mu0_over_4pi * (3.0_wp * rx * rz * invr5)
    K(2,1) = K(1,2)
    K(2,2) = mu0_over_4pi * (3.0_wp * ry * ry * invr5 - invr3)
    K(2,3) = mu0_over_4pi * (3.0_wp * ry * rz * invr5)
    K(3,1) = K(1,3)
    K(3,2) = K(2,3)
    K(3,3) = mu0_over_4pi * (3.0_wp * rz * rz * invr5 - invr3)
  end subroutine dipole_kernel_block

  subroutine write_prediction_csv( &
    path, rsphere_m_in, nobs_in, lon_deg, lat_deg, radius_m, &
    obs_bx_t, obs_by_t, obs_bz_t, pred_bx_t, pred_by_t, pred_bz_t &
  )
    character(len=*), intent(in) :: path
    real(wp), intent(in) :: rsphere_m_in
    integer(int32), intent(in) :: nobs_in
    real(wp), intent(in) :: lon_deg(:), lat_deg(:), radius_m(:), obs_bx_t(:), obs_by_t(:), obs_bz_t(:)
    real(wp), intent(in) :: pred_bx_t(:), pred_by_t(:), pred_bz_t(:)

    integer(int32) :: unit, i
    real(wp) :: btot_obs_nt, btot_pred_nt

    open(newunit=unit, file=trim(path), status='replace', action='write')
    write(unit,'(A)') &
      'longitude_deg,latitude_deg,radius_m,altitude_m,' // &
      'bx_obs_nt,by_obs_nt,bz_obs_nt,' // &
      'bx_pred_nt,by_pred_nt,bz_pred_nt,' // &
      'bx_resid_nt,by_resid_nt,bz_resid_nt,' // &
      'btot_obs_nt,btot_pred_nt,btot_resid_nt'
    do i = 1, nobs_in
      btot_obs_nt = T_to_nT * sqrt(obs_bx_t(i)**2 + obs_by_t(i)**2 + obs_bz_t(i)**2)
      btot_pred_nt = T_to_nT * sqrt(pred_bx_t(i)**2 + pred_by_t(i)**2 + pred_bz_t(i)**2)
      write(unit,'(*(g0,:,","))') &
        wrap180(lon_deg(i)), lat_deg(i), radius_m(i), radius_m(i) - rsphere_m_in, &
        T_to_nT * obs_bx_t(i), T_to_nT * obs_by_t(i), T_to_nT * obs_bz_t(i), &
        T_to_nT * pred_bx_t(i), T_to_nT * pred_by_t(i), T_to_nT * pred_bz_t(i), &
        T_to_nT * (pred_bx_t(i) - obs_bx_t(i)), &
        T_to_nT * (pred_by_t(i) - obs_by_t(i)), &
        T_to_nT * (pred_bz_t(i) - obs_bz_t(i)), &
        btot_obs_nt, btot_pred_nt, btot_pred_nt - btot_obs_nt
    end do
    close(unit)
  end subroutine write_prediction_csv

  subroutine write_dipole_csv( &
    path, rsphere_m_in, min_obs_radius_m_in, nsrc_in, nlat_in, nlon_in, nr_in, &
    src_lon_deg, src_lat_deg, src_radius_m, model_in &
  )
    character(len=*), intent(in) :: path
    real(wp), intent(in) :: rsphere_m_in, min_obs_radius_m_in
    integer(int32), intent(in) :: nsrc_in, nlat_in, nlon_in, nr_in
    real(wp), intent(in) :: src_lon_deg(:), src_lat_deg(:), src_radius_m(:), model_in(:)

    integer(int32) :: unit, i, layer_id
    real(wp) :: mx, my, mz

    open(newunit=unit, file=trim(path), status='replace', action='write')
    write(unit,'(A)') &
      'source_id,layer_id,lon_deg,lat_deg,radius_m,altitude_m,' // &
      'depth_below_min_obs_km,mx_am2,my_am2,mz_am2,moment_norm_am2'
    do i = 1, nsrc_in
      layer_id = ((i - 1) / (nlat_in * nlon_in)) + 1
      mx = model_in(3*(i-1) + 1)
      my = model_in(3*(i-1) + 2)
      mz = model_in(3*(i-1) + 3)
      write(unit,'(*(g0,:,","))') &
        i, layer_id, wrap180(src_lon_deg(i)), src_lat_deg(i), &
        src_radius_m(i), src_radius_m(i) - rsphere_m_in, &
        (min_obs_radius_m_in - src_radius_m(i)) / 1000.0_wp, &
        mx, my, mz, sqrt(mx*mx + my*my + mz*mz)
    end do
    close(unit)
  end subroutine write_dipole_csv

  pure real(wp) function wrap180(lon_deg)
    real(wp), intent(in) :: lon_deg
    wrap180 = modulo(lon_deg + 180.0_wp, 360.0_wp) - 180.0_wp
  end function wrap180

  pure integer(int32) function first_nonzero_int(a, b, c)
    integer(int32), intent(in) :: a, b, c
    first_nonzero_int = a
    if (first_nonzero_int <= 0) first_nonzero_int = b
    if (first_nonzero_int <= 0) first_nonzero_int = c
  end function first_nonzero_int

  pure real(wp) function unwrap_lon_to_ref(lon_deg, lon_ref_deg)
    real(wp), intent(in) :: lon_deg, lon_ref_deg
    unwrap_lon_to_ref = wrap180(lon_deg)
    do while (unwrap_lon_to_ref - lon_ref_deg >= 180.0_wp)
      unwrap_lon_to_ref = unwrap_lon_to_ref - 360.0_wp
    end do
    do while (unwrap_lon_to_ref - lon_ref_deg < -180.0_wp)
      unwrap_lon_to_ref = unwrap_lon_to_ref + 360.0_wp
    end do
  end function unwrap_lon_to_ref

  pure real(wp) function bytes_to_mib(nbytes)
    integer(int64), intent(in) :: nbytes
    bytes_to_mib = real(nbytes, wp) / (1024.0_wp * 1024.0_wp)
  end function bytes_to_mib

  function btot_vector(bx, by, bz) result(out)
    real(wp), intent(in) :: bx(:), by(:), bz(:)
    real(wp), allocatable :: out(:)
    allocate(out(size(bx)))
    out(:) = sqrt(bx*bx + by*by + bz*bz)
  end function btot_vector

  subroutine solve_linear_system(A, b, x, n, info)
    integer(int32), intent(in) :: n
    real(wp), intent(inout) :: A(n,n)
    real(wp), intent(inout) :: b(n)
    real(wp), intent(out) :: x(n)
    integer(int32), intent(out) :: info

    integer(int32) :: i, j, k, piv
    real(wp) :: maxv, tmp, factor, matrix_scale, pivot_tol
    real(wp), allocatable :: rowtmp(:)

    allocate(rowtmp(n))
    info = 0
    matrix_scale = maxval(abs(A))
    if (matrix_scale <= 0.0_wp) then
      info = 1
      deallocate(rowtmp)
      return
    end if
    pivot_tol = epsilon(1.0_wp) * matrix_scale * real(max(1_int32, n), wp)

    do k = 1, n - 1
      piv = k
      maxv = abs(A(k,k))
      do i = k + 1, n
        if (abs(A(i,k)) > maxv) then
          maxv = abs(A(i,k))
          piv = i
        end if
      end do

      if (maxv < pivot_tol) then
        info = 1
        deallocate(rowtmp)
        return
      end if

      if (piv /= k) then
        rowtmp(:) = A(k,:)
        A(k,:) = A(piv,:)
        A(piv,:) = rowtmp(:)
        tmp = b(k)
        b(k) = b(piv)
        b(piv) = tmp
      end if

      do i = k + 1, n
        factor = A(i,k) / A(k,k)
        A(i,k) = 0.0_wp
        do j = k + 1, n
          A(i,j) = A(i,j) - factor * A(k,j)
        end do
        b(i) = b(i) - factor * b(k)
      end do
    end do

    if (abs(A(n,n)) < pivot_tol) then
      info = 1
      deallocate(rowtmp)
      return
    end if

    x(n) = b(n) / A(n,n)
    do i = n - 1, 1, -1
      tmp = b(i)
      do j = i + 1, n
        tmp = tmp - A(i,j) * x(j)
      end do
      if (abs(A(i,i)) < pivot_tol) then
        info = 1
        deallocate(rowtmp)
        return
      end if
      x(i) = tmp / A(i,i)
    end do

    deallocate(rowtmp)
  end subroutine solve_linear_system

end program gravmag_sphere_dipole_grid_fit
