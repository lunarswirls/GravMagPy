program gravmag_sphere_gauss
  use, intrinsic :: iso_fortran_env, only: int32, real64
  implicit none

  !***********************************************************************
  ! gravmag_sphere_gauss
  !
  ! spectral alternative to the direct volume solver:
  ! - build elemental dipoles (magnetic) or point masses (gravity)
  ! - fit spherical-harmonic coefficients from radial component on a fitting sphere
  ! - reconstruct field components on the observation grid
  !
  ! old-vs-new note:
  ! - this is not the original legacy gauss-legendre integration path!
  ! - this file implements a modern least-squares spherical-harmonic fit
  ! - for compact/sharp sources, low lmax can smooth and spread anomalies
  !     in those cases prefer gravmag_sphere_bxyz (direct volume sum)
  !***********************************************************************

  integer, parameter :: wp = real64
  real(wp), parameter :: pi = 3.1415926535897932384626433832795_wp
  real(wp), parameter :: deg2rad = pi/180.0_wp
  real(wp), parameter :: mu0_over_4pi = 1.0e-7_wp
  real(wp), parameter :: grav_G = 6.67430e-11_wp
  real(wp), parameter :: T_to_nT = 1.0e9_wp
  real(wp), parameter :: mps2_to_mgal = 1.0e5_wp

  character(len=256) :: infile, outfile, tmp
  integer(int32) :: ios, argc

  real(wp) :: rsphere_km, rsphere_m

  ! solver controls (optional CLI overrides)
  integer(int32) :: lmax, refine_factor, ntheta_fit, nphi_fit
  integer(int32) :: src_nlat_cli, src_nlon_cli, src_nr_cli
  real(wp) :: reg_lambda, reg_power

  ! per-body input cards
  character(len=256) :: title
  integer(int32) :: knk

  real(wp) :: lat0_deg, lon0_deg, dlat_deg, dlon_deg, elvo_km
  integer(int32) :: nlat, nlon

  integer(int32) :: nr, ntheta, nphi, nblim

  real(wp) :: htheta, hphi, phi1, phi2, rho_dummy

  integer(int32) :: ifield, nrem
  real(wp) :: M_amp_Apm, inc_deg, dec_deg
  real(wp) :: rho_kgm3

  integer(int32) :: iprint, nfile, nopt
  real(wp) :: first, conint

  integer(int32) :: npts
  real(wp) :: depth_top_km, depth_bot_km
  real(wp) :: lat_max, lat_min, lon_max, lon_min
  real(wp), allocatable :: poly_lat(:), poly_lon(:)

  ! geometry & arrays
  real(wp) :: ro_m, r_top_m, r_bot_m, dr_m
  real(wp) :: cosI, sinI, cosD, sinD
  real(wp) :: Mx, My, Mz

  real(wp), allocatable :: xo(:,:), yo(:,:), zo(:,:)
  real(wp), allocatable :: Bx(:,:), By(:,:), Bz(:,:)

  real(wp) :: lat_lo, lat_hi, lon_lo, lon_hi, lon_ref
  integer(int32) :: nlat_s, nlon_s, nlat_seed, nlon_seed, nr_s
  real(wp) :: dlat_s_deg, dlon_s_deg

  integer(int32) :: ns, cap
  real(wp), allocatable :: xs_list(:), ys_list(:), zs_list(:)
  real(wp), allocatable :: mx_list(:), my_list(:), mz_list(:)
  real(wp), allocatable :: mass_list(:)

  integer(int32) :: i, j, is, it, ir, l, m
  real(wp) :: lat_deg, lon_deg, lat_rad, lon_rad
  real(wp) :: xs, ys, zs, dV, r_mid

  ! Gauss coefficient fit
  integer(int32), allocatable :: idx_g(:,:), idx_h(:,:)
  integer(int32) :: ncoef
  real(wp), allocatable :: ATA(:,:), ATb(:), coef(:), basis(:)

  real(wp) :: theta, phi, st, ct, cp, sp
  real(wp) :: xb, yb, zb
  real(wp) :: dBx, dBy, dBz
  real(wp) :: br_fit, w
  integer(int32) :: itheta, iphi

  real(wp), allocatable :: P(:,:), gh_g(:,:), gh_h(:,:)
  real(wp) :: fac, common, dP, plm1, s_safe
  real(wp) :: br, bt, bp, bxv, byv, bzv, btot_nt, gtot_mgal
  real(wp) :: reg_w, lnorm

  ! ----------------------------
  ! CLI
  ! ----------------------------
  argc = command_argument_count()
  if (argc < 3) then
    write(*,*) 'Usage: gravmag_sphere_gauss <R_sphere_km> <input.in> <output.txt> '// &
               '[lmax] [refine_factor] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] '// &
               '[source_nlat] [source_nlon] [source_nr]'
    stop 2
  end if

  call get_command_argument(1, tmp)
  read(tmp,*,iostat=ios) rsphere_km
  if (ios /= 0) stop 'ERROR: could not parse R_sphere_km'

  call get_command_argument(2, infile)
  call get_command_argument(3, outfile)

  lmax = 24
  refine_factor = 2
  ntheta_fit = 72
  nphi_fit = 144
  reg_lambda = 0.2_wp
  reg_power = 4.0_wp
  src_nlat_cli = 0_int32
  src_nlon_cli = 0_int32
  src_nr_cli = 0_int32

  if (argc >= 4) then
    call get_command_argument(4, tmp)
    read(tmp,*,iostat=ios) lmax
    if (ios /= 0 .or. lmax < 1) stop 'ERROR: invalid lmax'
  end if
  if (argc >= 5) then
    call get_command_argument(5, tmp)
    read(tmp,*,iostat=ios) refine_factor
    if (ios /= 0 .or. refine_factor < 1) stop 'ERROR: invalid refine_factor'
  end if
  if (argc >= 6) then
    call get_command_argument(6, tmp)
    read(tmp,*,iostat=ios) ntheta_fit
    if (ios /= 0 .or. ntheta_fit < 4) stop 'ERROR: invalid ntheta_fit'
  end if
  if (argc >= 7) then
    call get_command_argument(7, tmp)
    read(tmp,*,iostat=ios) nphi_fit
    if (ios /= 0 .or. nphi_fit < 8) stop 'ERROR: invalid nphi_fit'
  end if
  if (argc >= 8) then
    call get_command_argument(8, tmp)
    read(tmp,*,iostat=ios) reg_lambda
    if (ios /= 0 .or. reg_lambda < 0.0_wp) stop 'ERROR: invalid reg_lambda'
  end if
  if (argc >= 9) then
    call get_command_argument(9, tmp)
    read(tmp,*,iostat=ios) reg_power
    if (ios /= 0 .or. reg_power < 0.0_wp) stop 'ERROR: invalid reg_power'
  end if
  if (argc >= 10) then
    call get_command_argument(10, tmp)
    read(tmp,*,iostat=ios) src_nlat_cli
    if (ios /= 0 .or. src_nlat_cli < 0) stop 'ERROR: invalid source_nlat'
  end if
  if (argc >= 11) then
    call get_command_argument(11, tmp)
    read(tmp,*,iostat=ios) src_nlon_cli
    if (ios /= 0 .or. src_nlon_cli < 0) stop 'ERROR: invalid source_nlon'
  end if
  if (argc >= 12) then
    call get_command_argument(12, tmp)
    read(tmp,*,iostat=ios) src_nr_cli
    if (ios /= 0 .or. src_nr_cli < 0) stop 'ERROR: invalid source_nr'
  end if

  open(4, file=trim(infile), status='old', iostat=ios)
  if (ios /= 0) stop 'ERROR: cannot open input file'

  open(2, file=trim(outfile), status='replace', iostat=ios)
  if (ios /= 0) stop 'ERROR: cannot open output file'

  rsphere_m = rsphere_km * 1000.0_wp
  knk = 0

  allocate(idx_g(0:lmax,0:lmax), idx_h(0:lmax,0:lmax))
  call build_index_maps(lmax, idx_g, idx_h, ncoef)
  allocate(ATA(ncoef,ncoef), ATb(ncoef), coef(ncoef), basis(ncoef))
  allocate(P(0:lmax,0:lmax), gh_g(0:lmax,0:lmax), gh_h(0:lmax,0:lmax))

  do
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
    write(*,'(A,I0,2A)') 'Body ', knk, ': ', trim(title)

    read(4,*,iostat=ios) lat0_deg, lon0_deg, dlat_deg, dlon_deg, elvo_km, nlat, nlon
    if (ios /= 0) stop 'ERROR reading Card 2'

    read(4,*,iostat=ios) nr, ntheta, nphi, nblim
    if (ios /= 0) stop 'ERROR reading Card 3'

    read(4,*,iostat=ios) htheta, hphi, phi1, phi2, rho_dummy
    if (ios /= 0) stop 'ERROR reading Card 4'

    read(4,*,iostat=ios) ifield, nrem, M_amp_Apm, inc_deg, dec_deg
    if (ios /= 0) stop 'ERROR reading Card 5'
    if (ifield /= 1 .and. ifield /= 2) stop 'ERROR: Card 5 ifield must be 1(gravity) or 2(magnetic)'
    if (ifield == 2 .and. nrem /= 1) stop 'ERROR: magnetics require nrem=1'

    read(4,*,iostat=ios) iprint, nfile, nopt, first, conint
    if (ios /= 0) stop 'ERROR reading Card 6'

    if (ifield == 2) then
      cosI = cos(inc_deg*deg2rad)
      sinI = sin(inc_deg*deg2rad)
      cosD = cos(dec_deg*deg2rad)
      sinD = sin(dec_deg*deg2rad)
      Mx = M_amp_Apm * (cosI*cosD)
      My = M_amp_Apm * (cosI*sinD)
      Mz = M_amp_Apm * (sinI)
      rho_kgm3 = 0.0_wp
    else
      Mx = 0.0_wp
      My = 0.0_wp
      Mz = 0.0_wp
      rho_kgm3 = M_amp_Apm
    end if

    if (allocated(poly_lat)) deallocate(poly_lat, poly_lon)

    if (nblim == 1) then
      read(4,*,iostat=ios) lat_max, lat_min, lon_max, lon_min, depth_top_km, depth_bot_km
      if (ios /= 0) stop 'ERROR reading fixed limits (Card 7)'
      npts = 5
      allocate(poly_lat(npts), poly_lon(npts))
      poly_lat = [lat_min, lat_max, lat_max, lat_min, lat_min]
      poly_lon = [lon_min, lon_min, lon_max, lon_max, lon_min]
    else
      read(4,*,iostat=ios) npts, depth_top_km, depth_bot_km
      if (ios /= 0) stop 'ERROR reading polygon header (Card 7)'
      if (npts < 3) stop 'ERROR: polygon needs >=3 vertices'
      allocate(poly_lat(npts), poly_lon(npts))
      do i=1,npts
        read(4,*,iostat=ios) poly_lat(i), poly_lon(i)
        if (ios /= 0) stop 'ERROR reading polygon vertex'
      end do
    end if

    lon_ref = sum(poly_lon) / real(size(poly_lon), wp)
    call wrap180_scalar(lon_ref)
    call unwrap_lons_inplace(poly_lon, lon_ref)

    lat_lo = minval(poly_lat); lat_hi = maxval(poly_lat)
    lon_lo = minval(poly_lon); lon_hi = maxval(poly_lon)

    if (depth_top_km < 0.0_wp .or. depth_bot_km < 0.0_wp) stop 'ERROR: depth must be >= 0'
    if (depth_bot_km <= depth_top_km) stop 'ERROR: need depth_bot > depth_top'

    r_top_m = rsphere_m - depth_top_km*1000.0_wp
    r_bot_m = rsphere_m - depth_bot_km*1000.0_wp
    if (r_top_m <= 0.0_wp .or. r_bot_m <= 0.0_wp) stop 'ERROR: depth too large'

    if (r_bot_m > r_top_m) then
      dr_m = r_bot_m
      r_bot_m = r_top_m
      r_top_m = dr_m
    end if
    dr_m = r_top_m - r_bot_m

    ro_m = (rsphere_km + elvo_km) * 1000.0_wp

    if (.not. allocated(xo)) then
      allocate(xo(nlat,nlon), yo(nlat,nlon), zo(nlat,nlon))
      allocate(Bx(nlat,nlon), By(nlat,nlon), Bz(nlat,nlon))
    else
      if (size(xo,1) /= nlat .or. size(xo,2) /= nlon) then
        deallocate(xo,yo,zo,Bx,By,Bz)
        allocate(xo(nlat,nlon), yo(nlat,nlon), zo(nlat,nlon))
        allocate(Bx(nlat,nlon), By(nlat,nlon), Bz(nlat,nlon))
      end if
    end if

    do i=1,nlat
      do j=1,nlon
        lat_deg = lat0_deg + real(i-1,wp)*dlat_deg
        lon_deg = lon0_deg + real(j-1,wp)*dlon_deg
        call wrap180_scalar(lon_deg)

        lat_rad = lat_deg*deg2rad
        lon_rad = lon_deg*deg2rad

        xo(i,j) = ro_m * cos(lat_rad) * cos(lon_rad)
        yo(i,j) = ro_m * cos(lat_rad) * sin(lon_rad)
        zo(i,j) = ro_m * sin(lat_rad)
      end do
    end do

    ! Source mesh is decoupled from Card 2 output grid:
    ! - default source controls come from Card 3 (ntheta, nphi, nr)
    ! - optional CLI overrides can replace them globally
    nlat_seed = ntheta
    nlon_seed = nphi
    nr_s = nr
    if (src_nlat_cli > 0) nlat_seed = src_nlat_cli
    if (src_nlon_cli > 0) nlon_seed = src_nlon_cli
    if (src_nr_cli > 0) nr_s = src_nr_cli
    if (nlat_seed <= 0) nlat_seed = max(1_int32, nlat)
    if (nlon_seed <= 0) nlon_seed = max(1_int32, nlon)
    if (nr_s <= 0) nr_s = 1

    nlat_s = max(1_int32, refine_factor * nlat_seed)
    nlon_s = max(1_int32, refine_factor * nlon_seed)

    dlat_s_deg = (lat_hi - lat_lo) / real(nlat_s, wp)
    dlon_s_deg = (lon_hi - lon_lo) / real(nlon_s, wp)

    cap = max(1024_int32, nlat_s*nlon_s*nr_s)
    ns  = 0
    if (allocated(xs_list)) deallocate(xs_list, ys_list, zs_list, mx_list, my_list, mz_list, mass_list)
    allocate(xs_list(cap), ys_list(cap), zs_list(cap))
    allocate(mx_list(cap), my_list(cap), mz_list(cap))
    allocate(mass_list(cap))

    do it=1, nlat_s
      lat_deg = lat_lo + (real(it,wp)-0.5_wp)*dlat_s_deg
      if (lat_deg <= -90.0_wp .or. lat_deg >= 90.0_wp) cycle

      do is=1, nlon_s
        lon_deg = lon_lo + (real(is,wp)-0.5_wp)*dlon_s_deg

        if (.not. point_in_poly_lonlat(lon_deg, lat_deg, poly_lon, poly_lat)) cycle

        lat_rad = lat_deg*deg2rad
        lon_rad = lon_deg*deg2rad

        do ir=1, nr_s
          r_mid = r_bot_m + ((real(ir,wp)-0.5_wp)/real(nr_s,wp))*dr_m

          dV = (r_mid*r_mid) * cos(lat_rad) * (dlat_s_deg*deg2rad) * (dlon_s_deg*deg2rad) * (dr_m/real(nr_s,wp))
          if (dV <= 0.0_wp) cycle

          xs = r_mid * cos(lat_rad) * cos(lon_rad)
          ys = r_mid * cos(lat_rad) * sin(lon_rad)
          zs = r_mid * sin(lat_rad)

          ns = ns + 1
          if (ns > cap) then
            call grow_sources(cap, xs_list, ys_list, zs_list, mx_list, my_list, mz_list, mass_list)
          end if

          xs_list(ns) = xs
          ys_list(ns) = ys
          zs_list(ns) = zs
          if (ifield == 2) then
            mx_list(ns) = Mx * dV
            my_list(ns) = My * dV
            mz_list(ns) = Mz * dV
            mass_list(ns) = 0.0_wp
          else
            mx_list(ns) = 0.0_wp
            my_list(ns) = 0.0_wp
            mz_list(ns) = 0.0_wp
            mass_list(ns) = rho_kgm3 * dV
          end if
        end do
      end do
    end do

    write(*,'(A,I0,A,I0,A,I0)') '  Source mesh cells (nlat,nlon,nr): ', nlat_s, ', ', nlon_s, ', ', nr_s
    write(*,'(A,I0)') '  Source elements: ', ns

    ! ------------------------------------------
    ! Fit Gauss coefficients from Br on fitting sphere
    ! ------------------------------------------
    ATA(:,:) = 0.0_wp
    ATb(:) = 0.0_wp

    do itheta=1, ntheta_fit
      theta = (real(itheta,wp)-0.5_wp) * pi / real(ntheta_fit,wp)
      st = sin(theta)
      ct = cos(theta)

      do iphi=1, nphi_fit
        phi = -pi + (real(iphi-1,wp) * 2.0_wp*pi/real(nphi_fit,wp))
        cp = cos(phi)
        sp = sin(phi)

        xb = ro_m * st * cp
        yb = ro_m * st * sp
        zb = ro_m * ct

        dBx = 0.0_wp; dBy = 0.0_wp; dBz = 0.0_wp
        do is=1, ns
          if (ifield == 2) then
            call add_dipole_field_cart64(xb, yb, zb, xs_list(is), ys_list(is), zs_list(is), &
                                         mx_list(is), my_list(is), mz_list(is), br, bt, bp)
          else
            call add_pointmass_field_cart64(xb, yb, zb, xs_list(is), ys_list(is), zs_list(is), &
                                            mass_list(is), br, bt, bp)
          end if
          dBx = dBx + br
          dBy = dBy + bt
          dBz = dBz + bp
        end do

        br_fit = dBx*(st*cp) + dBy*(st*sp) + dBz*ct

        call fill_basis_br(theta, phi, lmax, idx_g, idx_h, basis)

        w = max(st, 1.0e-10_wp)
        do i=1, ncoef
          ATb(i) = ATb(i) + w * basis(i) * br_fit
          do j=1, ncoef
            ATA(i,j) = ATA(i,j) + w * basis(i) * basis(j)
          end do
        end do
      end do
    end do

    ! Degree-dependent Tikhonov regularization:
    !   ATA <- ATA + lambda * W_l
    ! with stronger damping at higher degree to suppress ringing
    if (reg_lambda > 0.0_wp) then
      do l=1, lmax
        lnorm = real(l,wp) / real(max(1_int32,lmax),wp)
        reg_w = max(lnorm, 1.0e-6_wp)**reg_power
        do m=0, l
          if (idx_g(l,m) > 0) ATA(idx_g(l,m),idx_g(l,m)) = ATA(idx_g(l,m),idx_g(l,m)) * (1.0_wp + reg_lambda*reg_w)
          if (idx_h(l,m) > 0) ATA(idx_h(l,m),idx_h(l,m)) = ATA(idx_h(l,m),idx_h(l,m)) * (1.0_wp + reg_lambda*reg_w)
        end do
      end do
    end if

    do i=1, ncoef
      ATA(i,i) = ATA(i,i) + 1.0e-22_wp
    end do

    call solve_linear_system(ATA, ATb, coef, ncoef, ios)
    if (ios /= 0) stop 'ERROR: failed to solve normal equations for Gauss coefficients'

    gh_g(:,:) = 0.0_wp
    gh_h(:,:) = 0.0_wp
    do i=1, lmax
      do j=0, i
        if (idx_g(i,j) > 0) gh_g(i,j) = coef(idx_g(i,j))
        if (idx_h(i,j) > 0) gh_h(i,j) = coef(idx_h(i,j))
      end do
    end do

    ! ------------------------------------------
    ! Evaluate field from coefficients on obs grid
    ! ------------------------------------------
    do i=1, nlat
      do j=1, nlon
        lat_deg = lat0_deg + real(i-1,wp)*dlat_deg
        lon_deg = lon0_deg + real(j-1,wp)*dlon_deg
        call wrap180_scalar(lon_deg)

        theta = (90.0_wp - lat_deg) * deg2rad
        phi = lon_deg * deg2rad
        st = sin(theta)
        ct = cos(theta)
        cp = cos(phi)
        sp = sin(phi)

        call legendre_point(lmax, ct, P)

        br = 0.0_wp
        bt = 0.0_wp
        bp = 0.0_wp

        s_safe = max(abs(st), 1.0e-10_wp)

        do it=1, lmax
          fac = (ro_m / ro_m)**(it+2)  ! =1 when evaluated at fitting radius
          do ir=0, it
            common = gh_g(it,ir)
            if (ir == 0) then
              br = br + real(it+1,wp) * fac * common * P(it,ir)
              if (it-1 >= ir) then
                plm1 = P(it-1,ir)
              else
                plm1 = 0.0_wp
              end if
              dP = (real(it,wp)*ct*P(it,ir) - real(it+ir,wp)*plm1) / s_safe
              bt = bt - fac * common * dP
            else
              common = gh_g(it,ir)*cos(real(ir,wp)*phi) + gh_h(it,ir)*sin(real(ir,wp)*phi)

              br = br + real(it+1,wp) * fac * common * P(it,ir)

              if (it-1 >= ir) then
                plm1 = P(it-1,ir)
              else
                plm1 = 0.0_wp
              end if
              dP = (real(it,wp)*ct*P(it,ir) - real(it+ir,wp)*plm1) / s_safe
              bt = bt - fac * common * dP

              bp = bp + fac * (real(ir,wp)/s_safe) * &
                   (gh_g(it,ir)*sin(real(ir,wp)*phi) - gh_h(it,ir)*cos(real(ir,wp)*phi)) * P(it,ir)
            end if
          end do
        end do

        bxv = br*st*cp + bt*ct*cp - bp*sp
        byv = br*st*sp + bt*ct*sp + bp*cp
        bzv = br*ct - bt*st

        Bx(i,j) = bxv
        By(i,j) = byv
        Bz(i,j) = bzv
      end do
    end do

    if (nfile /= 0) then
      write(2,'(A)') '# solver=spherical_harmonic_spectral'
      write(2,'(A,I0,A,ES13.5,A,F8.3)') '# lmax=', lmax, ' reg_lambda=', reg_lambda, ' reg_power=', reg_power
      if (ifield == 2) then
        write(2,'(A)') '# body_id lon_deg[-180,180) lat_deg Bx_nT By_nT Bz_nT Btot_nT'
      else
        write(2,'(A)') '# body_id lon_deg[-180,180) lat_deg gx_mGal gy_mGal gz_mGal gtot_mGal'
      end if
      do i=1,nlat
        do j=1,nlon
          lat_deg = lat0_deg + real(i-1,wp)*dlat_deg
          lon_deg = lon0_deg + real(j-1,wp)*dlon_deg
          call wrap180_scalar(lon_deg)

          if (ifield == 2) then
            btot_nt = sqrt((Bx(i,j)*T_to_nT)**2 + (By(i,j)*T_to_nT)**2 + (Bz(i,j)*T_to_nT)**2)
            write(2,'(I6,1X,F12.6,1X,F12.6,1X,ES15.7,1X,ES15.7,1X,ES15.7,1X,ES15.7)') &
              knk, lon_deg, lat_deg, Bx(i,j)*T_to_nT, By(i,j)*T_to_nT, Bz(i,j)*T_to_nT, btot_nt
          else
            gtot_mgal = sqrt((Bx(i,j)*mps2_to_mgal)**2 + (By(i,j)*mps2_to_mgal)**2 + (Bz(i,j)*mps2_to_mgal)**2)
            write(2,'(I6,1X,F12.6,1X,F12.6,1X,ES15.7,1X,ES15.7,1X,ES15.7,1X,ES15.7)') &
              knk, lon_deg, lat_deg, Bx(i,j)*mps2_to_mgal, By(i,j)*mps2_to_mgal, Bz(i,j)*mps2_to_mgal, gtot_mgal
          end if
        end do
      end do
    end if

  end do

  close(4)
  close(2)

contains

  ! keep longitudes on the canonical branch [-180,180)
  subroutine wrap180_scalar(lon_deg)
    real(wp), intent(inout) :: lon_deg
    lon_deg = modulo(lon_deg + 180.0_wp, 360.0_wp) - 180.0_wp
  end subroutine wrap180_scalar

  ! shift longitudes to a continuous branch near lon_ref_deg so polygon
  ! edges do not jump across the dateline during point-in-polygon tests
  subroutine unwrap_lons_inplace(lon_deg, lon_ref_deg)
    real(wp), intent(inout) :: lon_deg(:)
    real(wp), intent(in) :: lon_ref_deg
    integer(int32) :: i
    real(wp) :: x
    do i=1, size(lon_deg)
      x = lon_deg(i)
      call wrap180_scalar(x)
      do while (x - lon_ref_deg >= 180.0_wp)
        x = x - 360.0_wp
      end do
      do while (x - lon_ref_deg < -180.0_wp)
        x = x + 360.0_wp
      end do
      lon_deg(i) = x
    end do
  end subroutine unwrap_lons_inplace

  ! standard ray-casting polygon inclusion in lon-lat space
  ! polygon longitudes are assumed already unwrapped to a continuous branch
  logical function point_in_poly_lonlat(lon, lat, poly_lon, poly_lat)
    real(wp), intent(in) :: lon, lat
    real(wp), intent(in) :: poly_lon(:), poly_lat(:)

    integer(int32) :: i, j, n
    real(wp) :: xi, yi, xj, yj, den, xint
    logical :: inside, intersect

    n = size(poly_lon)
    if (n < 3) then
      point_in_poly_lonlat = .false.
      return
    end if

    inside = .false.
    j = n
    do i=1, n
      xi = poly_lon(i); yi = poly_lat(i)
      xj = poly_lon(j); yj = poly_lat(j)

      intersect = ((yi > lat) .neqv. (yj > lat))
      if (intersect) then
        den = yj - yi
        if (abs(den) < 1.0e-30_wp) den = sign(1.0e-30_wp, den + 1.0e-30_wp)
        xint = (xj - xi) * (lat - yi) / den + xi
        if (lon < xint) inside = .not. inside
      end if

      j = i
    end do

    point_in_poly_lonlat = inside
  end function point_in_poly_lonlat

  ! dynamic growth helper for source-element arrays
  subroutine grow_sources(cap, xs, ys, zs, mx, my, mz, mass)
    integer(int32), intent(inout) :: cap
    real(wp), allocatable, intent(inout) :: xs(:), ys(:), zs(:)
    real(wp), allocatable, intent(inout) :: mx(:), my(:), mz(:)
    real(wp), allocatable, intent(inout) :: mass(:)

    integer(int32) :: newcap
    real(wp), allocatable :: xs2(:), ys2(:), zs2(:), mx2(:), my2(:), mz2(:), mass2(:)

    ! grow by 2x or by a fixed block to reduce reallocation churn
    newcap = max(2*cap, cap + 100000_int32)
    allocate(xs2(newcap), ys2(newcap), zs2(newcap), mx2(newcap), my2(newcap), mz2(newcap), mass2(newcap))

    xs2(1:cap) = xs; ys2(1:cap) = ys; zs2(1:cap) = zs
    mx2(1:cap) = mx; my2(1:cap) = my; mz2(1:cap) = mz
    mass2(1:cap) = mass

    call move_alloc(xs2, xs)
    call move_alloc(ys2, ys)
    call move_alloc(zs2, zs)
    call move_alloc(mx2, mx)
    call move_alloc(my2, my)
    call move_alloc(mz2, mz)
    call move_alloc(mass2, mass)

    cap = newcap
  end subroutine grow_sources

  ! dipole field kernel in double precision
  subroutine add_dipole_field_cart64(xo, yo, zo, xs, ys, zs, mx, my, mz, dBx, dBy, dBz)
    real(wp), intent(in) :: xo, yo, zo
    real(wp), intent(in) :: xs, ys, zs
    real(wp), intent(in) :: mx, my, mz
    real(wp), intent(out) :: dBx, dBy, dBz

    real(wp) :: rx, ry, rz, r2, invr, invr3, invr5, mdotr, c

    rx = xo - xs
    ry = yo - ys
    rz = zo - zs

    r2 = rx*rx + ry*ry + rz*rz
    if (r2 <= 0.0_wp) then
      dBx = 0.0_wp; dBy = 0.0_wp; dBz = 0.0_wp
      return
    end if

    invr = 1.0_wp / sqrt(r2)
    invr3 = invr*invr*invr
    invr5 = invr3*invr*invr

    mdotr = mx*rx + my*ry + mz*rz
    c = 3.0_wp * mdotr * invr5

    dBx = mu0_over_4pi * (c*rx - mx*invr3)
    dBy = mu0_over_4pi * (c*ry - my*invr3)
    dBz = mu0_over_4pi * (c*rz - mz*invr3)
  end subroutine add_dipole_field_cart64

  ! point-mass gravity kernel in double precision
  subroutine add_pointmass_field_cart64(xo, yo, zo, xs, ys, zs, mass_kg, dgx, dgy, dgz)
    real(wp), intent(in) :: xo, yo, zo
    real(wp), intent(in) :: xs, ys, zs
    real(wp), intent(in) :: mass_kg
    real(wp), intent(out) :: dgx, dgy, dgz

    real(wp) :: rx, ry, rz, r2, invr, invr3, c

    rx = xs - xo
    ry = ys - yo
    rz = zs - zo

    r2 = rx*rx + ry*ry + rz*rz
    if (r2 <= 0.0_wp) then
      dgx = 0.0_wp; dgy = 0.0_wp; dgz = 0.0_wp
      return
    end if

    invr = 1.0_wp / sqrt(r2)
    invr3 = invr*invr*invr
    c = grav_G * mass_kg * invr3

    dgx = c * rx
    dgy = c * ry
    dgz = c * rz
  end subroutine add_pointmass_field_cart64

  ! map (l,m) coefficient pairs into packed 1-D vector indices
  ! g_lm and h_lm each get their own index stream except m=0
  ! where h_l0 is not used
  subroutine build_index_maps(lm, ig, ih, ncoef_out)
    integer(int32), intent(in) :: lm
    integer(int32), intent(out) :: ig(0:lm,0:lm), ih(0:lm,0:lm)
    integer(int32), intent(out) :: ncoef_out

    integer(int32) :: l, m, idx

    ig = 0
    ih = 0
    idx = 0
    do l=1, lm
      do m=0, l
        idx = idx + 1
        ig(l,m) = idx
        if (m > 0) then
          idx = idx + 1
          ih(l,m) = idx
        end if
      end do
    end do
    ncoef_out = idx
  end subroutine build_index_maps

  ! compute associated Legendre functions P_lm(x) for one x=cos(theta)
  subroutine legendre_point(lm, x, Pout)
    integer(int32), intent(in) :: lm
    real(wp), intent(in) :: x
    real(wp), intent(out) :: Pout(0:lm,0:lm)

    integer(int32) :: l, m
    real(wp) :: sq

    Pout(:,:) = 0.0_wp
    Pout(0,0) = 1.0_wp

    sq = sqrt(max(0.0_wp, 1.0_wp - x*x))

    do m=1, lm
      Pout(m,m) = -(2.0_wp*real(m,wp)-1.0_wp) * sq * Pout(m-1,m-1)
    end do

    do m=0, lm-1
      Pout(m+1,m) = (2.0_wp*real(m,wp)+1.0_wp) * x * Pout(m,m)
    end do

    do m=0, lm
      do l=m+2, lm
        Pout(l,m) = ((2.0_wp*real(l,wp)-1.0_wp)*x*Pout(l-1,m) - (real(l+m-1,wp))*Pout(l-2,m)) / real(l-m,wp)
      end do
    end do
  end subroutine legendre_point

  ! fill design-vector basis for Br at one (theta,phi) sample
  ! basis matches the packed coefficient ordering from build_index_maps
  subroutine fill_basis_br(theta, phi, lm, ig, ih, basis)
    real(wp), intent(in) :: theta, phi
    integer(int32), intent(in) :: lm
    integer(int32), intent(in) :: ig(0:lm,0:lm), ih(0:lm,0:lm)
    real(wp), intent(out) :: basis(:)

    real(wp) :: x
    real(wp), allocatable :: Pl(:,:)
    integer(int32) :: l, m, idx

    allocate(Pl(0:lm,0:lm))
    x = cos(theta)
    call legendre_point(lm, x, Pl)

    basis(:) = 0.0_wp
    do l=1, lm
      do m=0, l
        idx = ig(l,m)
        basis(idx) = real(l+1,wp) * Pl(l,m) * cos(real(m,wp)*phi)
        if (m > 0) then
          idx = ih(l,m)
          basis(idx) = real(l+1,wp) * Pl(l,m) * sin(real(m,wp)*phi)
        end if
      end do
    end do

    deallocate(Pl)
  end subroutine fill_basis_br

  ! dense Gaussian elimination with partial pivoting
  ! this keeps external dependencies minimal for a self-contained build
  subroutine solve_linear_system(A, b, x, n, info)
    integer(int32), intent(in) :: n
    real(wp), intent(inout) :: A(n,n)
    real(wp), intent(inout) :: b(n)
    real(wp), intent(out) :: x(n)
    integer(int32), intent(out) :: info

    integer(int32) :: i, j, k, piv
    real(wp) :: maxv, tmp, factor
    real(wp), allocatable :: rowtmp(:)

    allocate(rowtmp(n))
    info = 0

    do k=1, n-1
      piv = k
      maxv = abs(A(k,k))
      do i=k+1, n
        if (abs(A(i,k)) > maxv) then
          maxv = abs(A(i,k))
          piv = i
        end if
      end do

      ! near-singular pivot indicates ill-conditioned normal equations
      if (maxv < 1.0e-30_wp) then
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

      do i=k+1, n
        factor = A(i,k) / A(k,k)
        A(i,k) = 0.0_wp
        do j=k+1, n
          A(i,j) = A(i,j) - factor*A(k,j)
        end do
        b(i) = b(i) - factor*b(k)
      end do
    end do

    if (abs(A(n,n)) < 1.0e-30_wp) then
      info = 1
      deallocate(rowtmp)
      return
    end if

    x(n) = b(n) / A(n,n)
    do i=n-1, 1, -1
      tmp = b(i)
      do j=i+1, n
        tmp = tmp - A(i,j)*x(j)
      end do
      if (abs(A(i,i)) < 1.0e-30_wp) then
        info = 1
        deallocate(rowtmp)
        return
      end if
      x(i) = tmp / A(i,i)
    end do

    deallocate(rowtmp)
  end subroutine solve_linear_system

end program gravmag_sphere_gauss
