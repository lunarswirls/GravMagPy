module gravmag_sphere_subs
  use, intrinsic :: iso_fortran_env, only: real32, int32, real64
  implicit none
  private
  public :: wrap180_scalar, unwrap_lons_inplace
  public :: point_in_poly_lonlat
  public :: gauss_legendre

contains

  subroutine wrap180_scalar(x)
    ! Wrap to [-180, 180)
    real(real32), intent(inout) :: x
    x = modulo(x + 180.0_real32, 360.0_real32)
    if (x < 0.0_real32) x = x + 360.0_real32
    x = x - 180.0_real32
    ! keep -180 inclusive, +180 exclusive
    if (x >= 180.0_real32) x = x - 360.0_real32
  end subroutine wrap180_scalar


  subroutine unwrap_lons_inplace(lon_deg, lon_ref_deg)
    ! Shift longitudes by +/-360 so they lie on a continuous branch near lon_ref_deg.
    real(real32), intent(inout) :: lon_deg(:)
    real(real32), intent(in)    :: lon_ref_deg
    integer(int32) :: i
    real(real32) :: x

    do i=1, size(lon_deg)
      x = lon_deg(i)
      call wrap180_scalar(x)

      do while (x - lon_ref_deg >= 180.0_real32)
        x = x - 360.0_real32
      end do
      do while (x - lon_ref_deg < -180.0_real32)
        x = x + 360.0_real32
      end do

      lon_deg(i) = x
    end do
  end subroutine unwrap_lons_inplace


  logical function point_in_poly_lonlat(lon, lat, poly_lon, poly_lat) result(inside)
    ! Ray casting in lon-lat plane.
    real(real32), intent(in) :: lon, lat
    real(real32), intent(in) :: poly_lon(:), poly_lat(:)
    integer(int32) :: n, i, j
    real(real32) :: xi, yi, xj, yj
    logical :: c

    n = size(poly_lon)
    if (n < 3) then
      inside = .false.
      return
    end if
    if (size(poly_lat) /= n) then
      inside = .false.
      return
    end if

    c = .false.
    j = n
    do i=1, n
      xi = poly_lon(i); yi = poly_lat(i)
      xj = poly_lon(j); yj = poly_lat(j)

      if ( ((yi > lat) .neqv. (yj > lat)) ) then
        if ( lon < (xj - xi) * (lat - yi) / (yj - yi + 1.0e-30_real32) + xi ) then
          c = .not. c
        end if
      end if
      j = i
    end do
    inside = c
  end function point_in_poly_lonlat


  subroutine gauss_legendre(n, x, w)
    ! Nodes/weights on [-1,1], standard Newton iteration.
    ! Uses real64 internally for robustness, returns real32 arrays.
    integer(int32), intent(in) :: n
    real(real32), intent(out)  :: x(:), w(:)

    integer(int32) :: i, j, m
    real(real64) :: z, z1, p1, p2, p3, pp
    real(real64), parameter :: pi = 3.1415926535897932384626433832795_real64
    real(real64) :: xi, wi

    if (size(x) /= n .or. size(w) /= n) stop "gauss_legendre: bad array sizes"
    if (n < 1) stop "gauss_legendre: n must be >= 1"

    m = (n + 1) / 2
    do i=1, m
      z = cos(pi * (real(i,real64) - 0.25_real64) / (real(n,real64) + 0.5_real64))
      do
        p1 = 1.0_real64
        p2 = 0.0_real64
        do j=1, n
          p3 = p2
          p2 = p1
          p1 = ((2.0_real64*real(j,real64)-1.0_real64)*z*p2 - (real(j,real64)-1.0_real64)*p3) / real(j,real64)
        end do
        pp = real(n,real64) * (z*p1 - p2) / (z*z - 1.0_real64)
        z1 = z
        z  = z1 - p1/pp
        if (abs(z - z1) < 1.0e-14_real64) exit
      end do

      xi = z
      wi = 2.0_real64 / ((1.0_real64 - z*z) * pp*pp)

      x(i)       = real(-xi, real32)
      x(n+1-i)   = real( xi, real32)
      w(i)       = real( wi, real32)
      w(n+1-i)   = real( wi, real32)
    end do
  end subroutine gauss_legendre

end module gravmag_sphere_subs