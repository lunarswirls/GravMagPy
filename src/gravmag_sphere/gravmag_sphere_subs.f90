module gravmag_sphere_subs
  use, intrinsic :: iso_fortran_env, only: real32, int32
  implicit none
  private
  public :: wrap180_scalar, unwrap_lons_inplace, point_in_poly_lonlat

  !***********************************************************************
  ! gravmag_sphere_subs
  !
  ! geometry helpers shared by the solver
  !
  ! old method note:
  ! - legacy sphere code organized polygon boundaries using sorted
  !     theta/phi tables and multiple bookkeeping subroutines
  !
  ! current method:
  ! - polygon handling is done directly in lon-lat with:
  !     wrap/unwrap longitude support and a ray-casting inclusion test
  ! - this is simpler for fixed-limit and arbitrary multi-vertex bodies
  !***********************************************************************

contains

  ! canonical wrap to [-180, 180)
  subroutine wrap180_scalar(lon_deg)
    real(real32), intent(inout) :: lon_deg
    lon_deg = modulo(lon_deg + 180.0_real32, 360.0_real32) - 180.0_real32
  end subroutine wrap180_scalar

  subroutine unwrap_lons_inplace(lon_deg, lon_ref_deg)
    ! unwrap longitudes to be continuous near lon_ref_deg
    real(real32), intent(inout) :: lon_deg(:)
    real(real32), intent(in)    :: lon_ref_deg
    integer(int32) :: i
    real(real32) :: x, dx

    do i=1, size(lon_deg)
      x = lon_deg(i)
      call wrap180_scalar(x)
      dx = x - lon_ref_deg
      if (dx >  180.0_real32) x = x - 360.0_real32
      if (dx <= -180.0_real32) x = x + 360.0_real32
      lon_deg(i) = x
    end do
  end subroutine unwrap_lons_inplace

  logical function point_in_poly_lonlat(lon, lat, poly_lon, poly_lat)
    ! Ray casting in lon-lat plane, polygon assumed unwrapped already
    real(real32), intent(in) :: lon, lat
    real(real32), intent(in) :: poly_lon(:), poly_lat(:)

    integer(int32) :: i, j, n
    real(real32) :: xi, yi, xj, yj
    logical :: inside
    logical :: intersect

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

      ! crossing test for one polygon edge (xi,yi)->(xj,yj)
      intersect = ((yi > lat) .neqv. (yj > lat)) .and. &
                  (lon < (xj - xi) * (lat - yi) / max((yj - yi), 1.0e-30_real32) + xi)
      if (intersect) inside = .not. inside
      j = i
    end do

    point_in_poly_lonlat = inside
  end function point_in_poly_lonlat

end module gravmag_sphere_subs
