module gravmag_sphere_state
  use, intrinsic :: iso_fortran_env, only: real32, int32
  implicit none
  private
  public :: ensure_alloc, clear_fields
  public :: bx, by, bz

  !***********************************************************************
  ! gravmag_sphere_state
  !
  ! lightweight field-array state for bx/by/bz work grids
  !
  ! old method note:
  ! - legacy implementations relied heavily on large COMMON blocks
  !
  ! current method:
  ! - allocatable module arrays + explicit resize/clear helpers
  ! - safer memory behavior and easier reuse when grid dimensions change
  !***********************************************************************

  real(real32), allocatable :: bx(:,:), by(:,:), bz(:,:)

contains

  ! ensure bx/by/bz are allocated to exactly match the requested grid
  ! if dimensions changed between bodies, arrays are reallocated
  subroutine ensure_alloc(nlat, nlon)
    integer(int32), intent(in) :: nlat, nlon
    if (.not. allocated(bx)) then
      allocate(bx(nlat,nlon), by(nlat,nlon), bz(nlat,nlon))
    else
      if (size(bx,1) /= nlat .or. size(bx,2) /= nlon) then
        deallocate(bx,by,bz)
        allocate(bx(nlat,nlon), by(nlat,nlon), bz(nlat,nlon))
      end if
    end if
  end subroutine ensure_alloc

  ! explicit clear helper so callers do not rely on stale state
  subroutine clear_fields()
    if (allocated(bx)) bx = 0.0_real32
    if (allocated(by)) by = 0.0_real32
    if (allocated(bz)) bz = 0.0_real32
  end subroutine clear_fields

end module gravmag_sphere_state
