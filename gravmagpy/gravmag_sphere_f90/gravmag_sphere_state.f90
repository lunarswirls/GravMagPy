module gravmag_sphere_state
  use, intrinsic :: iso_fortran_env, only: real32, int32
  implicit none
  private
  public :: ensure_alloc, clear_fields
  public :: bx, by, bz

  real(real32), allocatable :: bx(:,:), by(:,:), bz(:,:)

contains

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

  subroutine clear_fields()
    if (allocated(bx)) bx = 0.0_real32
    if (allocated(by)) by = 0.0_real32
    if (allocated(bz)) bz = 0.0_real32
  end subroutine clear_fields

end module gravmag_sphere_state