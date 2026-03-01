module gravmag_sphere_physics
  use, intrinsic :: iso_fortran_env, only: real32, int32
  implicit none
  private
  public :: add_dipole_field_cart

contains

  subroutine add_dipole_field_cart(xo, yo, zo, xs, ys, zs, mx, my, mz, dBx, dBy, dBz)
    ! Dipole moment m = (mx,my,mz) in A*m^2
    ! Returns B contribution in Tesla.
    real(real32), intent(in)  :: xo, yo, zo      ! obs position (m)
    real(real32), intent(in)  :: xs, ys, zs      ! source position (m)
    real(real32), intent(in)  :: mx, my, mz      ! dipole moment (A*m^2)
    real(real32), intent(out) :: dBx, dBy, dBz   ! Tesla

    real(real32), parameter :: mu0_over_4pi = 1.0e-7_real32  ! (mu0/4pi) in SI
    real(real32) :: Rx, Ry, Rz
    real(real32) :: R2, invR, invR2, invR3, invR5
    real(real32) :: mdotr, c

    Rx = xo - xs
    Ry = yo - ys
    Rz = zo - zs

    R2 = Rx*Rx + Ry*Ry + Rz*Rz
    if (R2 <= 0.0_real32) then
      dBx = 0.0_real32; dBy = 0.0_real32; dBz = 0.0_real32
      return
    end if

    invR  = 1.0_real32 / sqrt(R2)
    invR2 = invR*invR
    invR3 = invR2*invR
    invR5 = invR3*invR2

    mdotr = mx*Rx + my*Ry + mz*Rz
    c     = 3.0_real32 * mdotr * invR5

    dBx = mu0_over_4pi * (c*Rx - mx*invR3)
    dBy = mu0_over_4pi * (c*Ry - my*invR3)
    dBz = mu0_over_4pi * (c*Rz - mz*invR3)
  end subroutine add_dipole_field_cart

end module gravmag_sphere_physics