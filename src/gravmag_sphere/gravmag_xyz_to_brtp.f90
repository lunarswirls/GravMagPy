program gravmag_xyz_to_brtp
  use, intrinsic :: iso_fortran_env, only: int32, real64
  implicit none

  !***********************************************************************
  ! gravmag_xyz_to_brtp
  !
  ! converts cartesian field components (x,y,z) on a lon-lat grid to
  ! spherical components (br,btheta,bphi) using the standard basis:
  !   - r: outward
  !   - theta: southward (increasing colatitude)
  !   - phi: eastward (increasing longitude)
  !
  ! input expected per row (with or without body_id):
  !   body_id lon lat Fx Fy Fz Ftot
  ! or
  !   lon lat Fx Fy Fz Ftot
  !
  ! output always writes:
  !   body_id lon lat Br Btheta Bphi Btot
  !
  ! units are preserved from the input components.
  !
  ! old-vs-new note:
  ! - legacy sphere workflows selected components inside the forward
  !     solver and wrote one component per run (messy)
  ! - this post-processing keeps the forward solver in cartesian XYZ
  !     and performs coordinate conversion as a separate deterministic step
  !***********************************************************************

  integer, parameter :: wp = real64
  real(wp), parameter :: pi = 3.1415926535897932384626433832795_wp
  real(wp), parameter :: deg2rad = pi/180.0_wp

  character(len=512) :: infile, outfile, line
  character(len=32)  :: unit_lbl
  integer(int32) :: argc, ios
  integer(int32) :: inunit, outunit
  integer(int32) :: lineno, nread, nwrite, nskip

  integer(int32) :: body_id
  real(wp) :: lon_deg, lat_deg
  real(wp) :: fx, fy, fz, ftot_in
  real(wp) :: br, btheta, bphi, ftot_out
  real(wp) :: theta, phi, st, ct, sp, cp
  logical :: wrote_header

  argc = command_argument_count()
  if (argc < 2) then
    write(*,*) 'Usage: gravmag_xyz_to_brtp <input_xyz.txt> <output_brtp.txt>'
    stop 2
  end if

  call get_command_argument(1, infile)
  call get_command_argument(2, outfile)

  inunit = 21
  outunit = 22

  open(inunit, file=trim(infile), status='old', action='read', iostat=ios)
  if (ios /= 0) stop 'ERROR: cannot open input file'

  open(outunit, file=trim(outfile), status='replace', action='write', iostat=ios)
  if (ios /= 0) stop 'ERROR: cannot open output file'

  unit_lbl = ''
  wrote_header = .false.
  lineno = 0
  nread = 0
  nwrite = 0
  nskip = 0

  do
    ! read one raw line at a time so we can preserve flexible row formats
    read(inunit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    lineno = lineno + 1

    if (len_trim(line) == 0) cycle

    ! pass over comments, but mine unit hints from the source header
    if (line(1:1) == '#') then
      if (index(line, '_nT') > 0) then
        unit_lbl = 'nT'
      else if (index(line, '_mGal') > 0 .or. index(line, '_MGAL') > 0) then
        unit_lbl = 'mGal'
      end if
      cycle
    end if

    ! try parsing 7-column format with body_id first
    read(line, *, iostat=ios) body_id, lon_deg, lat_deg, fx, fy, fz, ftot_in
    if (ios /= 0) then
      ! fallback: 6-column format without body_id
      body_id = 0
      read(line, *, iostat=ios) lon_deg, lat_deg, fx, fy, fz, ftot_in
      if (ios /= 0) then
        nskip = nskip + 1
        cycle
      end if
    end if

    if (.not. wrote_header) then
      ! write one self-describing header once, then stream data rows
      write(outunit,'(A)') '# converted from cartesian xyz to spherical br,btheta,bphi'
      if (trim(unit_lbl) == 'nT') then
        write(outunit,'(A)') '# body_id lon_deg[-180,180) lat_deg Br_nT Btheta_nT Bphi_nT Btot_nT'
      else if (trim(unit_lbl) == 'mGal') then
        write(outunit,'(A)') '# body_id lon_deg[-180,180) lat_deg Br_mGal Btheta_mGal Bphi_mGal Btot_mGal'
      else
        write(outunit,'(A)') '# body_id lon_deg[-180,180) lat_deg Br Btheta Bphi Btot'
      end if
      wrote_header = .true.
    end if

    phi = lon_deg * deg2rad
    theta = (90.0_wp - lat_deg) * deg2rad

    st = sin(theta)
    ct = cos(theta)
    sp = sin(phi)
    cp = cos(phi)

    ! dot product with local spherical basis:
    !   e_r     = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
    !   e_theta = (cos(theta)cos(phi), cos(theta)sin(phi),-sin(theta))
    !   e_phi   = (-sin(phi),         cos(phi),            0)
    br     =  fx*st*cp + fy*st*sp + fz*ct
    btheta =  fx*ct*cp + fy*ct*sp - fz*st
    bphi   = -fx*sp    + fy*cp

    ftot_out = sqrt(br*br + btheta*btheta + bphi*bphi)

    write(outunit,'(I8,1X,F12.6,1X,F12.6,1X,ES15.7,1X,ES15.7,1X,ES15.7,1X,ES15.7)') &
      body_id, lon_deg, lat_deg, br, btheta, bphi, ftot_out

    nread = nread + 1
    nwrite = nwrite + 1
  end do

  close(inunit)
  close(outunit)

  write(*,'(A,1X,A)') 'Input :', trim(infile)
  write(*,'(A,1X,A)') 'Output:', trim(outfile)
  write(*,'(A,I0)')   'Rows converted: ', nwrite
  if (nskip > 0) write(*,'(A,I0)') 'Rows skipped (parse failures): ', nskip

end program gravmag_xyz_to_brtp
