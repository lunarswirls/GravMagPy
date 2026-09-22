module gravmag_paths
  implicit none
  private
  public :: default_output_path
contains
  function parent_directory(path) result(directory)
    character(len=*), intent(in) :: path
    character(len=:), allocatable :: directory
    integer :: slash
    slash = scan(trim(path), '/', back=.true.)
    directory = '.'
    if (slash == 1) directory = '/'
    if (slash > 1) directory = path(:slash-1)
  end function parent_directory

  function case_directory(path) result(directory)
    character(len=*), intent(in) :: path
    character(len=:), allocatable :: directory, current, name, parent
    integer :: slash
    directory = parent_directory(path)
    current = directory
    do
      slash = scan(current, '/', back=.true.)
      name = current(slash+1:)
      if (name == 'output' .or. name == 'figs') then
        directory = parent_directory(current)
        return
      end if
      parent = parent_directory(current)
      if (parent == current) exit
      current = parent
    end do
  end function case_directory

  function shell_quote(value) result(quoted)
    character(len=*), intent(in) :: value
    character(len=:), allocatable :: quoted
    integer :: index
    quoted = achar(39)
    do index = 1, len_trim(value)
      if (value(index:index) == achar(39)) then
        quoted = quoted // achar(39) // achar(34) // achar(39) // achar(34) // achar(39)
      else
        quoted = quoted // value(index:index)
      end if
    end do
    quoted = quoted // achar(39)
  end function shell_quote

  subroutine default_output_path(input_path, output_path, suffix, strip_xyz)
    character(len=*), intent(in) :: input_path, suffix
    character(len=*), intent(inout) :: output_path
    logical, intent(in), optional :: strip_xyz
    character(len=:), allocatable :: filename, directory, destination
    integer :: slash, dot, status, command_status
    if (len_trim(output_path) == 0) then
      slash = scan(trim(input_path), '/', back=.true.)
      filename = trim(input_path(slash+1:))
      dot = scan(filename, '.', back=.true.)
      if (dot > 1) filename = filename(:dot-1)
      if (present(strip_xyz)) then
        if (strip_xyz .and. len(filename) > 4) then
          if (filename(len(filename)-3:) == '_xyz') filename = filename(:len(filename)-4)
        end if
      end if
      directory = case_directory(input_path)
      destination = directory // '/output/' // filename // suffix
      if (len(destination) > len(output_path)) error stop 'ERROR: default output path is too long'
      output_path = destination
    end if
    if (trim(input_path) == trim(output_path)) error stop 'ERROR: input and output paths must differ'
    directory = parent_directory(trim(output_path))
    ! quote every path character so spaces and shell metacharacters remain literal
    call execute_command_line('mkdir -p -- ' // shell_quote(directory), exitstat=status, cmdstat=command_status)
    if (command_status /= 0) error stop 'ERROR: could not launch output directory creation'
    if (status /= 0) error stop 'ERROR: could not create output directory'
  end subroutine default_output_path
end module gravmag_paths
