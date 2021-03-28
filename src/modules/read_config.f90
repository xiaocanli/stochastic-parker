!*******************************************************************************
!< Module for reading configuration file. The configuration file looks like
!< A = value or A: value. So this module defines a subroutine to read the value
!< with different delimiter.
!*******************************************************************************
module read_config
    implicit none
    private
    public get_variable, get_variable_int

    contains

    !---------------------------------------------------------------------------
    !< Function to the value of one variable with name var_name.
    !< Inputs:
    !<   fh: file handler.
    !<   var_name: the variable name.
    !<   delimiter: The value of the variable is after the delimiter.
    !< Returns:
    !<   var_value: the variable value.
    !---------------------------------------------------------------------------
    function get_variable(fh, var_name, delimiter) result(var_value)
        use constants, only: dp
        implicit none
        integer, intent(in) :: fh
        character(*), intent(in) :: var_name, delimiter
        real(dp) :: var_value
        character(len=150) :: single_line
        integer :: len1, IOstatus
        do while (index(single_line, var_name) == 0)
            read(fh, '(A)', IOSTAT=IOstatus) single_line
            if (IOStatus < 0) then
                exit
            endif
        enddo
        var_value = -1.0
        if (IOStatus == 0) then
            len1 = len(trim(single_line)) + 1
            ! For C or C++ code
            if (index(single_line, ";") /= 0) then
                len1 = index(single_line, ";")
            endif
            read(single_line(index(single_line, delimiter)+1:len1-1), *) var_value
        endif
    end function

    !---------------------------------------------------------------------------
    !< Get the variable as an integer.
    !---------------------------------------------------------------------------
    function get_variable_int(fh, var_name, delimiter) result(var_value_int)
        use constants, only: dp
        implicit none
        integer, intent(in) :: fh
        character(*), intent(in) :: var_name, delimiter
        real(dp) :: var_value
        integer :: var_value_int

        var_value = get_variable(fh, var_name, delimiter)
        var_value_int = int(var_value)
    end function get_variable_int

end module read_config
