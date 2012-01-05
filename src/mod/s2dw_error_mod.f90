!------------------------------------------------------------------------------
! s2dw_error_mod  -- S2DW library error class
! 
!! Functionality to handle errors that may occur in the s2dw library.  Public
!! s2dw error codes are defined, with corresponding private error comments and 
!! default halt execution status.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 October 2007
!
! Revisions:
!   October 2007 - Written by Jason McEwen
!------------------------------------------------------------------------------

module s2dw_error_mod

  use s2dw_types_mod, only: STRING_LEN

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: s2dw_error


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  integer, parameter :: S2DW_ERROR_NUM = 13

  integer, public, parameter :: &
       S2DW_ERROR_NONE = 0, &
       S2DW_ERROR_INIT = 1, &
       S2DW_ERROR_NOT_INIT = 2, &
       S2DW_ERROR_INIT_FAIL = 3, &
       S2DW_ERROR_MEM_ALLOC_FAIL = 4, &
       S2DW_ERROR_ARTH = 5, &
       S2DW_ERROR_SIZE_WARNING = 6, &
       S2DW_ERROR_SIZE_INVALID = 7, &
       S2DW_ERROR_SIZE_NOT_DEF = 8, &
       S2DW_ERROR_ARG_INVALID = 9, &
       S2DW_ERROR_ARG_WARNING = 10, &
       S2DW_ERROR_QUAD_STEP_EXCEED = 11, &
       S2DW_ERROR_ADMISS_FAIL = 12, &
       S2DW_ERROR_FILEIO = 13

  ! Each element of the error_comment array must have the same length, thus
  ! space with trailing space characters.  When come to use trim to remove 
  ! trailing spaces.
  !! Comment associated with each error type.
  character(len=STRING_LEN), parameter :: &
       error_comment(S2DW_ERROR_NUM+1) = &
       (/ & 
       'No error                                                                 ', &
       'Attempt to initialise object that has already been initialised           ', &
       'Object not initialised                                                   ', &
       'Object initialisation failed                                             ', &
       'Memory allocation failed                                                 ', &
       'Arithmetic exception                                                     ', &
       'Warning: Sizes not in recommended range                                  ', &
       'Invalid sizes                                                            ', &
       'Sizes not defined                                                        ', &
       'Arguments invalid                                                        ', &
       'Argument warning                                                         ', &
       'Exceeded number of steps limit whem computing quadrature                 ', &
       'Admissibility test failed (resolution of identity not satisfied)         ', &
       'File IO error                                                            ' &
       /) 

  !! Default program halt status of each error type.
  logical, parameter :: &
       halt_default(S2DW_ERROR_NUM+1) = &
       (/ &
       .false., &
       .true.,  &
       .true.,  &
       .true.,  &
       .true.,  &
       .true.,  &
       .false., &
       .true.,  &
       .true.,  &
       .true.,  &	
       .false.,  &	
       .true.,  &	
       .true.,  &	
       .true.  /)


  !----------------------------------------------------------------------------

contains


  !--------------------------------------------------------------------------
  ! s2dw_error
  !
  !! Displays error message corresponding to error_code and halt program 
  !! execution if required.
  !!
  !! Variables:
  !!   - error_code: Integer error code.
  !!   - [procedure]: Procedure name where s2dw_error called from.  Displayed 
  !!     when error message printed to screen.
  !!   - [comment_add]: If present, additional comment to append to default 
  !!     error comment.
  !!   - [comment_out]: If present the error comment is copied to comment_out
  !!     on output.
  !!   - [halt_in]: If present overrides default halt value.
  !
  !! @author J. D. McEwen
  !! @version 0.1 August 2004
  !
  ! Revisions:
  !   August 2004 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  subroutine s2dw_error(error_code, procedure, comment_add, &
       comment_out, halt_in)

    integer, intent(in) :: error_code
    character(len=*), intent(in), optional :: procedure, comment_add
    character(len=*), intent(inout), optional :: comment_out
    logical, intent(in), optional :: halt_in

    logical :: halt
    character(len=*), parameter :: comment_prefix = 'S2DW_ERROR: '

    !---------------------------------------
    ! Display error message
    !---------------------------------------

    if(present(procedure)) then

       if(present(comment_add)) then
          write(*,'(a,a,a,a,a,a,a,a)') comment_prefix, 'Error ''', &
               trim(error_comment(error_code+1)), &
               ''' occured in procedure ''', &
               trim(procedure), &
               '''', &
               ' - ', trim(comment_add)
       else
          write(*,'(a,a,a,a,a,a)') comment_prefix, 'Error ''', &
               trim(error_comment(error_code+1)), &
               ''' occured in procedure ''', &
               trim(procedure), &
               ''''
       end if

    else

       if(present(comment_add)) then
          write(*,'(a,a,a,a)') comment_prefix, &
               trim(error_comment(error_code+1)), &
               ' - ', trim(comment_add)
       else
          write(*,'(a,a)') comment_prefix, trim(error_comment(error_code+1))
       end if

    end if

    ! Copy error comment if comment_out present.
    if(present(comment_out)) comment_out = error_comment(error_code+1)

    !---------------------------------------
    ! Halt program execution if required
    !---------------------------------------

    if( present(halt_in) ) then
       halt = halt_in
    else
       halt = halt_default(error_code+1)
    end if

    if( halt ) then
       write(*,'(a,a,a,a,a)') comment_prefix, &
            '  Halting program execution ', &
            'due to error ''', trim(error_comment(error_code+1)), ''''
       stop
    end if

  end subroutine s2dw_error


end module s2dw_error_mod
