!------------------------------------------------------------------------------
! s2dw_wav2stats
!
!! Compute wavelet coefficient statistics from already computed
!! wavelet coefficients.
!! 
!! Usage: s2dw_wav2stats
!!   - [-help]: Displays usage information.
!!   - [-inp filename_in]: Name of input S2DW formatted fits/matlab file 
!!     containing wavelet and scaling coefficients.
!!   - [-out filename_out]: Name of output statistic file.
!!   - [-file_type file_type (fits; m)]: String specifying type of input 
!!     S2DW file to read  (fits or matlab m file) [default=fits].
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   Febraury 2012 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2dw_wav2stats

  use s2dw_types_mod
  use s2dw_error_mod
  use s2dw_fileio_mod
  use s2dw_core_mod
  use s2dw_stat_mod
  use s2_sky_mod

  implicit none

  character(len=STRING_LEN) :: filename_in, filename_out
  character(len=*), parameter ::  FILE_TYPE_FITS = 'fits'
  character(len=*), parameter ::  FILE_TYPE_MAT = 'm'
  character(len=STRING_LEN) :: file_type = FILE_TYPE_FITS
  character(len=STRING_LEN) :: description(0:0)
  type(s2dw_wav_abg), allocatable :: wavdyn(:)
  complex(dpc), allocatable :: scoeff(:,:)
  integer :: J
  integer :: B
  integer :: N
  integer :: bl_scoeff
  real(dp) :: alpha
  integer :: fail = 0
  integer :: nsim
  real(dp), allocatable :: mu(:,:), var(:,:), skew(:,:), kur(:,:)

  ! Parse options from command line.
  call parse_options()

  ! Read wavelet and scaling coefficients from file.
  select case (trim(file_type))

  case (FILE_TYPE_FITS)
     call s2dw_fileio_fits_wav_read(wavdyn, scoeff, J, B, N, bl_scoeff, &
          alpha, filename_in)

  case (FILE_TYPE_MAT)
     call s2dw_fileio_matlab_wav_read(wavdyn, scoeff, J, B, N, &
          bl_scoeff, alpha, filename_in)

  case default
     call s2dw_error(S2DW_ERROR_FILEIO, 's2dw_wav2stats', &
          comment_add='Invalid file type option')

  end select

  ! Allocate memory for statistics.
  allocate(mu(0:0,0:J), stat=fail)
  allocate(var(0:0,0:J), stat=fail)
  allocate(skew(0:0,0:J), stat=fail)
  allocate(kur(0:0,0:J), stat=fail)
  if(fail /= 0) then
     call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_wav2stats')
  end if

  ! Compute statistics.
  call s2dw_stat_moments(wavdyn, J, B, N, alpha, &
       mu(0,0:J), var(0,0:J), skew(0,0:J), kur(0,0:J))

  ! Save statistics.
  nsim = 1
  description(0) = trim(filename_in)
  call s2dw_stat_moments_write(filename_out, nsim, description, &
       J, mu(0:0,0:J), var(0:0,0:J), skew(0:0,0:J), kur(0:0,0:J))

  ! Free memory.
  call s2dw_core_free_wavdyn(wavdyn)
  deallocate(scoeff)
  deallocate(mu, var, skew, kur)


  !----------------------------------------------------------------------------

contains


  !---------------------------------------------------------------------
  ! parse_options
  !
  !! Parses the options passed when program called.
  !
  !! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
  !! @version 0.1 - November 2007
  !
  ! Revisions:
  !   November 2007 - Written by Jason McEwen 
  !---------------------------------------------------------------------

  subroutine parse_options()

    use extension, only: getArgument, nArguments

    implicit none

    integer :: n, i
    character(len=STRING_LEN) :: opt
    character(len=STRING_LEN) :: arg

    n = nArguments()

    do i = 1,n,2

       call getArgument(i,opt)

       if (i == n .and. trim(opt) /= '-help') then
          write(*,*) 'Error: Option ', trim(opt), ' has no argument'
          stop
       end if

       if(trim(opt) /= '-help') call getArgument(i+1,arg)

       ! Read each argument in turn
       select case (trim(opt))

       case ('-help')
          write(*,'(a)') 'Usage: s2dw_wav2stats [-inp filename_in]'
          write(*,'(a)') '                      [-out filename_out]'
          write(*,'(a)') '                      [-file_type file_type (fits; m)]'
          stop

       case ('-inp')
          filename_in = trim(arg)

       case ('-out')
          filename_out = trim(arg)

       case ('-file_type')
          file_type = trim(arg)

       case default
          print '("unknown option ",a," ignored")', trim(opt)            

       end select
    end do

  end subroutine parse_options


end program s2dw_wav2stats
