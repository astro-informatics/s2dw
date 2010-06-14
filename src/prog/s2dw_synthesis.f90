!------------------------------------------------------------------------------
! s2dw_synthesis
!
!! Reconstructs a Healpix sky map from S2DW wavelet and scaling coefficients.
!!
!! Notes:
!!   - The Healpix nside resolution parameter of the output sky map is set 
!!     by nside=B/2.
!!   - The spherical harmonic transform and inverse provided by Healpix 
!!     is not exact, hence these reduce the numerical precision of the 
!!     `exact' reconstruction of the original real space map from its wavelet 
!!     coefficients.
!! 
!! Usage: s2dw_synthesis
!!   - [-help]: Displays usage information.
!!   - [-inp filename_in]: Name of input S2DW formatted fits/matlab file 
!!     containing wavelet and scaling coefficients.
!!   - [-out filename_out]: Name of output Healpix sky fits map.
!!   - [-file_type file_type (fits; m)]: String specifying type of input 
!!     S2DW file to read  (fits or matlab m file) [default=fits].
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2007
!
! Revisions:
!   November 2007 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2dw_synthesis

	use s2dw_types_mod
	use s2dw_error_mod
	use s2dw_fileio_mod
	use s2dw_core_mod
	use s2_sky_mod

  implicit none

	complex(dpc), allocatable :: flm(:,:)
	complex(spc), allocatable :: flm_temp(:,:)
	real(dp), allocatable :: K_gamma(:,:)
	real(dp), allocatable :: Phi2(:)
	complex(dpc), allocatable :: Slm(:,:)
	real(dp), allocatable :: admiss(:)
  type(s2dw_wav_abg), allocatable :: wavdyn(:)
	complex(dpc), allocatable :: scoeff(:,:)
	integer :: J
	integer :: B
	integer :: N
	integer :: bl_scoeff
	real(dp) :: alpha
	integer :: nside
	logical :: admiss_pass
	type(s2_sky) :: sky
	integer :: fail = 0
	character(len=STRING_LEN) :: filename_in, filename_out

  character(len=*), parameter ::  FILE_TYPE_FITS = 'fits'
  character(len=*), parameter ::  FILE_TYPE_MAT = 'm'
  character(len=STRING_LEN) :: file_type = FILE_TYPE_FITS

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
       call s2dw_error(S2DW_ERROR_FILEIO, 's2dw_synthesis', &
            comment_add='Invalid file type option')

  end select

	! Allocate memory.
	allocate(flm(0:B-1,0:B-1), stat=fail)
	allocate(flm_temp(0:B-1,0:B-1), stat=fail)
	allocate(K_gamma(0:J,0:B-1), stat=fail)
	allocate(Phi2(0:B-1), stat=fail)
	allocate(Slm(0:B-1,0:N-1), stat=fail)
	allocate(admiss(0:B-1), stat=fail)
	if(fail /= 0) then
		call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_synthesis')
	end if

	! Compute kernels, scaling function and directionality coefficients.
	call s2dw_core_init_kernels(K_gamma, Phi2, bl_scoeff, J, B, alpha)
	call s2dw_core_init_directionality(Slm, B, N)
	admiss_pass = s2dw_core_admiss(admiss, K_gamma, Phi2, B, J)
	if(.not. admiss_pass) then
		call s2dw_error(S2DW_ERROR_ADMISS_FAIL, 's2dw_synthesis')
	end if

	! Perform S2DW synthesis.
	call s2dw_core_synthesis_wav2flm_dynamic(flm, wavdyn, scoeff, K_gamma, Phi2, &
		Slm, J, B, N, bl_scoeff, alpha)

	! Reconstruct real space sky from spherical harmonic coefficients.
	nside = B/2
	flm_temp(0:B-1,0:B-1) = flm(0:B-1,0:B-1)
	sky = s2_sky_init(flm_temp(0:B-1,0:B-1), B-1, B-1)
	call s2_sky_compute_map(sky, nside)

	! Save sky.
	call s2_sky_write_map_file(sky, filename_out)

	! Free memory.
	deallocate(flm, K_gamma, Phi2, Slm, admiss, scoeff)
	call s2_sky_free(sky)


	!----------------------------------------------------------------------------

  contains


    !---------------------------------------------------------------------
    ! parse_options
    !
    !! Parses the options passed when program called.
    !
    !!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
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
            write(*,'(a)') 'Usage: s2dw_synthesis [-inp filename_in]'
            write(*,'(a)') '                     [-out filename_out]'
            write(*,'(a)') '                     [-file_type file_type (fits; m)]'
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-file_type')
            file_type = trim(arg)

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program s2dw_synthesis
