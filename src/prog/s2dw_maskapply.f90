!------------------------------------------------------------------------------
! s2dw_maskapply
!
!! Apply a mask (containing only ones and zeros) to a set of wavelet 
!! coefficients.  The output is the product of the mask and data coefficients.
!! If the display status is set then the masked output coefficients are 
!! overwritten with a magic number that appears grey when plotted. Output 
!! coefficients produced with the display option set should *only* be used for 
!! display purposes, and *not* for any subsequent analysis.
!!
!! Usage: s2dw_maskapply
!!   - [-help]: Display usage information.
!!   - [-data filename_data]: Name of input file containing wavelet 
!!     coefficients of data to mask.
!!   - [-mask filename_mask]: Name of input file containing wavelet 
!!     coefficients of mask to apply.
!!   - [-out filename_out]: Name of output file to write containing masked 
!!     wavelet coefficients.
!!   - [-file_type file_type (fits; m)]: String specifying type of output 
!!     S2DW file to write  (fits or matlab m) [default=fits].
!!   - [-display display]:  Logical specify whether to produce output
!!     coefficients for display purposes only, in which case the masked
!!     values are set to a magic number so that they appear grey when plotted.
!
!! @author S. M. Feeney (stephen.feeney.09@ucl.ac.uk)
!! @version 0.1 - February 2012
!
! Revisions:
!   February 2012 - Written by Stephen Feeney
!------------------------------------------------------------------------------

program s2dw_maskapply

  use s2dw_types_mod
  use s2dw_error_mod
  use s2dw_core_mod
  use s2dw_fileio_mod
  use healpix_types, only: HPX_DBADVAL

  implicit none

  real(dp), parameter :: FITS_DISPLAY_GREY_MAGIC_NUMBER = HPX_DBADVAL !-1.6375e30

  type(s2dw_wav_abg), allocatable :: wavdyn_data(:), wavdyn_mask(:)
  complex(dpc), allocatable :: scoeff_data(:,:), scoeff_mask(:,:)
  integer :: J, J_check, jj
  integer :: L, L_check
  integer :: N, N_check
  integer :: bl_scoeff
  real(dp) :: alpha, alpha_check

  character(STRING_LEN) :: filename_data
  character(STRING_LEN) :: filename_mask
  character(STRING_LEN) :: filename_out
  character(len=*), parameter ::  FILE_TYPE_FITS = 'fits'
  character(len=*), parameter ::  FILE_TYPE_MAT = 'm'
  character(len=STRING_LEN) :: file_type = FILE_TYPE_FITS
  integer :: S2DW_FITS_FILENAME_EXT_LEN = 4

  logical :: display = .true.

  ! Parse input parameters.
  call parse_options()
  
  ! Read wavelet coefficients of data and mask.
  select case (trim(file_type))

  case (FILE_TYPE_FITS)
     call s2dw_fileio_fits_wav_read(wavdyn_data, scoeff_data, &
                                    J, L, N, bl_scoeff, &
                                    alpha, filename_data)
     call s2dw_fileio_fits_wav_read(wavdyn_mask, scoeff_mask, &
                                    J_check, L_check, N_check, bl_scoeff, &
                                    alpha_check, filename_mask)
     S2DW_FITS_FILENAME_EXT_LEN = 5

  case (FILE_TYPE_MAT)
     call s2dw_fileio_matlab_wav_read(wavdyn_data, scoeff_data, &
                                      J, L, N, bl_scoeff, &
                                      alpha, filename_data)
     call s2dw_fileio_matlab_wav_read(wavdyn_mask, scoeff_mask, &
                                      J_check, L_check, N_check, bl_scoeff, &
                                      alpha_check, filename_mask)
     S2DW_FITS_FILENAME_EXT_LEN = 2

  case default
     call s2dw_error(S2DW_ERROR_FILEIO, 's2dw_maskapply', &
                     comment_add='Invalid file type option')

  end select
  
  ! Verify data and mask transform parameters are consistent.
  if (J /= J_check .or. L /= L_check .or.&
       N /= N_check .or. alpha /= alpha_check) then
     call s2dw_error(S2DW_ERROR_ARG_INVALID, 's2dw_maskapply', &
                     comment_add='Inconsistent wavelet transform parameters')
  end if
  
  ! Apply mask. If display set, then overwrite masked coefficients
  ! of product with magic number that is displayed as grey.
  if(display) then
     
     do jj = 0, J
        
        wavdyn_data(jj)%coeff = wavdyn_data(jj)%coeff * wavdyn_mask(jj)%coeff + &
                                (1.0d0 - wavdyn_mask(jj)%coeff) * &
                                FITS_DISPLAY_GREY_MAGIC_NUMBER
        
     end do
     
  else
     
     do jj = 0, J
        
        wavdyn_data(jj)%coeff = wavdyn_data(jj)%coeff * wavdyn_mask(jj)%coeff
        
     end do
     
  end if

  ! Save wavelet and scaling coefficients.
  select case (trim(file_type))

  case (FILE_TYPE_FITS)
     call s2dw_fileio_fits_wav_write(wavdyn_data, scoeff_data, &
                                     J, L, N, bl_scoeff, alpha, &
                                     filename_out)

  case (FILE_TYPE_MAT)
     call s2dw_fileio_matlab_wav_write(wavdyn_data, scoeff_data, &
                                       J, L, N, bl_scoeff, alpha, &
                                       filename_out)

  case default
     call s2dw_error(S2DW_ERROR_FILEIO, 's2dw_maskapply', &
                     comment_add='Invalid file type option')

  end select
  
  ! Free memory.
  call s2dw_core_free_wavdyn(wavdyn_data)
  call s2dw_core_free_wavdyn(wavdyn_mask)
  deallocate(scoeff_data, scoeff_mask)


 !----------------------------------------------------------------------------

  contains


    !---------------------------------------------------------------------
    ! parse_options
    !
    !! Parse the options passed when program called.
    !
    !! @author S. M. Feeney (stephen.feeney.09@ucl.ac.uk)
    !! @version 0.1 - February 2012
    !
    ! Revisions:
    !   February 2012 - Written by Stephen Feeney 
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
          write(*,'(a,a,a)') 'Option ', trim(opt), ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a)') 'Usage: s2dw_maskapply [-data filename_data]'
            write(*,'(a)') '                      [-mask filename_mask]'
            write(*,'(a)') '                      [-out filename_out]'
            write(*,'(a)') '                      [-file_type file_type (fits; m)]'
            write(*,'(a)') '                      [-display display]'
            stop
          
          case ('-data')
            filename_data = trim(arg)

          case ('-mask')
            filename_mask = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-display')
            read(arg,*) display

          case ('-file_type')
            file_type = trim(arg)

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2dw_maskapply
