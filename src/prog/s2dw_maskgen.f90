!------------------------------------------------------------------------------
! s2dw_maskgen
!
!! Computes an extended mask in wavelet space from a sky mask (defined
!! on a Healpix grid) and the parameters describing a wavelet transform.
!! 
!! Notes:
!!   - The extended mask is constructed for the wavelet coefficients
!!     only and not the scaling coefficients.
!!   - The extended mask is constructed by taking the wavelet
!!     coefficients of the original mask and masking regions with
!!     coefficients greater than 100*thres_factor percent of the
!!     maximum wavelet coefficent, for each scale.  Regions with
!!     wavelet coefficients below this value are not masked.  The
!!     original mask (in the wavelet domain) is then reapplied to
!!     create the final extended mask.
!!
!! Usage: s2dw_maskgen
!!   - [-help]: Displays usage information.
!!   - [-inp filename_in]: Name of input Healpix mask fits map.
!!   - [-out filename_out]: Name of output S2DW formatted fits/matlab file 
!!     containing extended mask.
!!   - [-file_type file_type (fits; m)]: String specifying type of output 
!!     S2DW file to write  (fits or matlab m) [default=fits].
!!   - [-thres_factor thres_factor]: Threshold factor to apply to
!!     maximum wavelet coefficient when constructing extended mask
!!     [default=0.1].
!!   - [-B B]: Harmonic band limit (if not specified set to 2*nside).
!!   - [-N N]: Azimuthal band limit [default=3].
!!   - [-alpha alpha]: Basis dilation factor [default=2].
!!   - [-J J]: Maximum analysis scale (optional) [default=Jmax].
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2007
!
! Revisions:
!   November 2007 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2dw_maskgen

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
  integer :: J_max
  integer :: B
  logical :: B_in = .false.
  integer :: N
  integer :: bl_scoeff
  real(dp) :: alpha
  integer :: nside
  logical :: admiss_pass
  type(s2_sky) :: sky
  integer :: fail = 0
  logical :: use_Jmax
  character(len=STRING_LEN) :: filename_in, filename_out
  integer :: file_extension = 1

  character(len=*), parameter ::  FILE_TYPE_FITS = 'fits'
  character(len=*), parameter ::  FILE_TYPE_MAT = 'm'
  character(len=STRING_LEN) :: file_type = FILE_TYPE_FITS
  character(len=STRING_LEN) :: error_string

  real(dp), allocatable :: max_abs_coeff(:)
  integer :: jj, bl_lo, bl_hi, aa, bb
  real(dp) :: thres, thres_factor
  real(dp), allocatable :: xtp(:,:)

  ! Set default parameter values.
  N = 3
  alpha = 2d0
  use_Jmax = .true.
  filename_in = 'wmap_ilc_3yr_v2_n64.fits'
  filename_out = 'wmap_ilc_3yr_v2_n64.s2dw'
  thres_factor = 0.10

  ! Parse options from command line.
  call parse_options()

  ! Read sky.
  sky = s2_sky_init(filename_in, file_extension)
  nside = s2_sky_get_nside(sky)
  if (.not. B_in) B = 2*nside
  J_max = ceiling(log(real(B,dp))/log(alpha) - TOL_CEIL)
  if(use_Jmax) J = J_max
  if(J > J_max) then
     J = J_max
     write(error_string,'(a,i4)') 'J too large, setting to maximum J = ', J_max
     call s2dw_error(S2DW_ERROR_ARG_WARNING, 's2dw_maskgen', &
          comment_add=trim(error_string))
  end if

  ! Allocate memory.
  allocate(flm_temp(0:B-1,0:B-1), stat=fail)
  allocate(flm(0:B-1,0:B-1), stat=fail)
  allocate(K_gamma(0:J,0:B-1), stat=fail)
  allocate(Phi2(0:B-1), stat=fail)
  allocate(Slm(0:B-1,0:N-1), stat=fail)
  allocate(admiss(0:B-1), stat=fail)
  if(fail /= 0) then
     call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_maskgen')
  end if

  ! Compute spherical harmonic coefficients.
  call s2_sky_compute_alm(sky, B-1, B-1)
  call s2_sky_get_alm(sky, flm_temp(0:B-1,0:B-1))
  flm(0:B-1,0:B-1) = flm_temp(0:B-1,0:B-1)

  ! Compute kernels, scaling function and directionality coefficients.
  call s2dw_core_init_kernels(K_gamma, Phi2, bl_scoeff, J, B, alpha)
  call s2dw_core_init_directionality(Slm, B, N)
  admiss_pass = s2dw_core_admiss(admiss, K_gamma, Phi2, B, J)
  if(.not. admiss_pass) then
     call s2dw_error(S2DW_ERROR_ADMISS_FAIL, 's2dw_maskgen')
  end if

  ! Allocate memory for scaling coefficients (cannot be done earlier 
  ! since don't know bl_scoeff).
  allocate(scoeff(0:bl_scoeff-1, 0:bl_scoeff-1), stat=fail)
  if(fail /= 0) then
     call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_maskgen')
  end if

  ! Perform S2DW analysis.
  call s2dw_core_analysis_flm2wav_dynamic(wavdyn, scoeff, flm, K_gamma, Slm, &
       J, B, N, bl_scoeff, alpha)

  ! Find maximum absolute value of wavelet coefficients on each scale.
  allocate(max_abs_coeff(0:J), stat=fail)
  if(fail /= 0) then
     call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_maskgen')
  end if
  do jj = 0, J
     max_abs_coeff(jj) = maxval(abs(wavdyn(jj)%coeff))
  end do

  ! Compute extended mask.
  do jj = 0,J

     ! Set band limits.
     if(s2dw_core_assimilate_int(B / (alpha**(jj-1)), bl_hi)) then
        bl_hi = min(bl_hi , B)
     else
        bl_hi = min(ceiling(B / (alpha**(jj-1)) ) , B)
     end if
     if(s2dw_core_assimilate_int(B / (alpha**(jj+1)), bl_lo)) then
        bl_lo = max(bl_lo, 0)
     else
        bl_lo = max(floor(B / (alpha**(jj+1)) ), 0)
     end if

     ! Set threshold.
     thres = thres_factor * max_abs_coeff(jj)

     ! Determine extended regions of mask by thresholding wavelet
     ! coefficients.
     where(abs(wavdyn(jj)%coeff(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1)) >= thres)
        wavdyn(jj)%coeff(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) = 0d0
     elsewhere
        wavdyn(jj)%coeff(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) = 1d0
     end where   

     ! Convert original mask to wavelet space (i.e. from Healpix to
     ! equiangular representation of wavelet space) and apply.
     allocate(xtp(0:2*bl_hi-1,0:2*bl_hi-2), stat=fail)
     call s2_sky_extract_ab_s2dw(sky, xtp, bl_hi)
     if(fail /= 0) then
        call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_maskgen')
     end if
     do aa = 0,2*bl_hi-2
        do bb = 0,2*bl_hi-1
           wavdyn(jj)%coeff(aa, bb, 0:N-1) = & 
                wavdyn(jj)%coeff(aa, bb, 0:N-1) &
                * xtp(bb,aa)
        end do
     end do
     deallocate(xtp)

  end do   

  ! Save wavelet and scaling coefficients.
  select case (trim(file_type))

  case (FILE_TYPE_FITS)
     call s2dw_fileio_fits_wav_write(wavdyn, scoeff, J, B, N, bl_scoeff, &
          alpha, filename_out)

  case (FILE_TYPE_MAT)
     call s2dw_fileio_matlab_wav_write(wavdyn, scoeff, J, B, N, &
          bl_scoeff, alpha, filename_out)

  case default
     call s2dw_error(S2DW_ERROR_FILEIO, 's2dw_analysis', &
          comment_add='Invalid file type option')

  end select

  ! Free memory.
  deallocate(flm_temp, flm, K_gamma, Phi2, Slm, admiss, scoeff)
  call s2dw_core_free_wavdyn(wavdyn)
  call s2_sky_free(sky)
  deallocate(max_abs_coeff)


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

    integer :: nn, i
    character(len=STRING_LEN) :: opt
    character(len=STRING_LEN) :: arg

    nn = nArguments()

    do i = 1,nn,2

       call getArgument(i,opt)

       if (i == nn .and. trim(opt) /= '-help') then
          write(*,*) 'Error: Option ', trim(opt), ' has no argument'
          stop
       end if

       if(trim(opt) /= '-help') call getArgument(i+1,arg)

       ! Read each argument in turn
       select case (trim(opt))

       case ('-help')
          write(*,'(a)') 'Usage: s2dw_maskgen [-inp filename_in]'
          write(*,'(a)') '                    [-out filename_out]'
          write(*,'(a)') '                    [-file_type file_type (fits; m)]'
          write(*,'(a)') '                    [-thres_factor thres_factor]'  
          write(*,'(a)') '                    [-B B]'  
          write(*,'(a)') '                    [-N N]'  
          write(*,'(a)') '                    [-alpha alpha]' 
          write(*,'(a)') '                    [-J J]'
          stop

       case ('-inp')
          filename_in = trim(arg)

       case ('-out')
          filename_out = trim(arg)

       case ('-file_type')
          file_type = trim(arg)

       case ('-thres_factor')
          read(arg,*) thres_factor

       case ('-B')
          read(arg,*) B
          B_in = .true.

       case ('-N')
          read(arg,*) N

       case ('-alpha')
          read(arg,*) alpha

       case ('-J')
          read(arg,*) J
          use_Jmax = .false.

       case default
          print '("unknown option ",a," ignored")', trim(opt)

       end select
    end do

  end subroutine parse_options


end program s2dw_maskgen
