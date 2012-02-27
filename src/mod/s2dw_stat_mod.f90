!------------------------------------------------------------------------------
! s2dw_stat_mod  -- S2DW library stat class
! 
!! Functionality to compute statistics of wavelet coefficients.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 February 2012
!
! Revisions:
!   February 2012 - Written by Jason McEwen
!------------------------------------------------------------------------------

module s2dw_stat_mod

  use s2dw_types_mod
  use s2dw_error_mod
  use s2dw_core_mod
  use omp_lib

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
       s2dw_stat_moments, &
       s2dw_stat_moments_write, &
       s2dw_stat_histogram, &
       s2dw_stat_histogram_write


  !----------------------------------------------------------------------------

contains


  !--------------------------------------------------------------------------
  ! s2dw_stat_moments
  !
  !! Compute moments of wavelets coefficients for each scale,
  !! including the mean, variance, skewness and kurtosis.
  !!
  !! Notes:
  !!  - Memory for the computed moments must already be allocated by the 
  !!    calling routine.
  !!
  !! Variables:
  !!  - wavdyn(0:J)%coeff: Wavelet coefficients for each scale [input].
  !!  - J: Maximum analysis scale depth [input].
  !!  - B: Harmonic band limit [input].
  !!  - N: Azimuthal band limit [input].
  !!  - bl_scoeff: Upper band limit for scaling coefficients [input].
  !!  - alpha: Basis dilation factor [input].
  !!  - mean(0:J): Mean of wavelet coefficients for each scale.
  !!  - var(0:J): Variance of wavelet coefficients for each scale.
  !!  - skew(0:J): Skewness of wavelet coefficients for each scale.
  !!  - kur(0:J): Kurtosis of wavelet coefficients for each scale.
  !
  !! @author J. D. McEwen
  !! @version 0.1 February 2012
  !
  ! Revisions:
  !   February 2012 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  subroutine s2dw_stat_moments(wavdyn, J, B, N, alpha, &
       mean, var, skew, kur)

    type(s2dw_wav_abg), intent(in) :: wavdyn(0:J)
    integer, intent(in) :: J
    integer, intent(in) :: B
    integer, intent(in) :: N
    real(dp), intent(in) :: alpha
    real(dp), intent(out) :: mean(0:J)
    real(dp), intent(out) :: var(0:J)
    real(dp), intent(out) :: skew(0:J)
    real(dp), intent(out) :: kur(0:J)

    real(dp), allocatable :: s(:,:,:), p(:,:,:)
    real(dp) :: ep, sdev
    integer :: jj, bl_hi, bl_lo, nj
    integer :: fail = 0

    !$omp parallel default(none) &
    !$omp shared(jj, wavdyn, B, N, J, alpha, mean, var, skew, kur) &
    !$omp private(bl_hi, bl_lo, nj, s, p, ep, sdev, fail)
    !$omp do schedule(dynamic,1) 
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

       ! Set number of samples for given analysis depth jj.
       nj = (2*bl_hi-1) * (2*bl_hi) * N        
       if (nj < 2) then
          call s2dw_error(S2DW_ERROR_SIZE_INVALID, 's2dw_stat_moments', &
               comment_add='Require at least two samples to compute moments')
       end if

       ! Allocate temporary memory.
       allocate(s(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1), stat=fail)
       allocate(p(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1), stat=fail)
       if(fail /= 0) then
          call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_stat_moments')
       end if

       ! Compute moments.

       mean(jj) = &
            sum(wavdyn(jj)%coeff(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1)) / nj

       s(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) = &
            wavdyn(jj)%coeff(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) - mean(jj)
       ep = sum(s(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1))

       p(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) = &
            s(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) &
            * s(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1)
       var(jj) = sum(p(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1))

       p(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) = &
            p(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) &
            * s(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1)
       skew(jj) = sum(p(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1))

       p(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) = &
            p(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) &
            * s(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1)
       kur(jj) = sum(p(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1))

       var(jj) = (var(jj) - ep**2/real(nj,dp)) / real(nj-1,dp)
       sdev = sqrt(var(jj))

       if (var(jj) < TOL_ZERO) then
          call s2dw_error(S2DW_ERROR_ARTH_WARNING, 's2dw_stat_moments', &
               'Skewness/kurtosis undefined when variance zero')
          skew(jj) = 0d0
          kur(jj) = 0d0
       else
          skew(jj) = skew(jj) / (nj * sdev**3)
          kur(jj) = kur(jj) / (nj * var(jj)**2) - 3d0
       end if

       ! Free temporary memory.
       deallocate(s, p)

    end do
    !$omp end do
    !$omp end parallel

  end subroutine s2dw_stat_moments


  !--------------------------------------------------------------------------
  ! s2dw_stat_moments_write
  !
  !! Write moments of wavelets coefficients for each scale, including
  !! the mean, variance, skewness and kurtosis.  Moments are written
  !! to a single line (ordered mean, variance, skewness and kurtosis)
  !! for each set of wavelet coefficients.
  !!
  !! Variables:
  !!  - filename: Name of file to write wavelet coefficient moments to 
  !!    [input].
  !!  - nsim: Number of sets of wavelet coefficient moments, i.e. number 
  !!    of skies/maps analysed [input].
  !!  - description(0:nsim-1): Description of each set of wavelet
  !!    coefficients stats, typically a filename [input].
  !!  - J: Maximum analysis scale depth [input].
  !!  - mean(0:nsim-1,0:J): Mean of wavelet coefficients for each map
  !!    and scale [input].
  !!  - var(0:nsim-1,0:J): Variance of wavelet coefficients for each
  !!    map and scale [input].
  !!  - skew(0:nsim-1,0:J): Skewness of wavelet coefficients for each
  !!    map and scale [input].
  !!  - kur(0:nsim-1,0:J): Kurtosis of wavelet coefficients for each
  !!    map and scale [input].
  !
  !! @author J. D. McEwen
  !! @version 0.1 February 2012
  !
  ! Revisions:
  !   February 2012 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  subroutine s2dw_stat_moments_write(filename, nsim, description, &
       J, mu, var, skew, kur)

    character(len=*), intent(in) :: filename
    integer, intent(in) :: nsim
    character(len=STRING_LEN), intent(in) :: description(0:nsim-1)
    integer, intent(in) :: J
    real(dp), intent(in) :: mu(0:nsim-1,0:J)
    real(dp), intent(in) :: var(0:nsim-1,0:J)
    real(dp), intent(in) :: skew(0:nsim-1,0:J)
    real(dp), intent(in) :: kur(0:nsim-1,0:J)

    integer :: fileid = 10, isim, jj

    ! Open file for writing.
    open(fileid, file=filename, form='formatted', &
         status='replace', action='write')

    ! Write statistics.
    do isim = 0, nsim-1
       write(fileid, '(a,a)', advance='no') trim(description(isim)), ', '
       write(fileid, '(i4,a)', advance='no') J, ', '
       do jj = 0, J
          write(fileid, '(e27.20,a)', advance='no')  mu(isim,jj), ', '
       end do
       do jj = 0, J
          write(fileid, '(e27.20,a)', advance='no')  var(isim,jj), ', '
       end do
       do jj = 0, J
          write(fileid, '(e27.20,a)', advance='no')  skew(isim,jj), ', '
       end do
       do jj = 0, J
          write(fileid, '(e27.20,a)', advance='no')  kur(isim,jj), ', '
       end do
       write(fileid,*)
    end do

    ! Close file.
    close(fileid)

  end subroutine s2dw_stat_moments_write


  !--------------------------------------------------------------------------
  ! s2dw_stat_histogram
  !
  !! Compute histogram of wavelets coefficients for each scale.
  !!
  !! Notes:
  !!  - Memory for the computed histograms must already be allocated by the 
  !!    calling routine.
  !!
  !! Variables:
  !!  - wavdyn(0:J)%coeff: Wavelet coefficients for each scale [input].
  !!  - J: Maximum analysis scale depth [input].
  !!  - B: Harmonic band limit [input].
  !!  - N: Azimuthal band limit [input].
  !!  - alpha: Basis dilation factor [input].
  !!  - hist_nbins: Number of histogram bins [input].
  !!  - hist_bins(0:J,0:hist_nbins-1): Central position of each bin.
  !!  - hist_vals(0:J,0:hist_nbins-1): Count in each bin.
  !
  !! @author S. M. Feeney
  !! @version 0.1 February 2012
  !
  ! Revisions:
  !   February 2012 - Written by Stephen Feeney
  !--------------------------------------------------------------------------

  subroutine s2dw_stat_histogram(wavdyn, J, B, N, alpha, hist_nbins, &
                                 hist_bins, hist_vals)

    type(s2dw_wav_abg), intent(in) :: wavdyn(0:J)
    integer, intent(in) :: J
    integer, intent(in) :: B
    integer, intent(in) :: N
    real(dp), intent(in) :: alpha
    integer, intent(in) :: hist_nbins
    real(dp), intent(out) :: hist_bins(0:J,0:hist_nbins-1)
    integer, intent(out) :: hist_vals(0:J,0:hist_nbins-1)

    real(dp) :: hist_bin_size
    integer :: ii, jj, kk, bl_hi, bl_lo, nj
    integer :: fail = 0
    
    
    !$omp parallel default(none) &
    !$omp shared(jj, wavdyn, B, N, J, alpha, hist_nbins, hist_bins, hist_vals) &
    !$omp private(bl_hi, bl_lo, nj, ii, fail, hist_bin_size)
    !$omp do schedule(dynamic,1) 
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

       ! Set number of samples for given analysis depth jj.
       nj = (2*bl_hi-1) * (2*bl_hi) * N
       
       ! Determine and validate the number of bins.
       if (hist_nbins .lt. 2) then
          call s2dw_error(S2DW_ERROR_SIZE_INVALID, 's2dw_stat_histogram', &
               comment_add='Require at least two bins to compute histogram')
       end if
       if (hist_nbins .ne. size(hist_vals, dim = 2)) then
          call s2dw_error(S2DW_ERROR_SIZE_INVALID, 's2dw_stat_histogram', &
               comment_add='Histogram bin and count array sizes do not match')
       end if
       
       ! Determine the bin size.
       hist_bin_size = hist_bins(jj,1) - hist_bins(jj,0)
       
       ! Compute histogram. Consider first bin separately to ensure the
       ! full range of data values are counted.
       hist_vals(jj,0) = count(wavdyn(jj)%coeff(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) .ge. &
                               hist_bins(jj,0) - 0.5d0 * hist_bin_size .and. &
                               wavdyn(jj)%coeff(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) .le. &
                               hist_bins(jj,1) + 0.5d0 * hist_bin_size)
       do ii = 1, hist_nbins - 1
          hist_vals(jj,ii) = count(wavdyn(jj)%coeff(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) .gt. &
                                   hist_bins(jj,ii) - 0.5d0 * hist_bin_size .and. &
                                   wavdyn(jj)%coeff(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) .le. &
                                   hist_bins(jj,ii+1) + 0.5d0 * hist_bin_size)
       end do
       
    end do
    !$omp end do
    !$omp end parallel

  end subroutine s2dw_stat_histogram


  !--------------------------------------------------------------------------
  ! s2dw_stat_histogram_write
  !
  !! Write histograms of wavelets coefficients for each scale. The first line
  !! of output contains the j = 0 bins, the second line contains the j = 0
  !! bin counts and so on.
  !!
  !! Variables:
  !!  - filename: Name of file to write wavelet coefficient moments to 
  !!    [input].
  !!  - nsim: Number of sets of wavelet coefficient moments, i.e. number 
  !!    of skies/maps analysed [input].
  !!  - J: Maximum analysis scale depth [input].
  !!  - hist_nbins: Number of histogram bins [input].
  !!  - hist_bins(0:J,0:hist_nbins-1): Central position of each bin [input].
  !!  - hist_vals(0:J,0:hist_nbins-1): Count in each bin [input].
  !
  !! @author S. M. Feeney
  !! @version 0.1 February 2012
  !
  ! Revisions:
  !   February 2012 - Written by Stephen Feeney
  !--------------------------------------------------------------------------

  subroutine s2dw_stat_histogram_write(filename, J, hist_nbins, hist_bins, &
                                       hist_vals)
    
    character(len=*), intent(in) :: filename
    integer, intent(in) :: J, hist_nbins
    real(dp), intent(in) :: hist_bins(0:J,0:hist_nbins-1)
    integer, intent(in) :: hist_vals(0:J,0:hist_nbins-1)

    integer :: fileid = 10, ii, jj
    character(len=8) :: repeats

    ! Open file for writing.
    open(fileid, file=filename, form='formatted', &
         status='replace', action='write')

    ! Write histograms.
    write(repeats, '(I8.1)') hist_nbins
    do jj = 0, J
       write(fileid, '(' // repeats // '(E19.12, X))'), &
             hist_bins(jj,:)
       write(fileid, '(' // repeats // '(I10.1, X))'), &
             hist_vals(jj,:)
    end do

    ! Close file.
    close(fileid)

  end subroutine s2dw_stat_histogram_write


end module s2dw_stat_mod
