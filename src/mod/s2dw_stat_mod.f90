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
       s2dw_stat_moments
     

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
          kur(jj) = -3d0
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


end module s2dw_stat_mod
