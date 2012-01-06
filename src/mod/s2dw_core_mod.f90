!------------------------------------------------------------------------------
! s2dw_core_mod  -- S2DW library core class
! 
!! Functionality to perform scale discretised wavelet transform on the sphere.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 November 2007
!
! Revisions:
!   November 2007 - Written by Jason McEwen
!------------------------------------------------------------------------------

module s2dw_core_mod

  use s2dw_types_mod
  use s2dw_error_mod
  use s2dw_dl_mod
  use omp_lib

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
       s2dw_core_init_kernels, &
       s2dw_core_init_directionality, &
       s2dw_core_comp_Jmax, &
       s2dw_core_admiss, &
       s2dw_core_analysis_flm2wav, &
       s2dw_core_analysis_flm2wav_dynamic, &
       s2dw_core_synthesis_wav2flm, &
       s2dw_core_synthesis_wav2flm_dynamic, &
       s2dw_core_free_wavdyn, &
       s2dw_core_assimilate_int


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  integer, parameter :: FFTW_ESTIMATE=64, FFTW_MEASURE=0
  integer, parameter :: FFTW_FORWARD=-1, FFTW_BACKWARD=1


  !---------------------------------------
  ! Data types
  !---------------------------------------

  !! Structure to store wavelet coefficients for given scale.  Required for 
  !! dynamic memory allocation of wavelet coefficients.
  type, public :: s2dw_wav_abg
     real(dp), allocatable :: coeff(:,:,:)
  end type s2dw_wav_abg


  !----------------------------------------------------------------------------

contains


  !--------------------------------------------------------------------------
  ! Initialisation routine
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! s2dw_core_init_kernels
  !
  !! Compute kernels, scaling functions and band limits.
  !!
  !! Notes:
  !!  - Kernels are computed from differences of scaling functions, which 
  !!    are computed directly.  It is important to compute kernels in this
  !!    manner to ensure that any errors in the numerical integration 
  !!    collapse when computing the resolution of the identity.
  !!
  !! Variables:
  !!  - K_gamma(0:J, 0:B-1): Kernel functions computed for each scale
  !!    [output].
  !!  - Phi2(0:B-1): Scaling function squared at maximum j=J [output].
  !!  - bl_scoeff: Upper band limit for scaling coefficients [output].
  !!  - J: Maximum analysis scale depth [input].
  !!  - B: Harmonic band limit [input].
  !   - alpha: Basis dilation factor [input].
  !
  !! @author J. D. McEwen
  !! @version 0.1 October 2007
  !
  ! Revisions:
  !   October 2007 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  subroutine s2dw_core_init_kernels(K_gamma, Phi2, bl_scoeff, &
       J, B, alpha)

    integer, intent(in) :: J, B
    real(dp), intent(in) :: alpha
    real(dp), intent(out) :: K_gamma(0:J,0:B-1)
    real(dp), intent(out) :: Phi2(0:B-1)
    integer, intent(out) :: bl_scoeff

    real(dp) :: Phi2_all(0:J+1,0:B-1)
    real(dp) :: C_Phi
    real(dp) :: Phi2_diff
    integer :: bw
    integer :: jj, el, el_lo

    if(.not. s2dw_core_params_valid(J, B, alpha)) then
       call s2dw_error(S2DW_ERROR_SIZE_INVALID, 's2dw_core_init_kernels')
    end if

    do jj=0,J+1      
       if(.not. s2dw_core_assimilate_int(B / (alpha**(jj-1)), bw)) then
          bw = ceiling(B / (alpha**(jj-1)))  
       end if
       if(.not. s2dw_core_assimilate_int(real(bw,dp)/alpha, el_lo)) then
          el_lo = floor(real(bw,dp)/alpha)
       end if
       !C_Phi = qsimp(s2dw_core_phi2_integrand, bw, alpha, &
       !	(real(bw,dp)/alpha)+TOL_LIMIT, real(bw,dp)-TOL_LIMIT)
       C_Phi = qsimp(s2dw_core_phi2_integrand, bw, alpha, &
            real(el_lo,dp)+TOL_LIMIT, real(bw,dp)-TOL_LIMIT)
       do el=0,B-1								
          if(el<=el_lo) then
             Phi2_all(jj,el)=1
          elseif(el>=bw) then
             Phi2_all(jj,el)=0
          else				
             Phi2_all(jj,el) = &
                  qsimp(s2dw_core_phi2_integrand, bw, alpha, &
                  real(el,dp), real(bw,dp)-TOL_LIMIT) / C_Phi
          end if
          if(jj>=1) then
             Phi2_diff = Phi2_all(jj-1,el) - Phi2_all(jj,el)
             if(Phi2_diff < 0d0) then ! Due to numerical precision
                K_gamma(jj-1,el) = 0d0
             else
                K_gamma(jj-1,el) = sqrt(Phi2_diff)
             end if
          end if
       end do
    end do

    ! Save final scaling function.
    Phi2(0:B-1) = Phi2_all(J+1,0:B-1)

    ! Compute band limit for scaling coefficients and wavelets.
    bl_scoeff = bw
    !if(.not. s2dw_core_assimilate_int(real(B,dp) / (alpha**J), bl_scoeff) then
    !	bl_scoeff = ceiling(real(B,dp) / (alpha**J) )
    !end if

  end subroutine s2dw_core_init_kernels


  !--------------------------------------------------------------------------
  ! s2dw_core_phi2_integrand
  !
  !! Integrand of definite integral to be calculated to compute the scaling 
  !! function.
  !!
  !! Notes:
  !!  - Integrand is given by K_Phi^2(k')/k'.
  !!
  !! Variables:
  !!  - k: Dummy (scale) variable of integration [input].
  !!  - B: Harmonic band limit [input].
  !!  - alpha: Basis dilation factor [input].
  !!  - integrand: Value of integrand computed [output].
  !
  !! @author J. D. McEwen
  !! @version 0.1 October 2007
  !
  ! Revisions:
  !   October 2007 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  function s2dw_core_phi2_integrand(k, bw, alpha) result(integrand)

    real(dp), intent(in) :: k
    integer, intent(in) :: bw
    real(dp), intent(in) :: alpha
    real(dp) :: integrand

    real(dp) :: t
    integer :: bw_lo

    if(.not. s2dw_core_assimilate_int(real(bw,dp)/alpha, bw_lo)) then
       bw_lo = floor(real(bw,dp)/alpha)
    end if
    t = 2.0_dp * ( real(k-bw_lo,dp) / real(bw-bw_lo,dp) ) - 1
    !t = 2.0_dp * ( (alpha*k-bw) / ( (alpha-1)*bw) ) - 1
    integrand = exp( -2.0_dp / (1.0_dp-t**2) ) / k

  end function s2dw_core_phi2_integrand


  !--------------------------------------------------------------------------
  ! s2dw_core_init_directionality
  !
  !! Compute wavelet directionality coefficients.
  !!
  !! Variables:
  !!  - Slm(0:B-1, 0:N-1): Directionality coefficients [output].
  !!  - B: Harmonic band limit [input].
  !!  - N: Azimuthal band limit [input].
  !
  !! @author J. D. McEwen
  !! @version 0.1 October 2007
  !
  ! Revisions:
  !   October 2007 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  subroutine s2dw_core_init_directionality(Slm, B, N)

    integer, intent(in) :: B
    integer, intent(in) :: N
    complex(dpc), intent(out) :: Slm(0:B-1,0:N-1)

    integer :: el, m, betah, g, binom
    complex(dpc) :: eta			

    if(N>B) then
       call s2dw_error(S2DW_ERROR_SIZE_INVALID, 's2dw_core_init_directionality')
    end if

    Slm(0:B-1,0:N-1) = cmplx(0d0, 0d0)
    Slm(0,0) = cmplx(0d0, 0d0)
    do el = 1,B-1
       do m = 0,min(el,N-1)
          betah = (1 - (-1)**(N+m)) / 2;   
          if(betah > 0) then
             g = min(N-1, el - (1 + (-1)**(N+el))/2);						
             binom = binomial(g, (g-m)/2)
             if(mod(N-1,2) == 0) then
                eta = cmplx(1d0,0d0)
             else
                eta = cmplx(0d0,1d0)
             end if
             Slm(el,m) = sqrt(real(binom,dp)/(2d0**g)) * eta
          end if
       end do
    end do

  end subroutine s2dw_core_init_directionality


  !--------------------------------------------------------------------------
  ! s2dw_core_admiss
  !
  !! Compute resolution of the identity from scaled kernels and scaling 
  !! function and check equal to unity for all el.
  !!
  !! Notes:
  !!  - Resolution of identity must be unity for all el to ensure 
  !!    admissibility satisfied, i.e. to ensure exact reconstruction is 
  !!    possible.
  !!  - Resolution of identity computed by
  !!      Phi2(el) + sum_{jj=0}^{J} K_gamma^2(jj,el)     
  !!
  !! Variables:
  !!  - admiss(0:B-1): Resolution of the identity computed for each el
  !!    (all values should be close to one if passed admissibility test) 
  !!    [output].
  !!  - K_gamma(0:J, 0:B-1): Kernel functions computed for each scale
  !!    [input].
  !!  - Phi2(0:B-1): Scaling function squared at maximum j=J [input].
  !!  - J: Maximum analysis scale depth [input].
  !!  - B: Harmonic band limit [input].
  !!  - pass: Logical specifying whether kernels and scaling function 
  !!    passed admissibilty test [output].
  !
  !! @author J. D. McEwen
  !! @version 0.1 October 2007
  !
  ! Revisions:
  !   October 2007 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  function s2dw_core_admiss(admiss, K_gamma, Phi2, B, J) result(pass)

    integer, intent(in) :: J, B
    real(dp), intent(in) :: K_gamma(0:J,0:B-1)
    real(dp), intent(in) :: Phi2(0:B-1)
    real(dp), intent(out) :: admiss(0:B-1)
    logical :: pass

    integer :: el, jj

    pass = .true.
    do el=0,B-1
       admiss(el) = Phi2(el)
       do jj=0,J
          admiss(el) = admiss(el) + K_gamma(jj,el)**2
       end do
       if(abs(admiss(el) - 1d0) > TOL_ADMISS) pass = .false.					
    end do

  end function s2dw_core_admiss


  !--------------------------------------------------------------------------
  ! Analysis routines
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! s2dw_core_analysis_flm2wav
  !
  !! Compute wavelet and scaling coefficients from harmonic coefficients of 
  !! original signal using fast algorithm.
  !!
  !! Notes:
  !!  - Wavelet and scaling coefficients contain all information required to 
  !!    reconstruct the original signal exactly.
  !!
  !! Variables:
  !!  - wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1): Wavelet coefficients [output].
  !!  - scoeff(0:bl_scoeff-1, 0:bl_scoeff-1): Scaling coefficients 
  !!    [output].
  !!  - flm(0:B-1, 0:B-1): Harmonic coefficients of signal [input].
  !!  - K_gamma(0:J, 0:B-1): Kernel functions computed for each scale
  !!    [input].
  !!  - Slm(0:B-1, 0:N-1): Directionality coefficients [input].
  !!  - J: Maximum analysis scale depth [input].
  !!  - B: Harmonic band limit [input].
  !!  - N: Azimuthal band limit [input].
  !!  - bl_scoeff: Upper band limit for scaling coefficients [input].
  !!  - alpha: Basis dilation factor [input].
  !
  !! @author J. D. McEwen
  !! @version 0.1 October 2007
  !
  ! Revisions:
  !   October 2007 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  subroutine s2dw_core_analysis_flm2wav(wav, scoeff, flm, K_gamma, Slm, &
       J, B, N, bl_scoeff, alpha)

    integer, intent(in) :: J
    integer, intent(in) :: B
    integer, intent(in) :: N
    integer, intent(in) :: bl_scoeff
    real(dp), intent(in) :: alpha
    real(dp), intent(out) :: wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1)
    complex(dpc), intent(out) :: scoeff(0:bl_scoeff-1, 0:bl_scoeff-1)
    complex(dpc), intent(in) :: flm(0:B-1, 0:B-1)
    real(dp), intent(in) :: K_gamma(0:J, 0:B-1)
    complex(dpc), intent(in) :: Slm(0:B-1, 0:N-1)

    complex(dpc), allocatable :: Tmmm(:,:,:), Tmmg(:,:,:), Tmbg(:,:,:)
    real(dp), allocatable :: dl(:,:)
    complex(dpc) :: sb
    complex(dpc) :: eta
    integer :: el, m, mm, mmm, jj, bb, gg
    integer :: bl_hi, bl_lo, k_indicator, fail
    real(dp) :: msign, mmmsign
    real(dp) :: beta_b, gamma_g
    integer*8 :: fftw_plan
    complex(dpc) :: tmp(0:B-1)

    fail = 0

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

       ! Compute Tmmm.
       allocate(Tmmm(0:bl_hi-1, 0:bl_hi-1, 0:2*N-2), stat=fail)
       allocate(dl(-(bl_hi-1):bl_hi-1, -(bl_hi-1):bl_hi-1), stat=fail)
       if(fail /= 0) then
          call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_core_analysis_flm2wav')
       end if
       Tmmm(0:bl_hi-1, 0:bl_hi-1, 0:2*N-2) = cmplx(0d0, 0d0)
       !$omp parallel default(none) &
       !$omp shared(jj, el, K_gamma, Slm, flm, N, bl_lo, bl_hi) &
       !$omp private(dl,m,mm,mmm,msign,mmmsign,sb) &
       !$omp reduction(+: Tmmm)
       !$omp do schedule(dynamic, 10)
       do el = bl_lo+1,bl_hi-1

          ! For each l value create the plane of the d-matrix.
          call s2dw_dl_beta_operator(dl(-el:el,-el:el), PION2, el)

          if (mod(min(N-1,el),2) == 0) then
             ! Even case.
             mmmsign = -1d0
          else
             ! Odd case.
             mmmsign = 1d0
          end if

          do mmm = -min(N-1,el),min(N-1,el)
             mmmsign = -mmmsign

             do mm = 0,el

                msign = -1d0
                do m = 0,el
                   msign = -msign

                   if ((m < 0) .and. (mmm < 0)) then
                      sb = msign * mmmsign &
                           * K_gamma(jj,el)*Slm(el,-mmm) * conjg(flm(el, -m))
                   else if (m < 0) then
                      sb = msign &
                           * conjg(K_gamma(jj,el)*Slm(el,mmm) * flm(el, -m))
                   else if (mmm < 0) then
                      sb = mmmsign &
                           * K_gamma(jj,el)*Slm(el,-mmm) * flm(el, m)
                   else
                      sb = conjg(K_gamma(jj,el)*Slm(el,mmm)) * flm(el, m)
                   end if

                   Tmmm(m, mm, mmm+N-1) = Tmmm(m, mm, mmm+N-1) &
                        + exp(I*(mmm-m)*PION2) * dl(mm,m) * dl(mm,mmm) * sb

                end do
             end do
          end	do
       end do
       !$omp end do
       !$omp end parallel
       deallocate(dl)

       ! Compute Tmbg using fast SoV

       ! Compute Tmmg
       allocate(Tmmg(0:bl_hi-1, 0:bl_hi-1, 0:N-1), stat=fail)
       if(fail /= 0) then
          call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_core_analysis_flm2wav')
       end if
       Tmmg(0:bl_hi-1, 0:bl_hi-1, 0:N-1) = cmplx(0d0, 0d0)
       !$omp parallel default(none) &
       !$omp shared(Tmmm, gg, mm, bl_hi, N, Tmmg) &
       !$omp private(m, gamma_g, mmm, k_indicator)
       !$omp do schedule(dynamic,10) collapse(2)
       do gg = 0,N-1
          do mm = 0,bl_hi-1
             gamma_g = PI*gg/real(N,dp)
             do m = 0,bl_hi-1
                do mmm = -(N-1),N-1
                   k_indicator = (1-(-1)**(N+mmm))/2 
                   if(k_indicator > 0) then
                      Tmmg(m,mm,gg) = Tmmg(m,mm,gg) + Tmmm(m, mm, mmm+N-1) &
                           * exp( I * (mmm*gamma_g) )
                   end if
                end do
             end do
          end do
       end do
       !$omp end do
       !$omp end parallel
       deallocate(Tmmm)

       ! Compute Tmbg from Tmmg
       allocate(Tmbg(0:bl_hi-1, 0:2*bl_hi-1, 0:N-1), stat=fail)
       if(fail /= 0) then
          call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_core_analysis_flm2wav')
       end if
       Tmbg(0:bl_hi-1, 0:2*bl_hi-1, 0:N-1) = cmplx(0d0, 0d0)
       !$omp parallel default(none) &
       !$omp shared(Tmmg, gg, bb, bl_hi, N, Tmbg) &
       !$omp private(m, mm, beta_b, gamma_g, eta)
       !$omp do schedule(dynamic,10) collapse(2)
       do gg = 0,N-1
          do bb = 0,2*bl_hi-1
             beta_b = PI*(2d0*bb+1d0)/(4d0*bl_hi)					
             gamma_g = PI*gg/real(N,dp)				

             do mm = 0,bl_hi-1
                do m = 0,bl_hi-1

                   if (mm==0) then
                      eta = 1d0
                   else if(mod(m+N-1,2)==0) then
                      eta = 2d0*cos(mm*beta_b)
                   else
                      eta = I * 2d0*sin(mm*beta_b)
                   end if

                   Tmbg(m,bb,gg) = Tmbg(m,bb,gg) &
                        + Tmmg(m, mm, gg) * eta

                end do
             end do
          end do
       end do
       !$omp end do
       !$omp end parallel
       deallocate(Tmmg)

       ! Compute Tabg from Tmbg using FFT  
       wav(jj, 0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) = 0d0
       tmp(0:bl_hi-1) = cmplx(0d0, 0d0)
       call dfftw_plan_dft_c2r_1d(fftw_plan, 2*bl_hi-1, &
            tmp(0:bl_hi-1), wav(jj,0:2*bl_hi-2,0,0), &
            FFTW_MEASURE)
       !$omp parallel default(none) &
       !$omp shared(wav, jj, fftw_plan, Tmbg, gg, bb, bl_hi, N)
       !$omp do schedule(dynamic,10) collapse(2)
       do gg = 0,N-1
          do bb = 0,2*bl_hi-1
             call dfftw_execute_dft_c2r(fftw_plan, &
                  Tmbg(0:bl_hi-1,bb,gg), wav(jj,0:2*bl_hi-2,bb,gg))
          end do
       end do
       !$omp end do
       !$omp end parallel
       call dfftw_destroy_plan(fftw_plan)
       deallocate(Tmbg)

    end do

    ! Save scaling coefficients.
    scoeff(0:bl_scoeff-1, 0:bl_scoeff-1) = flm(0:bl_scoeff-1, 0:bl_scoeff-1)

  end subroutine s2dw_core_analysis_flm2wav


  !--------------------------------------------------------------------------
  ! s2dw_core_analysis_flm2wav_dynamic
  !
  !! Compute wavelet and scaling coefficients from harmonic coefficients of 
  !! original signal using fast algorithm, dynamically allocating the memory
  !! required for the wavelet coefficients dynamically.
  !!
  !! Notes:
  !!  - The memory required to store the wavelet coefficients is allocated
  !!    dynamically herein.
  !!  - Wavelet and scaling coefficients contain all information required to 
  !!    reconstruct the original signal exactly.
  !!
  !! Variables:
  !!  - wavdyn(0:J)%coeff: Wavelet coefficients for each scale (memory  
  !!    allocated herein) [output].
  !!  - scoeff(0:bl_scoeff-1, 0:bl_scoeff-1): Scaling coefficients 
  !!    [output].
  !!  - flm(0:B-1, 0:B-1): Harmonic coefficients of signal [input].
  !!  - K_gamma(0:J, 0:B-1): Kernel functions computed for each scale
  !!    [input].
  !!  - Slm(0:B-1, 0:N-1): Directionality coefficients [input].
  !!  - J: Maximum analysis scale depth [input].
  !!  - B: Harmonic band limit [input].
  !!  - N: Azimuthal band limit [input].
  !!  - bl_scoeff: Upper band limit for scaling coefficients [input].
  !!  - alpha: Basis dilation factor [input].
  !
  !! @author J. D. McEwen
  !! @version 0.1 February 2008
  !
  ! Revisions:
  !   October 2007 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  subroutine s2dw_core_analysis_flm2wav_dynamic(wavdyn, scoeff, flm, K_gamma, Slm, &
       J, B, N, bl_scoeff, alpha)

    integer, intent(in) :: J
    integer, intent(in) :: B
    integer, intent(in) :: N
    integer, intent(in) :: bl_scoeff
    real(dp), intent(in) :: alpha
    type(s2dw_wav_abg), intent(out), allocatable :: wavdyn(:)
    complex(dpc), intent(out) :: scoeff(0:bl_scoeff-1, 0:bl_scoeff-1)
    complex(dpc), intent(in) :: flm(0:B-1, 0:B-1)
    real(dp), intent(in) :: K_gamma(0:J, 0:B-1)
    complex(dpc), intent(in) :: Slm(0:B-1, 0:N-1)

    complex(dpc), allocatable :: Tmmm(:,:,:), Tmmg(:,:,:), Tmbg(:,:,:)
    real(dp), allocatable :: dl(:,:)
    real(dp) :: Ta(0:2*B-2)    ! Temp memory required to solve FFTW bug
    complex(dpc) :: sb
    complex(dpc) :: eta
    integer :: el, m, mm, mmm, jj, bb, gg
    integer :: bl_hi, bl_lo, k_indicator, fail
    real(dp) :: msign, mmmsign
    real(dp) :: beta_b, gamma_g
    integer*8 :: fftw_plan

    fail = 0

    allocate(wavdyn(0:J), stat=fail)
    if(fail /= 0) then
       call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_core_analysis_flm2wav_dynamic')
    end if

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

       ! Compute Tmmm.
       allocate(Tmmm(0:bl_hi-1, 0:bl_hi-1, 0:2*N-2), stat=fail)
       allocate(dl(-(bl_hi-1):bl_hi-1, -(bl_hi-1):bl_hi-1), stat=fail)
       if(fail /= 0) then
          call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_core_analysis_flm2wav_dynamic')
       end if
       Tmmm(0:bl_hi-1, 0:bl_hi-1, 0:2*N-2) = cmplx(0d0, 0d0)
       do el = bl_lo+1,bl_hi-1

          ! For each l value create the plane of the d-matrix.
          call s2dw_dl_beta_operator(dl(-el:el,-el:el), PION2, el)

          msign = -1d0
          do m = 0,el
             msign = -msign

             do mm = 0,el

                if (mod(min(N-1,el),2) == 0) then
                   ! Even case.
                   mmmsign = -1d0
                else
                   ! Odd case.
                   mmmsign = 1d0
                end if

                do mmm = -min(N-1,el),min(N-1,el)
                   mmmsign = -mmmsign

                   if ((m < 0) .and. (mmm < 0)) then
                      sb = msign * mmmsign &
                           * K_gamma(jj,el)*Slm(el,-mmm) * conjg(flm(el, -m))
                   else if (m < 0) then
                      sb = msign &
                           * conjg(K_gamma(jj,el)*Slm(el,mmm) * flm(el, -m))
                   else if (mmm < 0) then
                      sb = mmmsign &
                           * K_gamma(jj,el)*Slm(el,-mmm) * flm(el, m)
                   else
                      sb = conjg(K_gamma(jj,el)*Slm(el,mmm)) * flm(el, m)
                   end if

                   Tmmm(m, mm, mmm+N-1) = Tmmm(m, mm, mmm+N-1) &
                        + exp(I*(mmm-m)*PION2) * dl(mm,m) * dl(mm,mmm) * sb

                end do
             end do
          end	do
       end do
       deallocate(dl)

       ! Compute Tmbg using fast SoV

       ! Compute Tmmg
       allocate(Tmmg(0:bl_hi-1, 0:bl_hi-1, 0:N-1), stat=fail)
       if(fail /= 0) then
          call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_core_analysis_flm2wav_dynamic')
       end if
       Tmmg(0:bl_hi-1, 0:bl_hi-1, 0:N-1) = cmplx(0d0, 0d0)
       do m = 0,bl_hi-1
          do mm = 0,bl_hi-1
             do gg = 0,N-1
                gamma_g = PI*gg/real(N,dp)				

                do mmm = -(N-1),N-1
                   k_indicator = (1-(-1)**(N+mmm))/2 
                   if(k_indicator > 0) then
                      Tmmg(m,mm,gg) = Tmmg(m,mm,gg) + Tmmm(m, mm, mmm+N-1) &
                           * exp( I * (mmm*gamma_g) )
                   end if
                end do

             end do
          end do
       end do
       deallocate(Tmmm)

       ! Compute Tmbg from Tmmg
       allocate(Tmbg(0:bl_hi-1, 0:2*bl_hi-1, 0:N-1), stat=fail)
       if(fail /= 0) then
          call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_core_analysis_flm2wav_dynamic')
       end if
       Tmbg(0:bl_hi-1, 0:2*bl_hi-1, 0:N-1) = cmplx(0d0, 0d0)
       do bb = 0,2*bl_hi-1
          beta_b = PI*(2d0*bb+1d0)/(4d0*bl_hi)					
          do gg = 0,N-1
             gamma_g = PI*gg/real(N,dp)				

             do m=0,bl_hi-1

                Tmbg(m,bb,gg) = Tmbg(m,bb,gg) &
                     + Tmmg(m, 0, gg) 

                do mm =1,bl_hi-1

                   if(mod(m+N-1,2)==0) then
                      eta = 2d0*cos(mm*beta_b)
                   else
                      eta = I * 2d0*sin(mm*beta_b)
                   end if

                   Tmbg(m,bb,gg) = Tmbg(m,bb,gg) &
                        + Tmmg(m, mm, gg) * eta

                end do

             end do
          end do
       end do
       deallocate(Tmmg)

       ! Compute Tabg from Tmbg using FFT  
       allocate(wavdyn(jj)%coeff(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1), stat=fail)
       if(fail /= 0) then
          call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_core_analysis_flm2wav_dynamic')
       end if
       wavdyn(jj)%coeff(0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1) = 0d0

       do bb = 0,2*bl_hi-1
          do gg = 0,N-1

             call dfftw_plan_dft_c2r_1d(fftw_plan, 2*bl_hi-1, &
                  Tmbg(0:bl_hi-1,bb,gg), Ta(0:2*bl_hi-2), FFTW_ESTIMATE)
             call dfftw_execute(fftw_plan)
             call dfftw_destroy_plan(fftw_plan)

             wavdyn(jj)%coeff(0:2*bl_hi-2,bb,gg) = Ta(0:2*bl_hi-2)

             ! Note that dfftw_plan_dft_c2r_1d *must* be called with a one 
             ! dimensional real output array.  For example the following 
             ! produces invalid output:
             !   call dfftw_plan_dft_c2r_1d(fftw_plan, 2*bl_hi-1, &
             !     Tmbg(0:bl_hi-1,bb,gg), wav(jj,0:2*bl_hi-2,bb,gg), FFTW_ESTIMATE)
             ! This bug may be due to the way arrays are passed between 
             ! C and Fortran.

          end do
       end do
       deallocate(Tmbg)

    end do

    ! Save scaling coefficients.
    scoeff(0:bl_scoeff-1, 0:bl_scoeff-1) = flm(0:bl_scoeff-1, 0:bl_scoeff-1)

  end subroutine s2dw_core_analysis_flm2wav_dynamic


  !--------------------------------------------------------------------------
  ! Synthesis routines 
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! s2dw_core_synthesis_wav2flm
  !
  !! Synthesis harmonic coefficients of signal from wavelet and scaling
  !! coefficients.
  !!
  !! Notes:
  !!  - This routines uses less memory than using s2dw_core_synthesis_wav2wig
  !!    followed by s2dw_core_synthesis_wig2flm since the Wigner coefficients
  !!    do not need to be stored for all j.
  !!
  !! Variables:
  !!  - flm(0:B-1, 0:B-1): Synthesised harmonic coefficients of signal 
  !!    [output].
  !!  - wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1): Wavelet coefficients [input].
  !!  - scoeff(0:bl_scoeff-1, 0:bl_scoeff-1): Scaling coefficients 
  !!    [input].
  !!  - K_gamma(0:J, 0:B-1): Kernel functions computed for each scale
  !!    [input].
  !!  - Phi2(0:B-1): Scaling function squared at maximum j=J [input].
  !!  - Slm(0:B-1, 0:N-1): Directionality coefficients [input].
  !!  - J: Maximum analysis scale depth [input].
  !!  - B: Harmonic band limit [input].
  !!  - N: Azimuthal band limit [input].
  !!  - bl_scoeff: Upper band limit for scaling coefficients [input].
  !!  - alpha: Basis dilation factor [input].
  !
  !! @author J. D. McEwen
  !! @version 0.1 October 2007
  !
  ! Revisions:
  !   October 2007 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  subroutine s2dw_core_synthesis_wav2flm(flm, wav, scoeff, K_gamma, Phi2, &
       Slm, J, B, N, bl_scoeff, alpha)

    integer, intent(in) :: J
    integer, intent(in) :: B
    integer, intent(in) :: N
    integer, intent(in) :: bl_scoeff
    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1)
    complex(dpc), intent(in) :: scoeff(0:bl_scoeff-1, 0:bl_scoeff-1)
    complex(dpc), intent(out) :: flm(0:B-1, 0:B-1)
    real(dp), intent(in) :: K_gamma(0:J, 0:B-1)
    real(dp), intent(in) :: Phi2(0:B-1)
    complex(dpc), intent(in) :: Slm(0:B-1, 0:N-1)

    integer :: k_indicator, bb, gg, jj, m, mm, mmm, el
    integer :: bl_hi, bl_lo
    integer :: fail
    real(dp) :: beta_b, gamma_g
    real(dp) :: w
    complex(dpc), allocatable :: wig(:,:,:), V(:,:,:)
    !complex(dpc), allocatable :: Ummp(:,:,:)
    complex(dpc), allocatable :: Umnb(:,:,:)

    real(dp) :: Wa(0:2*B-2) 
    real(dp), allocatable :: dl(:,:)
    integer*8 :: fftw_plan

    fail = 0
    flm(0:B-1, 0:B-1) = cmplx(0d0, 0d0)

    do jj=0,J

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

       ! Step 1: Integrate over alpha using FFT 
       ! (i.e. compute Umbp but store in V)
       allocate(V(0:bl_hi-1, 0:2*bl_hi-1, 0:N-1), stat=fail)
       if(fail /= 0) then
          call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, &
               's2dw_core_synthesis_wav2flm')
       end if
       V(0:bl_hi-1, 0:2*bl_hi-1, 0:N-1) = cmplx(0d0,0d0)
       call dfftw_plan_dft_r2c_1d(fftw_plan, 2*bl_hi-1, &
            Wa(0:2*bl_hi-2), V(0:bl_hi-1,0,0), FFTW_MEASURE)
       !$omp parallel default(none) &
       !$omp shared(gg, bb, fftw_plan, wav, jj, V, bl_hi, N)
       !$omp do schedule(dynamic,10) collapse(2)
       do gg = 0,N-1
          do bb = 0,2*bl_hi-1
             call dfftw_execute_dft_r2c(fftw_plan, &
                  wav(jj,0:2*bl_hi-2,bb,gg), V(0:bl_hi-1,bb,gg))
             V(0:bl_hi-1,bb,gg) = V(0:bl_hi-1,bb,gg) &
                  * (4d0*PI**2/real(N,dp)) / (2d0*bl_hi-1d0) 
          end do
       end do
       !$omp end do
       !$omp end parallel
       call dfftw_destroy_plan(fftw_plan)
       ! Could delete wavelet coefficients here.

       ! Step 2a: Integerate over gamma (i.e. compute Umnb)	
       allocate(Umnb(0:bl_hi-1, 0:N-1, 0:2*bl_hi-1), stat=fail)
       if(fail /= 0) then
          call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_core_synthesis_wav2flm')
       end if
       Umnb(0:bl_hi-1, 0:N-1, 0:2*bl_hi-1) = cmplx(0d0, 0d0)
       !$omp parallel default(none) &
       !$omp shared(gg, mmm, bl_hi, N, Umnb, V) &
       !$omp private(gamma_g, k_indicator)
       !$omp do schedule(dynamic,1) collapse(2)
       do gg = 0,N-1
          do mmm = 0,N-1  	
             gamma_g = PI*gg/real(N,dp)
             k_indicator = (1-(-1)**(N+mmm))/2 
             if(k_indicator > 0) then
                Umnb(0:bl_hi-1,mmm,0:2*bl_hi-1) = &
                     Umnb(0:bl_hi-1,mmm,0:2*bl_hi-1) &
                     + V(0:bl_hi-1,0:2*bl_hi-1,gg) * exp(-I*mmm*gamma_g)
             end if
          end do
       end do
       !$omp end do
       !$omp end parallel

       ! Step 2b: Integrate over beta	(i.e. compute Vmmm)
       V(0:bl_hi-1, 0:2*bl_hi-2, 0:N-1) = cmplx(0d0, 0d0)
       !$omp parallel default(none) &
       !$omp shared(mmm, bb, bl_hi, N, Umnb) &
       !$omp private(beta_b, w, k_indicator, mm) &
       !$omp reduction(+: V)
       !$omp do schedule(dynamic,10) collapse(2)
       do bb = 0,2*bl_hi-1
          do mmm = 0,N-1
             beta_b = PI*(2d0*bb+1d0)/(4d0*bl_hi)
             w = quad_weights(beta_b, bl_hi)
             k_indicator = (1-(-1)**(N+mmm))/2 
             if(k_indicator > 0) then
                do mm = -(bl_hi-1),bl_hi-1  
                   V(0:bl_hi-1,mm+bl_hi-1,mmm) = V(0:bl_hi-1,mm+bl_hi-1,mmm) &
                        + Umnb(0:bl_hi-1, mmm, bb) * w * exp(-I*mm*beta_b)
                end do
             end if
          end do
       end do
       !$omp end do
       !$omp end parallel
       deallocate(Umnb)

       ! Step 3a: Combine Vmmm terms for m' and -m' outside of el sum
       !$omp parallel default(none) &
       !$omp shared(mmm, m, bl_hi, N, V) &
       !$omp private(k_indicator, mm)
       !$omp do schedule(dynamic,10) collapse(2)       
       do mmm = 0,N-1
          do mm = 1,bl_hi-1
             k_indicator = (1-(-1)**(N+mmm))/2	
             if(k_indicator > 0) then
                do m = 0,bl_hi-1
                   V(m, mm+bl_hi-1, mmm) = V(m, mm+bl_hi-1, mmm) &
                        + (-1)**(m+mmm) *  V(m, -mm+bl_hi-1, mmm) 
                end do
             end if
          end do
       end do
       !$omp end do
       !$omp end parallel

       ! Step 3b: Compute Wigner coefficients
       allocate(wig(0:bl_hi-1, 0:bl_hi-1, 0:N-1), stat=fail)
       allocate(dl(-(bl_hi-1):bl_hi-1, -(bl_hi-1):bl_hi-1), stat=fail)
       if(fail /= 0) then
          call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, &
               's2dw_core_synthesis_wav2flm')
       end if
       wig(0:bl_hi-1, 0:bl_hi-1, 0:N-1) = cmplx(0d0, 0d0)
       !$omp parallel default(none) &
       !$omp shared(el, bl_lo, bl_hi, N, V, wig) &
       !$omp private(dl, k_indicator, mmm, m)
       !$omp do schedule(dynamic,10)
       do el = bl_lo+1,bl_hi-1
          call s2dw_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
          do mmm = 0,min(N-1,el)			
             do m = 0,el
                k_indicator = (1-(-1)**(N+mmm))/2 	
                if(k_indicator > 0) then
                   wig(el,m,mmm) = sum(V(m, bl_hi-1:el+bl_hi-1, mmm)  &
                        * dl(0:el,m) * dl(0:el,mmm)) * exp(-I*(mmm-m)*PION2)
                end if
             end do
          end do
       end do
       !$omp end do
       !$omp end parallel
       deallocate(dl, V)

       ! Compute harmonic component (for this jj) from Wigner coefficients
       !$omp parallel default(none) &
       !$omp shared(el, bl_lo, bl_hi, N, wig, K_gamma, Slm, jj, flm) &
       !$omp private(m, mmm, k_indicator)
       !$omp do schedule(dynamic,10)
       do el = bl_lo+1,bl_hi-1
          do mmm = 0,min(N-1,el)
             do m = 0,el
                if (mmm == 0) then
                   flm(el,m) = flm(el,m) + &
                        (real(2d0*el+1,dp)/(8d0*pi**2)) &
                        * wig(el,m,0) * K_gamma(jj,el) * Slm(el,0)
                else
                   k_indicator = (1-(-1)**(N+mmm))/2
                   if(k_indicator > 0) then
                      flm(el,m) = flm(el,m) + &
                           2d0 * (real(2d0*el+1,dp)/(8*pi**2)) &
                           * wig(el,m,mmm) * K_gamma(jj,el) * Slm(el,mmm) 
                   end if
                end if
             end do
          end do
       end do
       !$omp end do
       !$omp end parallel
       deallocate(wig)

    end do

    ! Add in scaling part of signal
    do el = 0,bl_scoeff-1
       flm(el,0:bl_scoeff-1) = flm(el,0:bl_scoeff-1) &
            + Phi2(el) * scoeff(el,0:bl_scoeff-1)
    end do

  end subroutine s2dw_core_synthesis_wav2flm


  !--------------------------------------------------------------------------
  ! s2dw_core_synthesis_wav2flm_dynamic
  !
  !! Synthesise harmonic coefficients of signal from wavelet and scaling
  !! coefficients, where the memory of the wavelet coefficients has been 
  !! allocated dynamically (destroying wavelet coefficients in the process).
  !!
  !! Notes:
  !!  - The wavelet coefficients are destroyed by this routine and their 
  !!    memory is freed.
  !!  - This routines uses less memory than using s2dw_core_synthesis_wav2wig
  !!    followed by s2dw_core_synthesis_wig2flm since the Wigner coefficients
  !!    do not need to be stored for all j.
  !!
  !! Variables:
  !!  - flm(0:B-1, 0:B-1): Synthesised harmonic coefficients of signal 
  !!    [output].
  !!  - wavdyn(0:J)%coeff: Wavelet coefficients for each scale (destroyed 
  !!    and freed on output) [input/output].
  !!  - scoeff(0:bl_scoeff-1, 0:bl_scoeff-1): Scaling coefficients 
  !!    [input].
  !!  - K_gamma(0:J, 0:B-1): Kernel functions computed for each scale
  !!    [input].
  !!  - Phi2(0:B-1): Scaling function squared at maximum j=J [input].
  !!  - Slm(0:B-1, 0:N-1): Directionality coefficients [input].
  !!  - J: Maximum analysis scale depth [input].
  !!  - B: Harmonic band limit [input].
  !!  - N: Azimuthal band limit [input].
  !!  - bl_scoeff: Upper band limit for scaling coefficients [input].
  !!  - alpha: Basis dilation factor [input].
  !
  !! @author J. D. McEwen
  !! @version 0.1 February 2008
  !
  ! Revisions:
  !   February 2008 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  subroutine s2dw_core_synthesis_wav2flm_dynamic(flm, wavdyn, scoeff, K_gamma, Phi2, &
       Slm, J, B, N, bl_scoeff, alpha)

    integer, intent(in) :: J
    integer, intent(in) :: B
    integer, intent(in) :: N
    integer, intent(in) :: bl_scoeff
    real(dp), intent(in) :: alpha
    type(s2dw_wav_abg), intent(inout), allocatable :: wavdyn(:)
    complex(dpc), intent(in) :: scoeff(0:bl_scoeff-1, 0:bl_scoeff-1)
    complex(dpc), intent(out) :: flm(0:B-1, 0:B-1)
    real(dp), intent(in) :: K_gamma(0:J, 0:B-1)
    real(dp), intent(in) :: Phi2(0:B-1)
    complex(dpc), intent(in) :: Slm(0:B-1, 0:N-1)

    integer :: k_indicator, bb, gg, jj, m, mm, mmm, el
    integer :: bl_hi, bl_lo
    integer :: fail
    real(dp) :: beta_b, gamma_g
    real(dp) :: w
    complex(dpc), allocatable :: wig(:,:,:), V(:,:,:)
    complex(dpc), allocatable :: Umnb(:,:,:)

    complex(dpc) :: Ua(0:B-1)  ! Temp memory required to solve FFTW bug
    real(dp) :: Wa(0:2*B-2)    ! Temp memory required to solve FFTW bug
    real(dp), allocatable :: dl(:,:)
    integer*8 :: fftw_plan

    fail = 0
    flm(0:B-1, 0:B-1) = cmplx(0d0, 0d0)

    do jj=0,J

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

       ! Step 1: Integrate over alpha using FFT 
       ! (i.e. compute Umbp but store in V)
       allocate(V(0:bl_hi-1, 0:2*bl_hi-1, 0:N-1), stat=fail)
       if(fail /= 0) then
          call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, &
               's2dw_core_synthesis_wav2flmv_dynamic')
       end if
       V(0:bl_hi-1, 0:2*bl_hi-1, 0:N-1) = cmplx(0d0,0d0)
       do bb = 0,2*bl_hi-1
          do gg = 0,N-1
             !Wa(0:2*bl_hi-2) = wav(jj,0:2*bl_hi-2,bb,gg)
             Wa(0:2*bl_hi-2) = wavdyn(jj)%coeff(0:2*bl_hi-2,bb,gg)
             call dfftw_plan_dft_r2c_1d(fftw_plan, 2*bl_hi-1, &
                  Wa(0:2*bl_hi-2), Ua(0:bl_hi-1), FFTW_ESTIMATE)
             call dfftw_execute(fftw_plan)
             call dfftw_destroy_plan(fftw_plan)
             V(0:bl_hi-1,bb,gg) = Ua(0:bl_hi-1) &
                  * (4d0*PI**2/real(N,dp)) / (2d0*bl_hi-1d0) 
          end do
       end do
       ! Free memory for wavelet coefficients for this scale.
       deallocate(wavdyn(jj)%coeff)

       ! Step 2a: Integerate over gamma (i.e. compute Umnb)	
       allocate(Umnb(0:bl_hi-1, 0:N-1, 0:2*bl_hi-1), stat=fail)
       if(fail /= 0) then
          call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_core_synthesis_wav2flmv_dynamic')
       end if
       Umnb(0:bl_hi-1, 0:N-1, 0:2*bl_hi-1) = cmplx(0d0, 0d0)
       do gg = 0,N-1
          gamma_g = PI*gg/real(N,dp)
          do mmm = 0,N-1  	
             k_indicator = (1-(-1)**(N+mmm))/2 
             if(k_indicator > 0) then
                Umnb(0:bl_hi-1,mmm,0:2*bl_hi-1) = &
                     Umnb(0:bl_hi-1,mmm,0:2*bl_hi-1) &
                     + V(0:bl_hi-1,0:2*bl_hi-1,gg) * exp(-I*mmm*gamma_g)
             end if
          end do
       end do

       ! Step 2b: Integrate over beta	(i.e. compute Vmmm)
       V(0:bl_hi-1, 0:2*bl_hi-2, 0:N-1) = cmplx(0d0, 0d0)
       do bb = 0,2*bl_hi-1
          beta_b = PI*(2d0*bb+1d0)/(4d0*bl_hi)
          w = quad_weights(beta_b, bl_hi)
          do mmm = 0,N-1
             k_indicator = (1-(-1)**(N+mmm))/2 
             if(k_indicator > 0) then
                do mm = -(bl_hi-1),bl_hi-1  
                   V(0:bl_hi-1,mm+bl_hi-1,mmm) = V(0:bl_hi-1,mm+bl_hi-1,mmm) &
                        + Umnb(0:bl_hi-1, mmm, bb) * w * exp(-I*mm*beta_b)
                end do
             end if
          end do
       end do
       deallocate(Umnb)

       ! Step 3a: Combine Vmmm terms for m' and -m' outside of el sum
       do mmm = 0,N-1
          k_indicator = (1-(-1)**(N+mmm))/2	
          if(k_indicator > 0) then
             do m = 0,bl_hi-1
                do mm = 1,bl_hi-1    
                   V(m, mm+bl_hi-1, mmm) = V(m, mm+bl_hi-1, mmm) &
                        + (-1)**(m+mmm) *  V(m, -mm+bl_hi-1, mmm) 
                end do
             end do
          end if
       end do

       ! Step 3b: Compute Wigner coefficients
       allocate(wig(0:bl_hi-1, 0:bl_hi-1, 0:N-1), stat=fail)
       allocate(dl(-(bl_hi-1):bl_hi-1, -(bl_hi-1):bl_hi-1), stat=fail)
       if(fail /= 0) then
          call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, &
               's2dw_core_synthesis_wav2flmv_dynamic')
       end if
       wig(0:bl_hi-1, 0:bl_hi-1, 0:N-1) = cmplx(0d0, 0d0)
       do el = bl_lo+1,bl_hi-1
          call s2dw_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
          do m = 0,el
             do mmm = 0,min(N-1,el)			
                k_indicator = (1-(-1)**(N+mmm))/2 	
                if(k_indicator > 0) then
                   wig(el,m,mmm) = sum(V(m, bl_hi-1:el+bl_hi-1, mmm)  &
                        * dl(0:el,m) * dl(0:el,mmm)) * exp(-I*(mmm-m)*PION2)
                end if
             end do
          end do
       end do
       deallocate(dl, V)

       ! Compute harmonic component (for this jj) from Wigner coefficients
       do el = bl_lo+1,bl_hi-1
          do m = 0,el
             ! n=0 term
             flm(el,m) = flm(el,m) + &
                  (real(2d0*el+1,dp)/(8d0*pi**2)) &
                  * wig(el,m,0) * K_gamma(jj,el) * Slm(el,0)
             do mmm = 1,min(N-1,el)
                k_indicator = (1-(-1)**(N+mmm))/2 						
                if(k_indicator > 0) then
                   flm(el,m) = flm(el,m) + &
                        2d0 * (real(2d0*el+1,dp)/(8*pi**2)) &
                        * wig(el,m,mmm) * K_gamma(jj,el) * Slm(el,mmm) 
                end if
             end do
          end do
       end do
       deallocate(wig)

    end do

    ! Add in scaling part of signal
    do el = 0,bl_scoeff-1
       flm(el,0:bl_scoeff-1) = flm(el,0:bl_scoeff-1) &
            + Phi2(el) * scoeff(el,0:bl_scoeff-1)
    end do

    ! Free memory for array of structures for wavelet coefficients.
    deallocate(wavdyn)

  end subroutine s2dw_core_synthesis_wav2flm_dynamic


  !--------------------------------------------------------------------------
  ! Utility routines
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! s2dw_core_free_wavdyn
  !
  !! Free memory corresponding to dynamically allocated wavelet coefficients.
  !!
  !! Variables:
  !!  - wavdyn(0:J)%coeff: Wavelet coefficients for each scale (destroyed 
  !!    and freed on output) [input/output].
  !
  !! @author J. D. McEwen
  !! @version 0.1 February 2008
  !
  ! Revisions:
  !   February 2008 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  subroutine s2dw_core_free_wavdyn(wavdyn)

    type(s2dw_wav_abg), intent(inout), allocatable :: wavdyn(:)

    integer :: jj

    do jj=0,size(wavdyn)-1
       deallocate(wavdyn(jj)%coeff)
    end do

    deallocate(wavdyn)

  end subroutine s2dw_core_free_wavdyn


  !--------------------------------------------------------------------------
  ! s2dw_core_params_valid
  !
  !! Check size parameters valid.
  !!
  !! Notes:
  !!  - Only needs to be called once when initialising Kernels, then array 
  !!    sizes ensure other objects of correct size.
  !!
  !! Variables:
  !!  - J: Maximum analysis scale depth [input].
  !!  - B: Harmonic band limit [input].
  !   - alpha: Basis dilation factor [input].
  !!  - N: Azimuthal band limit [input].
  !!  - valid: Logical specifying whether passed parameter check [output].
  !
  !! @author J. D. McEwen
  !! @version 0.1 October 2007
  !
  ! Revisions:
  !   October 2007 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  function s2dw_core_params_valid(J, B, alpha, N) result(valid)

    integer, intent(in) :: J, B
    real(dp), intent(in) :: alpha
    integer, intent(in), optional :: N
    logical :: valid

    integer :: J_max

    valid = .true.

    J_max = s2dw_core_comp_Jmax(B, alpha)

    if(J>J_max .or. J<1) valid = .false.
    if(present(N)) then
       if(N>B) valid = .false.
    end if

  end function s2dw_core_params_valid


  !--------------------------------------------------------------------------
  ! s2dw_core_comp_Jmax
  !
  !! Compute maximum analysis depth allowed.
  !!
  !! Variables:
  !!  - B: Harmonic band limit [input].
  !   - alpha: Basis dilation factor [input].
  !!  - J_max: Maximum analysis scale depth allowed [output].
  !
  !! @author J. D. McEwen
  !! @version 0.1 November 2007
  !
  ! Revisions:
  !   November 2007 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  function s2dw_core_comp_Jmax(B, alpha) result(J_max)

    integer, intent(in) :: B
    real(dp), intent(in) :: alpha
    integer :: J_max

    !J_max = ceiling(log(real(B,dp))/log(alpha) - TOL_CEIL)
    if(.not. s2dw_core_assimilate_int(log(real(B,dp))/log(alpha), J_max)) then
       J_max = ceiling(log(real(B,dp))/log(alpha))
    end if

  end function s2dw_core_comp_Jmax


  !--------------------------------------------------------------------------
  ! s2dw_core_assimilate_int
  !
  !! If a real value is very close to an integer then set to integer, 
  !! otherwise leave unaltered.  Necessary to avoid numerical precision 
  !! problems when using ceiling and floor functions.
  !!
  !! Variables:
  !!  - x: Input value to (possibly) assimilate to integer [input].
  !!  - y: Value x is assimilated to (or 0 if not assimilated) [output].
  !!  - assimilated: Logical indicating where assimilation occured [output].
  !
  !! @author J. D. McEwen
  !! @version 0.1 November 2007
  !
  ! Revisions:
  !   November 2007 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  function s2dw_core_assimilate_int(x, y) result(assimilated)

    real(dp), intent(in) :: x
    integer, intent(out)  :: y
    logical:: assimilated

    if(abs( x - floor(x) ) < TOL_CEIL) then
       y = floor(x)
       assimilated = .true.
    elseif(abs( x - ceiling(x) ) < TOL_CEIL) then
       y = ceiling(x)
       assimilated = .true.
    else
       y = 0
       assimilated = .false.
    end if

  end function s2dw_core_assimilate_int


  !--------------------------------------------------------------------------
  ! quad_weights
  !
  !! Compute quadrature weights for exact quadrature with measure 
  !! dcos(theta).  Weights are derived by Driscoll and Healy.
  !!
  !! Variables:
  !!  - theta_t: Theta value to compute weight for [input].
  !!  - B: Harmonic band limit [input].
  !!  - w: Corresponding weight [output]
  !
  !! @author J. D. McEwen
  !! @version 0.1 October 2007
  !
  ! Revisions:
  !   October 2007 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  function quad_weights(theta_t, B) result(w)

    real(dp), intent(in) :: theta_t
    integer, intent(in) :: B
    real(dp) :: w	

    integer :: k

    w = 0d0
    do k = 0,B-1
       w = w + sin((2d0*k+1d0)*theta_t) / real(2d0*k+1d0,dp)
    end do
    w = (2d0/real(B,dp)) * sin(theta_t) * w

  end function quad_weights


  !--------------------------------------------------------------------------
  ! binomial
  !
  !! Compute Binomial coefficient nCr.
  !!
  !! Variables:
  !!  - n: Upper Binomial coefficient index [input].
  !!  - r: Lower Binomial coefficient index [input].
  !!  - b: Binomial coefficient value [output].
  !
  !! @author J. D. McEwen
  !! @version 0.1 October 2007
  !
  ! Revisions:
  !   October 2007 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  function binomial(n,r) result(b)

    integer, intent(in) :: n
    integer, intent(in) :: r
    integer :: b

    b = floor(0.5d0 + exp(logfact(n) - logfact(r) - logfact(n-r)))

  end function binomial


  !--------------------------------------------------------------------------
  ! logfact
  !
  !! Computes the natural logarithm of an (integer) factorial.
  !!
  !! Variables:
  !!  - n: Integer to compute factorial of [input].
  !!  - logfactn: Natural logarithm of factorial value computed [output].
  !
  !! @author J. D. McEwen
  !! @version 0.1 October 2007
  !
  ! Revisions:
  !   October 2007 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  function logfact(n) result(logfactn)
    integer, intent(in) :: n
    real(dp) :: logfactn

    real(dp) :: y, temp, sum, c(6), loggamma, x
    integer :: nn

    if (n < 0) then

       call s2dw_error(S2DW_ERROR_ARTH, 'logfact', &
            comment_add='Factorial argument negative')

    else

       ! The engine of this function actually calculates the gamma function,
       ! for which the real argument is x = n + 1.

       x = real(n, dp) + 1.0

       ! Table of fitting constants.

       c(1) = 76.18009172947146
       c(2) = - 86.50532032941677
       c(3) = 24.01409824083091
       c(4) = - 1.231739572450155
       c(5) = 0.1208650973866179e-2
       c(6) = - 0.5395239384953e-5

       ! Add up fit.

       temp = x + 5.5 - (x + 0.5) * log(x + 5.5);
       sum = 1.000000000190015
       y = x

       do nn = 1, 6
          y = y + 1.0;
          sum = sum + c(nn) / y;
       end do

       loggamma = - temp + log(2.5066282746310005 * sum / x);

    end if

    ! Finally make explicit the conversion back to log of the factorial.

    logfactn = loggamma

  end function logfact


  !--------------------------------------------------------------------------
  ! Numerical integration routines
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! trapzd
  !
  !! Computes nth stage of refinement of extended trapezoidal rule.
  !! Adapted from numerical recipes for computing scaling functions squared.
  !!
  !! Notes:
  !!   - Numerical recipies comment:
  !!     This routine computes the nth stage of refinement of an extended 
  !!     trapezoidal rule. func is input as the name of the function to be 
  !!     integrated between limits a and b, also input. When called with
  !!     n=1, the routine returns as s the crudest estimate of 
  !!     int_b^a f(x)dx. Subsequent calls with n=2,3,... (in that sequential
  !!     order) will improve the accuracy of s by adding 2n-2 additional
  !!     interior points. s should not be modified between sequential calls.
  !!   - Adapted for use in evaluating integrals required to compute
  !!     scaling functions squared (i.e. also takes parameters required 
  !!     to compute these integrals).
  !!
  !! Variables:
  !!   - func: "Pointer" to integrand function [input].
  !!   - B: Harmonic band limit parameter for integrand function [input].
  !!   - alpha:Basis dilation factor for integrand function [input].
  !!   - a_lim: Lower limit to evalue definite integral for [input].
  !!   - b_lim: Upper limit to evalue definite integral for [input].
  !!   - s: Value of integral computed to current order [input].
  !!   - n: Refinement order [output].
  !
  !! @author J. D. McEwen
  !! @version 0.1 June 2005
  !
  ! Revisions:
  !   June 2005 - Adapted from Numerical Recipes by Jason McEwen
  !--------------------------------------------------------------------------

  subroutine trapzd(func, B, alpha, a_lim, b_lim, s, n)

    integer, intent(in) :: B
    real(dp), intent(in) :: alpha
    real(dp), intent(IN) :: a_lim, b_lim
    real(dp), intent(INOUT) :: s
    integer, intent(IN) :: n
    interface
       function func(k, B, alpha) result(integrand)
         use s2dw_types_mod
         real(dp), intent(in) :: k, alpha
         integer, intent(in) :: B
         real(dp) :: integrand
       end function func
    end interface

    real(dp) :: del, fsum, x
    integer :: it, j

    if (n == 1) then
       s = 0.5d0 * (b_lim-a_lim)*(func(a_lim, B, alpha) + func(b_lim, B, alpha))
    else
       it = 2**(n-2)
       del = (b_lim-a_lim) / real(it, dp)  
       ! This is the spacing of the points to be added.
       x = a_lim + 0.5d0*del
       fsum = 0d0
       do j = 1,it
          fsum = fsum + func(x, B, alpha)
          x = x + del
       end do
       s = 0.5d0 * (s + (b_lim-a_lim)*fsum/real(it,dp)) 
       ! This replaces s by its refined value.
    end if

  end subroutine trapzd


  !--------------------------------------------------------------------------
  ! qtrap
  !
  !! Computes the integral of the function func from a to b using the 
  !! trapezoid method.  Adapted from numerical recipes for computing scaling 
  !! functions squared.
  !!
  !! Notes:
  !!   - Numerical recipies comment:
  !!     Returns the integral of the function func from a to b. 
  !!     The parameter EPS should be set to the desired fractional accuracy
  !!     and JMAX so that 2 to the power JMAX-1 is the maximum allowed
  !!     number of steps. Integration is performed by the trapezoidal rule.
  !!   - Adapted for use in evaluating integrals required to compute
  !!     scaling functions squared (i.e. also takes parameters required 
  !!     to compute these integrals).
  !!
  !! Variables:
  !!   - func: "Pointer" to integrand function [input].
  !!   - B: Harmonic band limit parameter for integrand function [input].
  !!   - alpha:Basis dilation factor for integrand function [input].
  !!   - a_lim: Lower limit to evalue definite integral for [input].
  !!   - b_lim: Upper limit to evalue definite integral for [input].
  !!   - integral: Value of the evaluated integral [output].
  !
  !! @author J. D. McEwen
  !! @version 0.1 June 2005
  !
  ! Revisions:
  !   June 2005 - Adapted from Numerical Recipes by Jason McEwen
  !--------------------------------------------------------------------------

  function qtrap(func, B, alpha, a_lim, b_lim) result(integral)

    integer, intent(in) :: B
    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: a_lim, b_lim
    interface
       function func(k, B, alpha) result(integrand)
         use s2dw_types_mod
         real(dp), intent(in) :: k, alpha
         integer, intent(in) :: B
         real(dp) :: integrand
       end function func
    end interface
    real(dp) :: integral

    integer, parameter :: JMAX=20
    real(dp) :: olds
    integer :: j

    olds = 0.0    !Initial value of olds is arbitrary.

    do j = 1,JMAX
       call trapzd(func, B, alpha, a_lim, b_lim, integral, j)
       if (j > 5) then     ! Avoid spurious early convergence.
          if (abs(integral-olds) < TOL_QUAD*abs(olds) .or. &
               (integral == 0.0 .and. olds == 0.0)) RETURN
       end if
       olds = integral
    end do

    ! If reach JMAX without finishing then call error.
    call s2dw_error(S2DW_ERROR_QUAD_STEP_EXCEED, 'qtrap')

  end function qtrap


  !--------------------------------------------------------------------------
  ! qsimp
  !
  !! Computes the integral of the function func from a to b using
  !! Simpson's rule.  Adapted from numerical recipes for computing scaling 
  !! functions squared.
  !!
  !! Notes:
  !!   - Numerical recipies comment:
  !!     Returns the integral of the function func from a to b.
  !!     The parameter EPS should be set to the desired fractional 
  !!     accuracy and JMAX so that 2 to the power JMAX-1 is the maximum
  !!     allowed number of steps. Integration is performed by Simpson's rule.
  !!   - Adapted for use in evaluating integrals required to compute
  !!     scaling functions squared (i.e. also takes parameters required 
  !!     to compute these integrals).
  !!
  !! Variables:
  !!   - func: "Pointer" to integrand function [input].
  !!   - B: Harmonic band limit parameter for integrand function [input].
  !!   - alpha:Basis dilation factor for integrand function [input].
  !!   - a_lim: Lower limit to evalue definite integral for [input].
  !!   - b_lim: Upper limit to evalue definite integral for [input].
  !!   - integral: Value of the evaluated integral [output].
  !
  !! @author J. D. McEwen
  !! @version 0.1 June 2005
  !
  ! Revisions:
  !   June 2005 - Adapted from Numerical Recipes by Jason McEwen
  !--------------------------------------------------------------------------

  function qsimp(func, B, alpha, a_lim, b_lim) result(integral)

    integer, intent(in) :: B
    real(dp), intent(in) :: alpha
    real(dp), intent(IN) :: a_lim, b_lim
    interface
       function func(k, B, alpha) result(integrand)
         use s2dw_types_mod
         real(dp), intent(in) :: k, alpha
         integer, intent(in) :: B
         real(dp) :: integrand
       end function func
    end interface
    real(dp) :: integral

    integer, parameter :: JMAX=40      
    integer :: j
    real(dp) :: os,ost,st

    ost=0.0
    os= 0.0
    do j=1,JMAX
       call trapzd(func, B, alpha, a_lim, b_lim, st, j)
       integral = (4d0*st-ost)/3d0     !Compare equation (4.2.4).
       if (j > 5) then                 !Avoid spurious early convergence.
          if (abs(integral-os) < TOL_QUAD*abs(os) .or. &
               (integral == 0.0 .and. os == 0.0)) RETURN
       end if
       os=integral
       ost=st
    end do

    ! If reach JMAX without finishing then call error.
    call s2dw_error(S2DW_ERROR_QUAD_STEP_EXCEED, 'qsimp')

  end function qsimp


end module s2dw_core_mod
