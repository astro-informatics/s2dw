!------------------------------------------------------------------------------
! s2dw_wstats
!
!! Compute wavelet coefficient statistics from a map (and then discard
!! wavelet coefficients).
!!
!! Usage: s2dw_wstats
!!   - [-help]: Display usage information.
!!   - [-inp filename_in]: Name of input file containing list of maps to 
!!     compute wavelet coefficients and then statistics of.
!!   - [-out filename_out]: Name of output statistic file.
!!   - [-B B]: Harmonic band limit (if not specified set to 2*nside).
!!   - [-N N]: Azimuthal band limit [default=3].
!!   - [-alpha alpha]: Basis dilation factor [default=2].
!!   - [-J J]: Maximum analysis scale [default=Jmax].
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   November 2012 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2dw_wstats

  use s2_sky_mod
  use s2dw_types_mod
  use s2dw_error_mod
  use s2dw_core_mod
  use s2dw_stat_mod
#ifdef MPI
  use mpi
#endif

  implicit none

  character(len=STRING_LEN) :: filename_in, filename_out
  character(len=STRING_LEN), allocatable :: description(:)
  character(len=STRING_LEN) :: line
  complex(dpc), allocatable :: flm(:,:)
  complex(spc), allocatable :: flm_sp(:,:)
  real(dp), allocatable :: K_gamma(:,:)
  real(dp), allocatable :: Phi2(:)
  complex(dpc), allocatable :: Slm(:,:)
  real(dp), allocatable :: admiss(:)
  type(s2dw_wav_abg), allocatable :: wavdyn(:)
  complex(dpc), allocatable :: scoeff(:,:)
  type(s2_sky) :: sky
  real(dp), allocatable :: mu(:,:), var(:,:), skew(:,:), kur(:,:)
  integer :: isim, nsim
  integer :: J
  integer :: J_max
  integer :: B
  integer :: N = 3
  integer :: bl_scoeff
  real(dp) :: alpha = 2d0
  logical :: admiss_pass
  logical :: use_Jmax
  integer :: fail = 0
  integer :: fileid, iostat
  character(len=STRING_LEN) :: error_string

#ifdef MPI
  integer :: status(1:MPI_STATUS_SIZE)
  integer :: rank, size, master, node, ierr, tag, start, end, nbox
  real(dp), allocatable :: boxed(:)
  integer :: ibox, jj
  integer :: MPI_REAL_DP
  call mpi_init(ierr)
  call mpi_type_create_f90_real(dp_p, dp_r, MPI_REAL_DP, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, size, ierr)
  master = 0
#endif

  ! Parse input parameters.
  call parse_options()

  ! Set J to maximum if not specified and check valid.
  J_max = ceiling(log(real(B,dp))/log(alpha) - TOL_CEIL)
  if(use_Jmax) J = J_max
  if(J > J_max) then
     J = J_max
     write(error_string,'(a,i4)') 'J too large, setting to maximum J = ', J_max
     call s2dw_error(S2DW_ERROR_ARG_WARNING, 's2dw_analysis', &
          comment_add=trim(error_string))
  end if

  ! Allocate memory.
  allocate(flm_sp(0:B-1,0:B-1), stat=fail)
  allocate(flm(0:B-1,0:B-1), stat=fail)
  allocate(K_gamma(0:J,0:B-1), stat=fail)
  allocate(Phi2(0:B-1), stat=fail)
  allocate(Slm(0:B-1,0:N-1), stat=fail)
  allocate(admiss(0:B-1), stat=fail)
  if(fail /= 0) then
     call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_analysis')
  end if

  ! Compute kernels, scaling function and directionality coefficients.
  call s2dw_core_init_kernels(K_gamma, Phi2, bl_scoeff, J, B, alpha)
  call s2dw_core_init_directionality(Slm, B, N)
  admiss_pass = s2dw_core_admiss(admiss, K_gamma, Phi2, B, J)
  if(.not. admiss_pass) then
     call s2dw_error(S2DW_ERROR_ADMISS_FAIL, 's2dw_analysis')
  end if

  ! Allocate memory for scaling coefficients (cannot be done earlier 
  ! since don't know bl_scoeff).
  allocate(scoeff(0:bl_scoeff-1, 0:bl_scoeff-1), stat=fail)
  if(fail /= 0) then
     call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_analysis')
  end if

  ! Open file containing list of maps.
  open(fileid, file=filename_in, form='formatted', status='old')

  ! Count number of simulations.
  nsim = 0
  do
     read(fileid,'(a)',iostat=iostat) line
     if (iostat < 0) exit
     nsim = nsim + 1
  end do
  rewind(fileid)

  ! Allocate memory for file descriptions.
  allocate(description(0:nsim-1), stat=fail)
  allocate(mu(0:nsim-1,0:J), stat=fail)
  allocate(var(0:nsim-1,0:J), stat=fail)
  allocate(skew(0:nsim-1,0:J), stat=fail)
  allocate(kur(0:nsim-1,0:J), stat=fail)
  if(fail /= 0) then
     call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_wstats')
  end if

  ! Read all filenames and close file.
  do isim = 0, nsim-1
     read(fileid,'(a)',iostat=iostat) description(isim)
  end do
  close(fileid)

  ! Read maps and compute moments. 
#ifdef MPI
  start = rank * nsim / size
  end = (rank+1) * nsim / size - 1
  if (rank == size-1) end = nsim - 1 
  do isim = start, end
#else
  do isim = 0, nsim-1
#endif

     ! Read sky.
     sky = s2_sky_init(description(isim), S2_SKY_FILE_TYPE_MAP)

     ! Compute spherical harmonic coefficients.
     call s2_sky_compute_alm(sky, B-1, B-1)
     call s2_sky_get_alm(sky, flm_sp(0:B-1,0:B-1))
     flm(0:B-1,0:B-1) = flm_sp(0:B-1,0:B-1)

     ! Compute wavelet transform.
     call s2dw_core_analysis_flm2wav_dynamic(wavdyn, scoeff, flm, &
          K_gamma, Slm, J, B, N, bl_scoeff, alpha)

     ! Compute wavelet coefficient statistics.
     call s2dw_stat_moments(wavdyn, J, B, N, alpha, &
          mu(isim,0:J), var(isim,0:J), skew(isim,0:J), kur(isim,0:J))

     ! Free temporary memory.
     call s2_sky_free(sky)
     call s2dw_core_free_wavdyn(wavdyn)

  end do

  ! Transfer and then save data.
#ifdef MPI
  if (rank /= master) then

     ! Box up data.
     nbox = 4 * (end-start+1) * (J+1)
     allocate(boxed(0:nbox-1), stat=fail)
     if(fail /= 0) then
        call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_wstats')
     end if
     ibox = 0
     do isim = start, end
        do jj = 0, J
           boxed(ibox) = mu(isim, jj)
           ibox = ibox + 1
           boxed(ibox) = var(isim, jj)
           ibox = ibox + 1
           boxed(ibox) = skew(isim, jj)
           ibox = ibox + 1
           boxed(ibox) = kur(isim, jj)
           ibox = ibox + 1
        end do
     end do

     ! Send wavelet statistics to master node.
     call mpi_send(boxed(0), nbox, MPI_REAL_DP, master, &
          rank, MPI_COMM_WORLD, ierr)

     ! Free boxed data.
     deallocate(boxed)

  else

     ! Read all descriptions on master.
     open(fileid, file=filename_in, form='formatted', status='old')
     do isim = 0, nsim-1
        read(fileid,'(a)',iostat=iostat) description(isim)
     end do
     close(fileid)

     ! Receive wavelet statistics on master node and then save.
     do node = 0, size-1
        if (node /= master) then

           ! Set up box to receive data.
           start = node * nsim / size
           end = (node+1) * nsim / size - 1
           if (node == size-1) end = nsim - 1 
           nbox = 4 * (end-start+1) * (J+1)
           allocate(boxed(0:nbox-1), stat=fail)
           if(fail /= 0) then
              call s2dw_error(S2DW_ERROR_MEM_ALLOC_FAIL, 's2dw_wstats')
           end if

           ! Receive data.
           call mpi_recv(boxed(0), nbox, MPI_REAL_DP, node, &
                node, MPI_COMM_WORLD, status, ierr)

           ! Unbox data.
           ibox = 0
           do isim = start, end
              do jj = 0, J
                 mu(isim, jj) = boxed(ibox)
                 ibox = ibox + 1
                 var(isim, jj) = boxed(ibox)
                 ibox = ibox + 1
                 skew(isim, jj) = boxed(ibox)
                 ibox = ibox + 1
                 kur(isim, jj) = boxed(ibox)
                 ibox = ibox + 1
              end do
           end do

           ! Free boxed data.
           deallocate(boxed)

        end if
     end do

#endif

     ! Save wavelet coefficient statistics.
     call s2dw_stat_moments_write(filename_out, nsim, description, J, &
          mu(0:nsim-1,0:J), var(0:nsim-1,0:J), &
          skew(0:nsim-1,0:J), kur(0:nsim-1,0:J))

#ifdef MPI
  end if
#endif

  ! Free memory.
  deallocate(description)
  deallocate(flm_sp, flm, K_gamma, Phi2, Slm, admiss, scoeff)
  deallocate(mu, var, skew, kur)

#ifdef MPI
  call mpi_finalize(ierr)
#endif


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
          write(*,'(a)') 'Usage: s2dw_wstat [-inp filename_in]'
          write(*,'(a)') '                  [-out filename_out]'
          write(*,'(a)') '                  [-B B]'  
          write(*,'(a)') '                  [-N N]'  
          write(*,'(a)') '                  [-alpha alpha]' 
          write(*,'(a)') '                  [-J J]'
#ifdef MPI
          call mpi_finalize(ierr)
#endif
          stop

       case ('-inp')
          filename_in = trim(arg)

       case ('-out')
          filename_out = trim(arg)

       case ('-B')
          read(arg,*) B

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


end program s2dw_wstats
