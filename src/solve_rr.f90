MODULE solve_rr_vals
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & rnd_seed, &
  & ndim,    & ! Size of Hilvert space
  & nz,      & ! Number of frequencies
  & nl,      & ! Number of Left vector
  & itermax, & ! Max. number of iteraction
  & iter_old   ! Number of iteraction of previous run
  !
  REAL(8),SAVE :: &
  & z_seed, & ! Seed frequency
  & threshold ! Convergence Threshold
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & ham(:,:), &
  & rhs(:), &
  & v12(:), v2(:), & ! (ndim): Working vector
  & r_l(:), & ! (nl) : Projeccted residual vector 
  & x(:,:),    & ! (nl,nz) : Projected result 
  & z(:)         ! (nz): Frequency
  !
  ! Variables for Restart
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & alpha(:), beta(:), & ! (iter_old) 
  & r_l_save(:,:) ! (nl,iter_old) Projected residual vectors
  !
END MODULE solve_rr_vals
!
! Routines
!
MODULE solve_rr_routines
  !
  IMPLICIT NONE
  !
CONTAINS
  !
SUBROUTINE input_size()
  !
  USE solve_rr_vals, ONLY : ndim, nl, nz, itermax, threshold, rnd_seed
  !
  IMPLICIT NONE
  !
  NAMELIST /input/ ndim, nz, itermax, threshold, rnd_seed
  !
  ndim = 5
  nz = 1
  itermax = 0
  threshold = 1d-8
  rnd_seed = 1
  !
  READ(*,input,err=100)
  nl = ndim
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Standard Inputs  #####"
  WRITE(*,*)
  WRITE(*,*) "   Dimension of vvector : ", ndim
  WRITE(*,*) "  Number of frequencies : ", nz
  WRITE(*,*) "  Number of left vector : ", nl
  WRITE(*,*) "       Max. iteractions : ", itermax
  WRITE(*,*) "              Threshold : ", threshold
  WRITE(*,*) "        Seed for Random : ", threshold
  !
  return
  !
100 write(*,*) "Stop in stdin. reading namelist file"
  !
  stop
  !
END SUBROUTINE input_size
!
! Input restart variables from file
!
SUBROUTINE input_restart()
  !
  USE solve_rr_vals, ONLY : iter_old, v2, v12, alpha, beta, z_seed, r_l_save, nl, ndim
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ierr
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Check Restart File  #####"
  WRITE(*,*)
  !
  open(fi, file = 'restart.dat',status="old", action = 'read',iostat = ierr)
  !
  IF(ierr == 0) THEN
     !
     WRITE(*,*) "  Restart file is found."
     !
     READ(fi,*) iter_old
     WRITE(*,*) "  iter_old : ", iter_old
     !
     ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl,iter_old))
     READ(fi,*) z_seed
     READ(fi,*) alpha(1:iter_old)
     READ(fi,*) beta(1:iter_old)
     READ(fi,*) r_l_save(1:nl, 1:iter_old)
     READ(fi,*) v2(1:ndim)
     READ(fi,*) v12(1:ndim)
     !
     close(fi)
     !
  ELSE
     !
     WRITE(*,*) "  Restart file is NOT found."
     iter_old = 0
     !
  END IF
  !
END SUBROUTINE input_restart
!
! Generate Equations
!
SUBROUTINE generate_system()
  !
  USE solve_rr_vals, ONLY : ndim, nz, ham, rhs, z, rnd_seed
  !
  IMPLICIT NONE
  !
  INTEGER :: idim, iz
  REAL(8) :: ham0(ndim,ndim), rnd(rnd_seed)
  CHARACTER(100) :: cndim, form
  !
  CALL RANDOM_NUMBER(rnd(1:rnd_seed))
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Generate Linear System  #####"
  WRITE(*,*)
  !
  READ(*,*) z(1:nz)
  WRITE(*,*) "  Frequency :"
  DO iz = 1, nz
     WRITE(*,*) iz, z(iz)
  END DO
  !
  !CALL RANDOM_NUMBER(ham(1:ndim, 1:ndim))
  ham(1:ndim, 1:ndim)=0d0
  CALL RANDOM_NUMBER(ham(1, 1))
  DO idim = 2, ndim
     CALL RANDOM_NUMBER(ham(idim, idim))
     CALL RANDOM_NUMBER(ham(idim, idim-1))
  END DO
  CALL RANDOM_NUMBER(rhs(1:ndim))
  !
  CALL dgemm("T", "N", ndim, ndim, ndim, 1d0, ham, ndim, ham, ndim, 0d0, ham0, ndim)
  !
  CALL dcopy(ndim*ndim,ham0,1,ham,1)
  !
  WRITE(cndim,*) ndim
  WRITE(form,'(a,a,a)') "(", TRIM(ADJUSTL(cndim)), "e15.5)"
  !
  WRITE(*,*) 
  WRITE(*,*) "  Right Hand Side Vector :"
  WRITE(*,form) rhs(1:ndim)
  !
  WRITE(*,*) 
  WRITE(*,*) "  Matrix :"
  DO idim = 1, ndim
     WRITE(*,form) ham(1:ndim,idim)
  END DO
  !
END SUBROUTINE generate_system
!
! Output variables for restart
!
SUBROUTINE output_restart()
  !
  USE solve_rr_vals, ONLY : iter_old, v2, v12, alpha, beta, z_seed, r_l_save, nl, ndim
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Output Restart File  #####"
  WRITE(*,*)
  !
  OPEN(fo, file = 'restart.dat', action = 'write')
  !
  WRITE(fo,*) iter_old
  WRITE(fo,*) z_seed
  WRITE(fo,*) alpha(1:iter_old)
  WRITE(fo,*) beta(1:iter_old)
  WRITE(fo,*) r_l_save(1:nl, 1:iter_old)
  WRITE(fo,*) v2(1:ndim)
  WRITE(fo,*) v12(1:ndim)
  !
  close(fo)
  !
  WRITE(*,*) "  Restart File is written."
  !
END SUBROUTINE output_restart
!
! Check Result
!
SUBROUTINE output_result()
  !
  USE mathlib, ONLY : dgemv
  USE solve_rr_vals, ONLY : v2, ndim, x, rhs, z, nz, ham
  !
  IMPLICIT NONE
  !
  INTEGER :: iz
  CHARACTER(100) :: cndim, form
  !
  WRITE(cndim,*) ndim
  WRITE(form,'(a,a,a)') "(", TRIM(ADJUSTL(cndim)), "e15.5)"
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Check Results  #####"
  WRITE(*,*)
  !
  DO iz = 1, nz
     !
     v2(1:ndim) = z(iz) * x(1:ndim,iz) - rhs(1:ndim)
     CALL dgemv("N", ndim, ndim, -1d0, Ham, ndim, x(1:ndim,iz), 1, 1d0, v2, 1)
     !
     write(*,form) v2(1:ndim)
     !
  END DO
  !
END SUBROUTINE output_result
!
END MODULE solve_rr_routines
!
!
!
PROGRAM solve_rr
  !
  USE shifted_krylov, ONLY : CG_R_init, CG_R_restart, CG_R_update, &
  &                          CG_R_getcoef, CG_R_getvec, CG_R_finalize
  USE solve_rr_routines, ONLY : input_size, input_restart, generate_system, &
  &                              output_restart, output_result
  USE solve_rr_vals, ONLY : alpha, beta, ndim, nz, nl, itermax, iter_old, ham, &
  &                         rhs, v12, v2, r_l, r_l_save, threshold, x, z, z_seed
  USE mathlib, ONLY : dgemv, ddot
  !
  IMPLICIT NONE
  !
  ! Variables for Restart
  !
  INTEGER :: &
  & itermin, & ! First iteration in this run
  & iter,    & ! Counter for Iteration
  & status(3)
  !
  ! Input Size of vectors
  !
  CALL input_size()
  !
  ALLOCATE(v12(ndim), v2(ndim), r_l(nl), x(nl,nz), z(nz), ham(ndim,ndim), rhs(ndim))
  !
  CALL generate_system()
  !
  ! Check: Whether the restart file is exist.
  !
  CALL input_restart()
  !
  WRITE(*,*)
  WRITE(*,*) "#####  CG Initialization  #####"
  WRITE(*,*)
  !
  IF(iter_old > 0) THEN
    !
    ! When restarting, counter
    !
    itermin = iter_old + 1
    CALL CG_R_restart(ndim, nl, nz, x, z, max(0,itermax), threshold, status, &
    &                 iter_old, v2, v12, alpha, beta, z_seed, r_l_save)
    !
    ! These vectors were saved in CG_R routine
    !
    DEALLOCATE(alpha, beta, r_l_save)
    !
    IF(status(1) /= 0) GOTO 10
    !
  ELSE
     !
     itermin = 1
     !
     ! Generate Right Hand Side Vector
     !
     v2(1:ndim) = rhs(1:ndim)
     !
     CALL CG_R_init(ndim, nl, nz, x, z, max(0,itermax), threshold)
     !
  END IF
  !
  ! CG_R Loop
  !
  WRITE(*,*)
  WRITE(*,*) "#####  CG Iteration  #####"
  WRITE(*,*)
  !
  DO iter = 1, abs(itermax)
     !
     ! Projection of Residual vector into the space
     ! spaned by left vectors
     !
     r_l(1:nl) = v2(1:ndim)
     !
     ! Matrix-vector product
     !
     CALL dgemv("N", ndim, ndim, 1d0, Ham, ndim, v2, 1, 0d0, v12, 1)
     !
     ! Update result x with CG_R
     !
     CALL CG_R_update(v12, v2, x, r_l, status)
     !
     WRITE(*,'(a,4i,e15.5)') "DEBUG : ", iter, status, v12(1)
     IF(status(1) /= 0) EXIT
     !
  END DO
  !
  IF(status(1) < 0) THEN
     WRITE(*,*) "  Not Converged in iteration ", -status(1)
  ELSE
     WRITE(*,*) "  Converged in iteration ", status(1)
  END IF
  iter_old = abs(status(1))
  !
  ! Get these vectors for restart in the Next run
  !
  IF(itermax > 0) THEN
     !
     ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl,iter_old))
     !
     CALL CG_R_getcoef(alpha, beta, z_seed, r_l_save)
     CALL CG_R_getvec(v12)
     !
     CALL output_restart()
     !
     DEALLOCATE(alpha, beta, r_l_save)
     !     
  END IF
  !
10 CONTINUE
  !
  ! Deallocate all intrinsic vectors
  !
  CALL CG_R_finalize()
  !
  ! Output to a file
  !
  CALL output_result()
  !
  DEALLOCATE(v12, v2, r_l, x, z)
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Done  #####"
  WRITE(*,*)
  !
END PROGRAM solve_rr
