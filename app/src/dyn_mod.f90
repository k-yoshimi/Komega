!
!    Copyright 2016 Mitsuaki Kawamura
!
!    This file is part of ISSP Math Library.
!
!    ISSP Math Library is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ISSP Math Library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ISSP Math Library.  If not, see <http://www.gnu.org/licenses/>.
!
!
! Modules for the calculation of the spectrum
!
MODULE dyn_mod
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Main routine for the spectrum
!
SUBROUTINE dyn()
  !
  USE shiftk_io, ONLY : input_restart_parameter, input_restart_vector, &
  &                     output_result, output_result_debug, output_restart_parameter, output_restart_vector
  USE shiftk_vals, ONLY : alpha, beta, calctype, ndim, nomega, maxloops, iter_old, lBiCG, nl, &
  &                       outrestart, v12, v2, v14, v4, rhs, r_l, r_l_save, stdout, threshold, x_l, z, z_seed
  !
#if defined(MPI)
  USE pshifted_cocg, ONLY : pCOCG_init, pCOCG_restart, pCOCG_update, &
  &                         pCOCG_getcoef, pCOCG_getvec, pCOCG_finalize
  USE pshifted_bicg, ONLY : pBiCG_init, pBiCG_restart, pBiCG_update, &
  &                         pBiCG_getcoef, pBiCG_getvec, pBiCG_finalize
#else
  USE shifted_cocg, ONLY : COCG_init, COCG_restart, COCG_update, &
  &                        COCG_getcoef, COCG_getvec, COCG_finalize
  USE shifted_bicg, ONLY : BiCG_init, BiCG_restart, BiCG_update, &
  &                        BiCG_getcoef, BiCG_getvec, BiCG_finalize
#endif
  USE lobpcg_mod, ONLY : zdotcMPI
#if defined (MPI)
  USE mpi, ONLY : MPI_COMM_WORLD
#endif
  !
  USE ham_prod_mod, ONLY : ham_prod
  !
  IMPLICIT NONE
  !
  INTEGER :: &
  & iter, & ! Counter for Iteration
  & status(3)
  !
#if defined(NO_PROJ)
  INTEGER :: jter
  COMPLEX(8),allocatable :: test_r(:,:,:) 
#endif
  !
  nl = 1
#if defined(NO_PROJ)
  nl = ndim
#endif
  !lBiCG = .TRUE.
  !
  ALLOCATE(v12(ndim), v2(ndim), r_l(nl), x_l(nl,nomega))
  IF(lBiCG) ALLOCATE(v14(ndim), v4(ndim))
#if defined(DEBUG)
  ALLOCATE(test_r(ndim,maxloops,2))
#endif
  !
  ! Restart or frm scratch
  !
  IF(TRIM(calctype) == "recalc" .OR. TRIM(calctype) == "restart") THEN
     !
     ! For restarting with the previous result
     !
     CALL input_restart_parameter()
     IF(TRIM(calctype) == "restart") CALL input_restart_vector()
     maxloops = MAX(maxloops, iter_old)
     !
     WRITE(stdout,*)
     WRITE(stdout,*) "##########  CG Restart  ##########"
     WRITE(stdout,*)
     !
     IF(outrestart .EQV. .TRUE.) THEN
        IF(lBiCG) THEN
#if defined (MPI)
           CALL pBiCG_restart(ndim, nl, nomega, x_l, z, maxloops, threshold, MPI_COMM_WORLD, &
           &                 status, iter_old, v2, v12, v4, v14, alpha, beta, z_seed, r_l_save)
#else
           CALL BiCG_restart(ndim, nl, nomega, x_l, z, maxloops, threshold, &
           &                 status, iter_old, v2, v12, v4, v14, alpha, beta, z_seed, r_l_save)
#endif
        ELSE
#if defined (MPI)
           CALL pCOCG_restart(ndim, nl, nomega, x_l, z, maxloops, threshold, MPI_COMM_WORLD, &
           &                 status, iter_old, v2, v12,          alpha, beta, z_seed, r_l_save)
#else
           CALL COCG_restart(ndim, nl, nomega, x_l, z, maxloops, threshold, &
           &                 status, iter_old, v2, v12,          alpha, beta, z_seed, r_l_save)
#endif
        END IF
     ELSE
        IF(lBiCG) THEN
#if defined (MPI)
           CALL pBiCG_restart(ndim, nl, nomega, x_l, z, 0,        threshold, MPI_COMM_WORLD, &
           &                 status, iter_old, v2, v12, v4, v14, alpha, beta, z_seed, r_l_save)
#else
           CALL BiCG_restart(ndim, nl, nomega, x_l, z, 0,        threshold, &
           &                 status, iter_old, v2, v12, v4, v14, alpha, beta, z_seed, r_l_save)
#endif
        ELSE
#if defined (MPI)
           CALL pCOCG_restart(ndim, nl, nomega, x_l, z, 0,        threshold, MPI_COMM_WORLD, &
           &                 status, iter_old, v2, v12,          alpha, beta, z_seed, r_l_save)
#else
           CALL COCG_restart(ndim, nl, nomega, x_l, z, 0,        threshold, &
           &                 status, iter_old, v2, v12,          alpha, beta, z_seed, r_l_save)
#endif
        END IF
     END IF
     DEALLOCATE(alpha, beta, r_l_save)
     !
     IF(iter_old == maxloops .OR. TRIM(calctype) == "recalc") GOTO 10
     !
  ELSE IF(TRIM(calctype) == "normal") THEN
     !
     ! Compute from scratch
     !
     WRITE(stdout,*)
     WRITE(stdout,*) "##########  CG Initialization  ##########"
     WRITE(stdout,*)
     !
     v2(1:ndim) = rhs(1:ndim)
     IF(lBiCG) v4(1:ndim) = rhs(1:ndim)
     !
     IF(outrestart .EQV. .TRUE.) THEN
        IF(lBiCG) THEN
#if defined (MPI)
           CALL pBiCG_init(ndim, nl, nomega, x_l, z, maxloops, threshold, MPI_COMM_WORLD)
#else
           CALL BiCG_init(ndim, nl, nomega, x_l, z, maxloops, threshold)
#endif
        ELSE
#if defined (MPI)
           CALL pCOCG_init(ndim, nl, nomega, x_l, z, maxloops, threshold, MPI_COMM_WORLD)
#else
           CALL COCG_init(ndim, nl, nomega, x_l, z, maxloops, threshold)
#endif
        END IF
     ELSE
        IF(lBiCG) THEN
#if defined (MPI)
           CALL pBiCG_init(ndim, nl, nomega, x_l, z, 0,        threshold, MPI_COMM_WORLD)
#else
           CALL BiCG_init(ndim, nl, nomega, x_l, z, 0,        threshold)
#endif
        ELSE
#if defined (MPI)
           CALL pCOCG_init(ndim, nl, nomega, x_l, z, 0,        threshold, MPI_COMM_WORLD)
#else
           CALL COCG_init(ndim, nl, nomega, x_l, z, 0,        threshold)
#endif
        END IF
     END IF
     !
  ELSE
     !
     WRITE(stdout,*) "ERROR ! calctype = ", TRIM(calctype)
     STOP
     !
  END IF
  !
  ! COCG/BiCG Loop
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "#####  BiCG Iteration  #####"
  WRITE(stdout,*)
  !
  WRITE(stdout,'(a)') "    iter status1 status2 status3      Residual       Proj. Res."
  !
  DO iter = 1, maxloops
     !
     ! Projection of Residual vector into the space
     ! spaned by left vectors
     !
     r_l(1) = zdotcMPI(ndim, rhs, v2)
     !
#if defined(DEBUG)
     test_r(1:ndim,iter,1) = v2(1:ndim)
     IF(lBiCG) THEN
        test_r(1:ndim,iter,2) = v4(1:ndim)
     ELSE
        test_r(1:ndim,iter,2) = v2(1:ndim)
     END IF
#endif
     !
#if defined(NO_PROJ)
     r_l(1:ndim) = v2(1:ndim)
#endif
     !
     ! Matrix-vector product
     !
     CALL ham_prod(v2, v12)
     IF(lBiCG) CALL ham_prod(v4, v14)
     !
     ! Update result x with COCG
     !
     IF(lBiCG) THEN
#if defined (MPI)
        CALL pBiCG_update(v12, v2, v14, v4, x_l, r_l, status)
#else
        CALL BiCG_update(v12, v2, v14, v4, x_l, r_l, status)
#endif
     ELSE
#if defined (MPI)
        CALL pCOCG_update(v12, v2,          x_l, r_l, status)
#else
        CALL COCG_update(v12, v2,          x_l, r_l, status)
#endif
     END IF
     !
     WRITE(stdout,'(i8,3i8,2e15.5)') iter, status, DBLE(v12(1)), ABS(r_l(1))
     IF(status(1) < 0) EXIT
     !
  END DO
  !
  IF(status(2) == 0) THEN
     WRITE(stdout,*) "  Converged in iteration ", ABS(status(1))
  ELSE IF(status(2) == 1) THEN
     WRITE(stdout,*) "  Not Converged in iteration ", ABS(status(1))
  ELSE IF(status(2) == 2) THEN
     WRITE(stdout,*) "  Alpha becomes infinity", ABS(status(1))
  ELSE IF(status(2) == 3) THEN
     WRITE(stdout,*) "  Pi_seed becomes zero", ABS(status(1))
  ELSE IF(status(2) == 4) THEN
     WRITE(stdout,*) "  Residual & Shadow residual are orthogonal", ABS(status(1))
  END IF
  iter_old = ABS(status(1))
  !
10 CONTINUE
  !
#if defined(DEBUG)
  !
  ! Check othogonality of residual vectors
  !
  DO iter = 1, iter_old
     DO jter = 1, iter_old
        WRITE(stdout,'(e15.5)',advance="no") &
        & ABS(zdotcMPI(ndim, test_r(1:ndim,jter,2), test_r(1:ndim,iter,1)) )
     END DO
     WRITE(stdout,*)
  END DO
#endif
  !
  ! Get these vectors for restart in the Next run
  !
  IF(outrestart .EQV. .TRUE.) THEN
     !
     ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl, iter_old))
     !
     IF(lBiCG) THEN
#if defined (MPI)
        CALL pBiCG_getcoef(alpha, beta, z_seed, r_l_save)
        CALL pBiCG_getvec(v12,v14)
#else
        CALL BiCG_getcoef(alpha, beta, z_seed, r_l_save)
        CALL BiCG_getvec(v12,v14)
#endif
     ELSE
#if defined (MPI)
        CALL pCOCG_getcoef(alpha, beta, z_seed, r_l_save)
        CALL pCOCG_getvec(v12)
#else
        CALL COCG_getcoef(alpha, beta, z_seed, r_l_save)
        CALL COCG_getvec(v12)
#endif
     END IF
     !
     CALL output_restart_parameter()
     CALL output_restart_vector()
     !
     DEALLOCATE(alpha, beta, r_l_save)
     !     
  END IF
  !
  ! Deallocate all intrinsic vectors
  !
  IF(lBiCG) THEN
#if defined (MPI)
     CALL pBiCG_finalize()
#else
     CALL BiCG_finalize()
#endif
  ELSE
#if defined (MPI)
     CALL pCOCG_finalize()
#else
     CALL COCG_finalize()
#endif
  END IF
  !
  ! Output to a file
  !
#if defined(NO_PROJ)
  CALL output_result_debug()
#else
  CALL output_result()
#endif
  !
  DEALLOCATE(v12, v2, r_l, x_l, z)
  IF(lBiCG) DEALLOCATE(v14, v4)
  !
END SUBROUTINE dyn
!
END MODULE dyn_mod
