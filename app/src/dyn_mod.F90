!
! ISSP Math Library - A library for solving linear systems in materials science
! Copyright (C) 2016 Mitsuaki Kawamura
! 
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
! 
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
! 
! For more details, See `COPYING.LESSER' in the root directory of this library.
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
  !$ use omp_lib
  !
  USE shiftk_io, ONLY : input_restart_parameter, input_restart_vector, &
  &                     output_result, output_result_debug, output_restart_parameter, output_restart_vector
  USE shiftk_vals, ONLY : alpha, beta, calctype, ndim, nomega, maxloops, iter_old, solver, nl, myrank, &
  &                       outrestart, v12, v2, v14, v4, rhs, r_l, r_l_save, stdout, threshold, x_l, z, z_seed,&
  &                       v_n, Av_n
  !
  USE komega_cocg, ONLY : komega_COCG_init, komega_COCG_restart, komega_COCG_update, komega_COCG_getcoef, &
  &                       komega_COCG_getresidual, komega_COCG_getvec, komega_COCG_finalize
  USE komega_bicg, ONLY : komega_BiCG_init, komega_BiCG_restart, komega_BiCG_update, komega_BiCG_getcoef, &
  &                       komega_BiCG_getresidual, komega_BiCG_getvec, komega_BiCG_finalize
  USE komega_shifted_qmr_sym, ONLY: komega_shifted_qmr_sym_init, komega_shifted_qmr_sym_init_restart, &
          & komega_shifted_qmr_sym_update, komega_shifted_qmr_sym_output_restart, komega_shifted_qmr_sym_finalize, &
          & komega_shifted_qmr_sym_getresidual
  USE komega_shifted_qmr_sym_b, ONLY: komega_shifted_qmr_sym_b_init, komega_shifted_qmr_sym_b_init_restart, &
          & komega_shifted_qmr_sym_b_update, komega_shifted_qmr_sym_b_output_restart, komega_shifted_qmr_sym_b_finalize, &
          & komega_shifted_qmr_sym_b_getresidual
  USE lobpcg_mod, ONLY : zdotcMPI
#if defined(__MPI)
  USE mpi, ONLY : MPI_COMM_WORLD, MPI_DOUBLE_COMPLEX, MPI_WTIME
#endif
  !
  USE ham_prod_mod, ONLY : ham_prod
  !
  IMPLICIT NONE
  !
  INTEGER :: &
  & iter, & ! Counter for Iteration
  & status(3), &
  & fres = 30, &
  & iomega, ierr, ix
  !
  REAL(8) :: res2norm(nomega), xx
#if defined(__MPI)
  DOUBLE PRECISION :: t1, t2, t3
#endif
  double precision :: t01,t02,t03,t11,t12,t13
  logical :: lcollect=.FALSE.
  !
#if defined(NO_PROJ)
  INTEGER :: jter
  COMPLEX(8),allocatable :: test_r(:,:,:) 
#endif
  !
  IF (TRIM(solver) == 'lBiCG' .OR. TRIM(solver) == 'COCG') THEN
      nl = 1
  ELSE IF (TRIM(solver) == 'shifted_qmr_sym' .OR. TRIM(solver) == 'shifted_qmr_sym_b') THEN
      nl = ndim ! currently, shifted_qmr only works for nl = ndim
  END IF
#if defined(NO_PROJ)
  nl = ndim
#endif
  !
  ALLOCATE(v12(ndim), v2(ndim), r_l(nl), x_l(nl,nomega))
  IF(TRIM(solver) == "lBiCG") THEN
    ALLOCATE(v14(ndim), v4(ndim))
  ELSE IF (TRIM(solver) == "shifted_qmr_sym" .OR. TRIM(solver) == 'shifted_qmr_sym_b') THEN
    ALLOCATE(v_n(ndim), Av_n(ndim))
  END IF
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
        IF(TRIM(solver) == "lBiCG") THEN
#if defined(__MPI)
        CALL komega_BiCG_restart(ndim, nl, nomega, x_l, z, maxloops, threshold, status, &
           &                 iter_old, v2, v12, v4, v14, alpha, beta, z_seed, r_l_save, MPI_COMM_WORLD)
#else
           CALL komega_BiCG_restart(ndim, nl, nomega, x_l, z, maxloops, threshold, status, &
           &                 iter_old, v2, v12, v4, v14, alpha, beta, z_seed, r_l_save)
#endif
        ELSE IF (TRIM(solver) == "COCG") THEN
#if defined(__MPI)
           CALL komega_COCG_restart(ndim, nl, nomega, x_l, z, maxloops, threshold, status, &
           &                 iter_old, v2, v12,          alpha, beta, z_seed, r_l_save, MPI_COMM_WORLD)
#else
           CALL komega_COCG_restart(ndim, nl, nomega, x_l, z, maxloops, threshold, status, &
           &                 iter_old, v2, v12,          alpha, beta, z_seed, r_l_save)
#endif
        ELSE IF (TRIM(solver) == "shifted_qmr_sym") THEN
            CALL komega_shifted_qmr_sym_init_restart(ndim, nl, nomega, x_l, z, rhs(:), maxloops, threshold,&
                    & Av_n, v_n, 10, "output/qmr_restart.dat", status)
        ELSE IF (TRIM(solver) == "shifted_qmr_sym_b") THEN
            CALL komega_shifted_qmr_sym_b_init_restart(ndim, nl, nomega, x_l, z, rhs(:), maxloops, threshold,&
                    & Av_n, v_n, 10, "output/qmr_restart.dat", status)
        ELSE
            STOP
        END IF
     ELSE
         IF(TRIM(solver) == "lBiCG") THEN
#if defined(__MPI)
           CALL komega_BiCG_restart(ndim, nl, nomega, x_l, z, 0,        threshold, status, &
           &                 iter_old, v2, v12, v4, v14, alpha, beta, z_seed, r_l_save, MPI_COMM_WORLD)
#else
           CALL komega_BiCG_restart(ndim, nl, nomega, x_l, z, 0,        threshold, status, &
           &                 iter_old, v2, v12, v4, v14, alpha, beta, z_seed, r_l_save)
#endif
        ELSE IF(TRIM(solver) == "COCG") THEN
#if defined(__MPI)
           CALL komega_COCG_restart(ndim, nl, nomega, x_l, z, 0,        threshold, status, &
           &                 iter_old, v2, v12,          alpha, beta, z_seed, r_l_save, MPI_COMM_WORLD)
#else
           CALL komega_COCG_restart(ndim, nl, nomega, x_l, z, 0,        threshold, status, &
           &                 iter_old, v2, v12,          alpha, beta, z_seed, r_l_save)
#endif
        ELSE IF(TRIM(solver) == "shifted_qmr_sym") THEN
           CALL komega_shifted_qmr_sym_init_restart(ndim, nl, nomega, x_l, z, rhs(:), maxloops, threshold,&
                     & Av_n, v_n, 10, "output/qmr_restart.dat", status)
        ELSE IF(TRIM(solver) == "shifted_qmr_sym_b") THEN
            CALL komega_shifted_qmr_sym_b_init_restart(ndim, nl, nomega, x_l, z, rhs(:), maxloops, threshold,&
                    & Av_n, v_n, 10, "output/qmr_restart.dat", status)
        ELSE
           STOP
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
     IF(TRIM(solver) == 'lBiCG') v4(1:ndim) = rhs(1:ndim)
     !
     IF(outrestart .EQV. .TRUE.) THEN
         IF(TRIM(solver) == 'lBiCG') THEN
#if defined(__MPI)
           CALL komega_BiCG_init(ndim, nl, nomega, x_l, z, maxloops, threshold, MPI_COMM_WORLD)
#else
           CALL komega_BiCG_init(ndim, nl, nomega, x_l, z, maxloops, threshold)
#endif
        ELSE IF(TRIM(solver) == 'COCG') THEN
#if defined(__MPI)
           CALL komega_COCG_init(ndim, nl, nomega, x_l, z, maxloops, threshold, MPI_COMM_WORLD)
#else
           CALL komega_COCG_init(ndim, nl, nomega, x_l, z, maxloops, threshold)
#endif
        ELSE IF (TRIM(solver) == 'shifted_qmr_sym') THEN
             CALL komega_shifted_qmr_sym_init(ndim, nl, nomega, x_l, v_n, z, rhs(1:ndim), maxloops, threshold)
        ELSE IF (TRIM(solver) == 'shifted_qmr_sym_b') THEN
            CALL komega_shifted_qmr_sym_b_init(ndim, nl, nomega, x_l, v_n, z, rhs(1:ndim), maxloops, threshold)
        ELSE
             WRITE(*,*) 'The specified solver is not implemented'
             STOP
        END IF
     ELSE
        IF(TRIM(solver) == 'lBiCG') THEN
#if defined(__MPI)
           CALL komega_BiCG_init(ndim, nl, nomega, x_l, z, 0,        threshold, MPI_COMM_WORLD)
#else
           CALL komega_BiCG_init(ndim, nl, nomega, x_l, z, 0,        threshold)
#endif
        ELSE IF(TRIM(solver) == 'COCG') THEN
#if defined(__MPI)
           CALL komega_COCG_init(ndim, nl, nomega, x_l, z, 0,        threshold, MPI_COMM_WORLD)
#else
           CALL komega_COCG_init(ndim, nl, nomega, x_l, z, 0,        threshold)
#endif
        ELSE IF (TRIM(solver) == 'shifted_qmr_sym') THEN
            CALL komega_shifted_qmr_sym_init(ndim, nl, nomega, x_l, v_n, z, rhs(1:ndim), maxloops, threshold)
        ELSE IF (TRIM(solver) == 'shifted_qmr_sym_b') THEN
            CALL komega_shifted_qmr_sym_b_init(ndim, nl, nomega, x_l, v_n, z, rhs(1:ndim), maxloops, threshold)
        ELSE
            STOP
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
  IF(myrank == 0) OPEN(fres, file = "residual.dat")
  !
  WRITE(stdout,'(a)') "    iter status1 status2 status3      Residual       Proj. Res."
  !
#if defined(__MPI)
  OPEN(10, file='timer.dat', status='replace')
  t3 = 0.0
  !$ t03 = 0.0D0
  !$ t13 = 0.0D0
#else
  !$ OPEN(10, file='timer.dat', status='replace')
  !$ t03 = 0.0D0
  !$ t13 = 0.0D0
#endif
  DO iter = 1, maxloops
     !
     ! Projection of Residual vector into the space
     ! spaned by left vectors
     !
     !r_l(1) = zdotcMPI(ndim, rhs, v2)
     r_l(1) = DOT_PRODUCT(rhs, v2)
     !
#if defined(DEBUG)
     test_r(1:ndim,iter,1) = v2(1:ndim)
     IF(TRIM(solver) == 'lBiCG') THEN
        test_r(1:ndim,iter,2) = v4(1:ndim)
     ELSE IF (TRIM(solver) == 'COCG') THEN
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
#if defined(__MPI)
     !xx = 0.0
     !DO ix = 1, ndim
     !   xx = xx + v2(ix)*CONJG(v2(ix))
     !END DO
     !WRITE(*,*) myrank, xx
     !STOP
     !call MPI_BCAST(v2, ndim, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
     !call MPI_BCAST(v4, ndim, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
#endif
     !
#if defined(__MPI)
     t1 = MPI_WTIME()
#endif
     !$ t01 = omp_get_wtime()
     lcollect=.TRUE.
     !if(lcollect) call start_collection("region01")
     !if(lcollect) call stop_collection("region01")
     !$ t02 = omp_get_wtime()
     !$ t03 = t03 + (t02 - t01)
     !$ t13 = t13 + (t12 - t11)
     lcollect=.FALSE.
     IF(TRIM(solver) == 'lBiCG') THEN
         CALL ham_prod(v2, v12, t11, t12, lcollect)
         CALL ham_prod(v4, v14, t11, t12, lcollect)
     ELSE IF(TRIM(solver) == 'COCG') THEN
             CALL ham_prod(v2, v12, t11, t12, lcollect)
     ELSE IF(TRIM(solver) == 'shifted_qmr_sym' .OR. TRIM(solver) == 'shifted_qmr_sym_b') THEN
         CALL ham_prod(v_n, Av_n, t11, t12, lcollect)
     ELSE
         STOP
     END IF
     !
     ! Update result x with COCG
#if defined(__MPI)
     t2 = MPI_WTIME()
     t3 = t3 + (t2 - t1)
#endif
     !
     IF(TRIM(solver) == 'lBiCG') THEN
        CALL komega_BiCG_update(v12, v2, v14, v4, x_l, r_l, status)
     ELSE IF (TRIM(solver) == 'COCG') THEN
        CALL komega_COCG_update(v12, v2,          x_l, r_l, status)
     ELSE IF (TRIM(solver) == 'shifted_qmr_sym') THEN
        CALL komega_shifted_qmr_sym_update(Av_n, v_n, x_l, status)
     ELSE IF (TRIM(solver) == 'shifted_qmr_sym_b') THEN
         CALL komega_shifted_qmr_sym_b_update(Av_n, v_n, x_l, status)
     ELSE
        STOP
     END IF
     !
     IF(TRIM(solver) == 'lBiCG') THEN
        CALL komega_BiCG_getresidual(res2norm)
     ELSE IF (TRIM(solver) == 'COCG') THEN
        CALL komega_COCG_getresidual(res2norm)
     ELSE IF (TRIM(solver) == 'shifted_qmr_sym') THEN
        CALL komega_shifted_qmr_sym_getresidual(res2norm)
     ELSE IF (TRIM(solver) == 'shifted_qmr_sym_b') THEN
         CALL komega_shifted_qmr_sym_b_getresidual(res2norm)
     ELSE
         STOP
     END IF
     !
     !WRITE(*,*) myrank, 'v2', v2(1:5)
     !CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     !STOP
     !
     IF(myrank == 0) THEN
        DO iomega = 1, nomega
           WRITE(fres, '(2i8,5e25.15)') iter, iomega, z(iomega), x_l(1,iomega), res2norm(iomega)
        END DO
        WRITE(fres,*)
     END IF
     !
#if defined(__MPI)
     if( myrank.eq.0 .and. mod(iter,100).eq.0 ) then
        WRITE(10,*) t3
        t3 = 0.0
        !$ write(10,*) "timer:",iter,t03,t13
        !$ t03 = 0.0D0
        !$ t13 = 0.0D0
     END IF
#else
     if( mod(iter,100).eq.0 ) then
        !$ write(10,*) "timer:",iter,t03,t13
        !$ t03 = 0.0D0
        !$ t13 = 0.0D0
     end if
#endif
     !
     IF (TRIM(solver) == 'lBiCG' .OR. TRIM(solver) == 'COCG') THEN
         WRITE(stdout,'(i8,3i8,2e15.5)') iter, status, DBLE(v12(1)), ABS(r_l(1))
     ELSE IF (TRIM(solver) == 'shifted_qmr_sym' .OR. TRIM(solver) == 'shifted_qmr_sym_b') THEN
         WRITE(stdout,'(i8,3i8,2e15.5)') iter, status, DBLE(v_n(1)), ABS(res2norm(1))
     END IF
     IF(status(1) < 0) EXIT
     !
  END DO
#if defined(__MPI)
  CLOSE(10)
#else
  CLOSE(10)
#endif
  !
  CLOSE(fres)
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
        !& ABS(zdotcMPI(ndim, test_r(1:ndim,jter,2), test_r(1:ndim,iter,1)) )
        & ABS(DOT_PRODUCT(test_r(1:ndim,jter,2), test_r(1:ndim,iter,1)) )
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
     IF(TRIM(solver) == "lBiCG") THEN
        CALL komega_BiCG_getcoef(alpha, beta, z_seed, r_l_save)
        CALL komega_BiCG_getvec(v12,v14)
     ELSE IF(TRIM(solver) == "COCG") THEN
        CALL komega_COCG_getcoef(alpha, beta, z_seed, r_l_save)
        CALL komega_COCG_getvec(v12)
     END IF
     !
     CALL output_restart_parameter()
     CALL output_restart_vector()
     IF (TRIM(solver) == 'shifted_qmr_sym') THEN
         CALL komega_shifted_qmr_sym_output_restart(Av_n, v_n, 11, "output/qmr_restart.dat")
     ELSE IF (TRIM(solver) == 'shifted_qmr_sym_b') THEN
         CALL komega_shifted_qmr_sym_b_output_restart(Av_n, v_n, 11, "output/qmr_restart.dat")
     END IF
     !
     DEALLOCATE(alpha, beta, r_l_save)
     !     
  END IF
  !
  ! Deallocate all intrinsic vectors
  !
  IF(TRIM(solver) == 'lBiCG') THEN
     CALL komega_BiCG_finalize()
  ELSE IF (TRIM(solver) == 'COCG') THEN
     CALL komega_COCG_finalize()
  ELSE IF (TRIM(solver) == 'shifted_qmr_sym') THEN
      CALL komega_shifted_qmr_sym_finalize()
  ELSE IF (TRIM(solver) == 'shifted_qmr_sym_b') THEN
      CALL komega_shifted_qmr_sym_b_finalize()
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
  IF(TRIM(solver) == "lBiCG") DEALLOCATE(v14, v4)
  !
END SUBROUTINE dyn
!
END MODULE dyn_mod
