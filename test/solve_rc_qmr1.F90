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
MODULE solve_rc_qmr1_vals
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & rnd_seed, &
  & ndim,    & ! Size of Hilvert space
  & nz,      & ! Number of frequencies
  & nl,      & ! Number of Left vector
  & itermax, & ! Max. number of iteration
  & iter_old   ! Number of iteration of previous run
  !
  REAL(8),SAVE :: &
  & threshold ! Convergence Threshold
  !
  LOGICAL,SAVE :: &
  & restart !< It is restart run or not
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & z(:)         ! (nz): Frequency
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & ham(:,:), &
  & rhs(:), &
  & x_n(:,:), & ! (nl,nz) : Projected result
  & Av_n(:), &
  & v_n(:)
  !
END MODULE solve_rc_qmr1_vals
!
! Routines
!
MODULE solve_rc_qmr1_routines
  !
  IMPLICIT NONE
  !
CONTAINS
  !
SUBROUTINE input_size()
  !
  USE solve_rc_qmr1_vals, ONLY : ndim, nl, nz, itermax, threshold, rnd_seed, restart
  !
  IMPLICIT NONE
  !
  NAMELIST /input/ ndim, nz, itermax, threshold, rnd_seed, nl, restart
  !
  ndim = 5
  nl = 5
  nz = 1
  itermax = 5
  threshold = 1d-8
  rnd_seed = 1
  restart = .false.
  !
  READ(*,input,err=100)
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Standard Inputs  #####"
  WRITE(*,*)
  WRITE(*,*) "   Dimension of vvector : ", ndim
  WRITE(*,*) "  Number of frequencies : ", nz
  WRITE(*,*) "  Number of left vector : ", nl
  WRITE(*,*) "       Max. iterations : ", itermax
  WRITE(*,*) "              Threshold : ", threshold
  WRITE(*,*) "        Seed for Random : ", rnd_seed
  WRITE(*,*) "                Restart : ", restart
  !
  return
  !
100 write(*,*) "Stop in stdin. reading namelist file"
  !
  stop
  !
END SUBROUTINE input_size
!
! Generate Equations
!
SUBROUTINE generate_system()
  !
  USE solve_rc_qmr1_vals, ONLY : ndim, nz, ham, rhs, z, rnd_seed
  USE mathlib, ONLY : zgemm, zcopy
  !
  IMPLICIT NONE
  !
  INTEGER :: idim, iz, seedsize
  INTEGER,ALLOCATABLE :: seed(:)
  REAL(8) :: ham_r(ndim,ndim), rhs_r(ndim), ham_i(ndim,ndim), rhs_i(ndim), rnd(rnd_seed)
  COMPLEX(8) :: ham0(ndim,ndim)
  CHARACTER(100) :: cndim, form
  !
  CALL random_seed(size=seedsize) !Get size
  ALLOCATE(seed(seedsize)) !Allocate array
  seed(:) = rnd_seed
  CALL random_seed(put=seed)
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
  ham_r(1:ndim, 1:ndim) = 0d0
  ham_i(1:ndim, 1:ndim) = 0d0
  CALL RANDOM_NUMBER(ham_r(1, 1))
  CALL RANDOM_NUMBER(ham_i(1, 1))
  DO idim = 2, ndim
     CALL RANDOM_NUMBER(ham_r(idim, idim))
     CALL RANDOM_NUMBER(ham_r(idim, idim-1))
     CALL RANDOM_NUMBER(ham_i(idim, idim))
     CALL RANDOM_NUMBER(ham_i(idim, idim-1))
  END DO
  ham(1:ndim,1:ndim) = CMPLX(ham_r(1:ndim, 1:ndim), 0d0, KIND(0d0))
  !
  CALL zgemm("C", "N", ndim, ndim, ndim, CMPLX(1d0, 0d0, KIND(0d0)), ham, ndim, ham, ndim, CMPLX(0d0, 0d0, KIND(0d0)), ham0, ndim)
  CALL zcopy(ndim*ndim,ham0,1,ham,1)
  !
  CALL RANDOM_NUMBER(rhs_r(1:ndim))
  CALL RANDOM_NUMBER(rhs_i(1:ndim))
  rhs(1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
  rhs(1) = CMPLX(1d0, 1d0, KIND(0d0))
  rhs(1:ndim) = CMPLX(rhs_r(1:ndim), rhs_i(1:ndim), KIND(0d0))
  !
  WRITE(cndim,*) ndim * 2
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
! Check Result
!
SUBROUTINE output_result()
  !
  USE mathlib, ONLY : dgemv
  USE solve_rc_qmr1_vals, ONLY : ndim, nl, x_n, rhs, z, nz, ham
  !
  IMPLICIT NONE
  !
  INTEGER :: iz
  CHARACTER(100) :: cnl, form
  !
  WRITE(cnl,*) nl * 2
  WRITE(form,'(a,a,a)') "(", TRIM(ADJUSTL(cnl)), "e15.5)"
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Check Results  #####"
  WRITE(*,*)
  !
  WRITE(*,*) "  Resulting Vector"
  DO iz = 1, nz
     write(*,form) x_n(1:nl,iz)
  END DO
  !
  WRITE(*,*) "  Residual Vector"
  IF(nl /= ndim) THEN
     WRITE(*,*) "    Skip.  nl /= ndim."
     RETURN
  END IF
  !
END SUBROUTINE output_result
!
END MODULE solve_rc_qmr1_routines
!
!
!
PROGRAM solve_rc_qmr1
  !
  USE komega_shifted_qmr_sym, ONLY : komega_shifted_qmr_sym_init, komega_shifted_qmr_sym_update, &
  &                       komega_shifted_qmr_sym_getcoef, komega_shifted_qmr_sym_getvec, &
          &komega_shifted_qmr_sym_finalize, komega_shifted_qmr_sym_init_restart, komega_shifted_qmr_sym_output_restart
  USE solve_rc_qmr1_routines, ONLY : input_size, generate_system, &
  &                             output_result
  USE solve_rc_qmr1_vals, ONLY : ndim, nz, nl, itermax, iter_old, ham, restart, &
          &                         rhs, threshold, x_n, v_n, Av_n, z
  USE mathlib, ONLY : zgemv
  !
  IMPLICIT NONE
  !
  ! Variables for Restart
  !
  INTEGER :: &
  & iter, jter, & ! Counter for Iteration
  & status(2)
  !
  ! Input Size of vectors
  !
  CALL input_size()
  !
  ALLOCATE(x_n(nl,nz), v_n(nl), Av_n(nl), z(nz), ham(ndim,ndim), rhs(ndim))
  !
  CALL generate_system()
  !
  ! Check: Whether the restart file is exist.
  !
  WRITE(*,*)
  WRITE(*,*) "#####  CG Initialization  #####"
  WRITE(*,*)
  !
  IF(restart) THEN
    !
    ! When restarting, counter
    CALL komega_shifted_qmr_sym_init_restart(ndim, nl, nz, x_n, z, rhs(1:ndim), itermax, threshold, &
            &Av_n, v_n, 10, "restart.dat", status)
    ! These vectors were saved in COCG routine
    !
    IF(status(1) /= 0) GOTO 10
    !
  ELSE
     !
     ! Generate Right Hand Side Vector
     !
     CALL komega_shifted_qmr_sym_init(ndim, nl, nz, x_n, v_n, z, rhs(1:ndim), itermax, threshold)
     !
  END IF
  !
  ! COCG Loop
  !
  WRITE(*,*)
  WRITE(*,*) "#####  CG Iteration  #####"
  WRITE(*,*)
  !
  DO iter = 1, itermax
     !
     ! Matrix-vector product
    CALL zgemv('N', nl, nl, CMPLX(1d0, 0d0, KIND(0d0)), ham, nl, v_n, 1, CMPLX(0d0, 0d0, KIND(0d0)), Av_n, 1)
     !
     ! Update result x with COCG
     !
     CALL komega_shifted_qmr_sym_update(Av_n, v_n, x_n, status)
     !
     WRITE(*,'(a,i8,2i5,e15.5)') "DEBUG : ", iter, status, DBLE(Av_n(1))
!     IF(status(1) < 0) EXIT
     !
  END DO
  !
  IF(status(2) == 0) THEN
     WRITE(*,*) "  Converged in iteration ", ABS(status(1))
  ELSE IF(status(2) == 1) THEN
     WRITE(*,*) "  Not Converged in iteration ", ABS(status(1))
  END IF
  iter_old = abs(status(1))
  !
  ! Get these vectors for restart in the Next run
  !
  CALL komega_shifted_qmr_sym_output_restart(Av_n, v_n, 20, "restart.dat")
  !
10 CONTINUE
  !
  ! Deallocate all intrinsic vectors
  !
  CALL komega_shifted_qmr_sym_finalize()
  !
  ! Output to a file
  !
  CALL output_result()
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Done  #####"
  WRITE(*,*)
  !
END PROGRAM solve_rc_qmr1
