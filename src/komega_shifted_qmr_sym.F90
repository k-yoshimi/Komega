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
! Routines for COCG
!
MODULE komega_shifted_qmr_sym
  !
  PRIVATE
  !
  PUBLIC komega_shifted_qmr_sym_init, komega_shifted_qmr_sym_restart, komega_shifted_qmr_sym_update, &
  & komega_shifted_qmr_sym_getcoef, komega_shifted_qmr_sym_getvec, komega_shifted_qmr_sym_finalize, &
  & komega_shifted_qmr_sym_getresidual, komega_shifted_qmr_sym_output_restart, &
  & komega_shifted_qmr_sym_init_restart, komega_shifted_qmr_getresidual
  !
CONTAINS
!
! Allocate & initialize variables
!
#if defined(__MPI)
SUBROUTINE komega_shifted_qmr_sym_init(ndim0, nl0, nz0, x, v_n, z0, b0, itermax0, threshold0, comm0) BIND(C)
#else
SUBROUTINE komega_shifted_qmr_sym_init(ndim0, nl0, nz0, x, v_n, z0, b0, itermax0, threshold0) BIND(C)
#endif
  !
  USE ISO_C_BINDING
  USE komega_parameter, ONLY : iter, itermax, ndim, nl, nz, &
  &                            threshold, iz_seed, lz_conv, lmpi, comm
  USE komega_vals_shifted_qmr_sym, ONLY : alpha_n, beta_n, beta_nmin1, c_n, c_nmin1, c_nmin2, g_n, &
         & g_npls1, p_n, p_nmin1, p_nmin2, r_nmin1, s_n, s_nmin1, s_nmin2, t_n_n, t_nmin1_n, t_nmin2_n, &
         & t_npls1_n, w_n, w_nmin1, x_nmin1, b, res_fast, t_nmin2_nmin2, t_nmin1_nmin1, &
         & v_nmin1, v_npls1,  z
  USE komega_math, ONLY : dcopy, zdotcMPI
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ndim0, nl0, nz0, itermax0
  REAL(8),INTENT(IN) :: threshold0
  COMPLEX(8),INTENT(IN) :: z0(nz0)
  COMPLEX(8),INTENT(IN) :: b0(nl0)
  COMPLEX(8),INTENT(OUT) :: x(nl0,nz0), v_n(nl0)
#if defined(__MPI)
  INTEGER(C_INT),INTENT(IN) :: comm0
#endif
  ndim = ndim0
  nl = nl0
  nz = nz0
  itermax = itermax0
  threshold = threshold0
  !
  comm = 0
  lmpi = .FALSE.
#if defined(__MPI)
  comm = comm0
  lmpi = .TRUE.
#endif
  !
  ALLOCATE(z(nz), c_n(nz), c_nmin1(nz), c_nmin2(nz), g_n(nz), g_npls1(nz), p_n(nl, nz), p_nmin1(nl, nz),&
          & p_nmin2(nl, nz), r_nmin1(nz), s_n(nz), &
          & s_nmin1(nz), s_nmin2(nz), t_n_n(nz), t_nmin1_n(nz), t_nmin2_n(nz), t_npls1_n(nz), &
          & w_n(nl, nz), w_nmin1(nl, nz), x_nmin1(nl, nz), lz_conv(nz), b(nl), res_fast(nz, itermax0), &
          & t_nmin2_nmin2(nz), t_nmin1_nmin1(nz), v_nmin1(nl), v_npls1(nl))
  CALL zcopy(nz,z0,1,z,1)

  c_n(1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  c_nmin1(1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  g_n(1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  g_npls1(1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  p_n(1:nl, 1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  r_nmin1(1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  s_n(1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  s_nmin1(1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  t_n_n(1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  t_nmin1_n(1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  t_nmin2_n(1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  t_npls1_n(1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  w_n(1:nl, 1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  w_nmin1(1:nl, 1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  x(1:nl, 1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  x_nmin1(1:nl, 1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  lz_conv(1:nz) = .FALSE.

  alpha_n = CMPLX(1d0, 0d0, KIND(0d0))
  beta_n = CMPLX(1d0, 0d0, KIND(0d0))
  beta_nmin1 = CMPLX(1d0, 0d0, KIND(0d0))

  iter = 0
  lz_conv(1:nz) = .FALSE.
  CALL zcopy(nl,b0, 1, b, 1)
  v_n(:) = b(:) / SQRT(zdotcMPI(nl, b, b))
  g_n(:) = SQRT(zdotcMPI(nl, b, b))
  !
END SUBROUTINE komega_shifted_qmr_sym_init
!
! Update x, p, r
!
SUBROUTINE komega_shifted_qmr_sym_update(Av_n, v_n, x_n, status)
  USE komega_parameter, ONLY : iter, itermax, ndim, nl, nz, &
  &                            threshold, almost0, lz_conv
  USE komega_vals_shifted_qmr_sym, ONLY : alpha_n, beta_n, beta_nmin1, c_n, c_nmin1, c_nmin2, g_n, &
          & g_npls1, p_n, p_nmin1, p_nmin2, r_nmin1, s_n, s_nmin1, s_nmin2, t_n_n, t_nmin1_n, t_nmin2_n, &
          & t_npls1_n, w_n, w_nmin1, x_nmin1, b, res_fast, t_nmin2_nmin2, t_nmin1_nmin1, v_nmin1, &
          & v_npls1, z
  USE komega_math, ONLY : zcopy, zdotcMPI
  !
  IMPLICIT NONE
  !
  COMPLEX(8), INTENT(IN) :: Av_n(nl)
  COMPLEX(8), INTENT(INOUT) :: v_n(nl), x_n(nl, nz)
  INTEGER,INTENT(INOUT) :: status(2)
  ! Working variables
  INTEGER :: iz
  COMPLEX(8), DIMENSION(nl) :: v_tld_npls1
  LOGICAL flag_converge
  COMPLEX(8), DIMENSION(2) :: ttmp
  !
  iter = iter + 1
  status(1:2) = 0

  alpha_n = zdotcMPI(nl, v_n, -Av_n)
  IF (iter /= 1) THEN
          v_tld_npls1(:) = -Av_n(:) - alpha_n * v_n(:) - beta_nmin1 * v_nmin1(:)
  ELSE
    v_tld_npls1(:) = -Av_n(:) - alpha_n * v_n(:)
  END IF
  beta_n = SQRT(zdotcMPI(nl, v_tld_npls1, v_tld_npls1))
  v_npls1(:) = v_tld_npls1(:) / beta_n

  IF (iter /= 1) THEN
    t_nmin1_n(:) = beta_nmin1
  END IF

  t_n_n(:) = alpha_n + z(:)
  t_npls1_n(:) = beta_n

  DO iz = 1, nz
    ! compute residuals
    IF (iter == 1) THEN
      w_n(:, iz) = v_n(:)
    ELSE
      w_n(:, iz) = -s_nmin1(iz) * w_nmin1(:, iz) + c_nmin1(iz) * v_n(:)
    END IF
    r_nmin1(iz) = ABS(g_n(iz)) * SQRT(zdotcMPI(nl, w_n(:, iz), w_n(:, iz)))
    res_fast(iz, iter) = r_nmin1(iz)

    IF (iter /= 1 .AND. r_nmin1(iz) / SQRT(DBLE(zdotcMPI(nl, b, b))) < threshold) THEN
      lz_conv(iz) = .TRUE.
    ELSE
      IF (iter == 2) THEN
        ttmp(1) = c_nmin1(iz) * t_nmin1_n(iz) + s_nmin1(iz) * t_n_n(iz)
        ttmp(2) = -CONJG(s_nmin1(iz)) * t_nmin1_n(iz) + c_nmin1(iz) * t_n_n(iz)
        t_nmin1_n(iz) = ttmp(1)
        t_n_n(iz) = ttmp(2)
      ELSE IF (iter >= 3) THEN
        ttmp(1) = s_nmin2(iz) * t_nmin1_n(iz)
        ttmp(2) = c_nmin2(iz) * t_nmin1_n(iz)
        t_nmin2_n(iz) = ttmp(1)
        t_nmin1_n(iz) = ttmp(2)
        ttmp(1) = c_nmin1(iz) * t_nmin1_n(iz) + s_nmin1(iz) * t_n_n(iz)
        ttmp(2) = -CONJG(s_nmin1(iz)) * t_nmin1_n(iz) + c_nmin1(iz) * t_n_n(iz)
        t_nmin1_n(iz) = ttmp(1)
        t_n_n(iz) = ttmp(2)
      END IF
      c_n(iz) = ABS(t_n_n(iz)) / SQRT(ABS(t_n_n(iz))**2 + ABS(t_npls1_n(iz))**2)
      s_n(iz) = CONJG(t_npls1_n(iz)/t_n_n(iz) * c_n(iz))
      t_n_n(iz) = c_n(iz) * t_n_n(iz) + s_n(iz) * t_npls1_n(iz)
      t_npls1_n(iz) = 0

      ttmp(1) = -CONJG(s_n(iz)) * g_n(iz)
      ttmp(2) = c_n(iz) * g_n(iz)
      g_npls1(iz) = ttmp(1)
      g_n(iz) = ttmp(2)
      IF (iter == 1) THEN
        p_n(:, iz) = v_n(:)
      ELSE IF (iter == 2) THEN
        p_n(:,iz) = v_n(:) - (t_nmin1_n(iz)/t_nmin1_nmin1(iz)) * p_nmin1(:, iz)
      ELSE
        p_n(:,iz) = v_n(:) - (t_nmin2_n(iz) / (t_nmin2_nmin2(iz))) * p_nmin2(:,iz) &
                & - (t_nmin1_n(iz)/t_nmin1_nmin1(iz)) * p_nmin1(:,iz)
      END IF
      IF (iter == 1) THEN
        x_n(:,iz) = (g_n(iz)/t_n_n(iz)) * p_n(:,iz)
      ELSE
        x_n(:,iz) = x_nmin1(:,iz) + (g_n(iz)/t_n_n(iz)) * p_n(:,iz)
      END IF
    END IF
  END DO

  flag_converge = .TRUE.
  DO iz = 1,nz
    IF (lz_conv(iz) .EQV. .FALSE.) THEN
      flag_converge = .FALSE.
    END IF
  END DO

  IF (flag_converge .EQV. .TRUE.) THEN
    status(1) = -iter
    status(2) = 0
  ELSE IF (iter == itermax) THEN
    status(1) = iter
    status(2) = 1
  ELSE
    status(1) = iter
    status(2) = 0
  END IF

  beta_nmin1 = beta_n
  c_nmin2(:) = c_nmin1(:)
  c_nmin1(:) = c_n(:)
  g_n(:) = g_npls1(:)
  IF (iter >= 2) THEN
    p_nmin2(:,:) = p_nmin1(:,:)
  END IF
  p_nmin1(:,:) = p_n(:,:)
  s_nmin2(:) = s_nmin1(:)
  s_nmin1(:) = s_n(:)
  IF (iter >= 2) THEN
    t_nmin2_nmin2(:) = t_nmin1_nmin1(:)
  END IF
  t_nmin1_nmin1(:) = t_n_n(:)
  v_nmin1(:) = v_n(:)
  v_n(:) = v_npls1(:)
  w_nmin1(:,:) = w_n(:,:)
  x_nmin1(:,:) = x_n(:,:)
  !
END SUBROUTINE komega_shifted_qmr_sym_update
!
! Deallocate private arrays
!
SUBROUTINE komega_shifted_qmr_sym_finalize()
  !
  USE komega_parameter, ONLY : itermax
  USE komega_vals_shifted_qmr_sym, ONLY: alpha_n, beta_n, beta_nmin1, c_n, c_nmin1, c_nmin2, g_n, &
          & g_npls1, p_n, p_nmin1, p_nmin2, r_nmin1, s_n, s_nmin1, s_nmin2, t_n_n, t_nmin1_n, t_nmin2_n, &
          & t_npls1_n, w_n, w_nmin1, x_nmin1, b, res_fast, t_nmin2_nmin2, t_nmin1_nmin1, &
          & v_nmin1, v_npls1,  z
  USE komega_parameter, ONLY : lz_conv
  !
  IMPLICIT NONE
  !
  IF(itermax > 0) THEN
    DEALLOCATE(z, c_n, c_nmin1, c_nmin2, g_n, g_npls1, p_n, p_nmin1,&
            & p_nmin2, r_nmin1, s_n, &
            & s_nmin1, s_nmin2, t_n_n, t_nmin1_n, t_nmin2_n, t_npls1_n, &
            & w_n, w_nmin1, x_nmin1, lz_conv, b, res_fast, &
            & t_nmin2_nmin2, t_nmin1_nmin1, v_nmin1, v_npls1)
  END IF
  !
END SUBROUTINE komega_shifted_qmr_sym_finalize

SUBROUTINE komega_shifted_qmr_sym_output_restart(Av_n, v_n, fo, filename)
  USE komega_parameter, ONLY : nl, nz, iter, itermax
  USE komega_vals_shifted_qmr_sym, ONLY: alpha_n, beta_n, beta_nmin1, c_n, c_nmin1, c_nmin2, g_n, &
          & g_npls1, p_n, p_nmin1, p_nmin2, r_nmin1, s_n, s_nmin1, s_nmin2, t_n_n, t_nmin1_n, t_nmin2_n, &
          & t_npls1_n, w_n, w_nmin1, x_nmin1, b, res_fast, t_nmin2_nmin2, t_nmin1_nmin1, &
          & v_nmin1, v_npls1,  z

  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(IN) :: Av_n(nl), v_n(nl)
  INTEGER, INTENT(IN) :: fo
  CHARACTER(*), INTENT(IN) :: filename

  WRITE(*,*)
  WRITE(*,*) "#####  Output Restart File  #####"
  WRITE(*,*)
  !
  OPEN(fo, file = filename, action = 'write')
  !
  WRITE(fo,*) iter
  WRITE(fo, *) alpha_n
  WRITE(fo, *) beta_n
  WRITE(fo, *) beta_nmin1
  WRITE(fo, *) c_n
  WRITE(fo, *) c_nmin1
  WRITE(fo, *) c_nmin2
  WRITE(fo, *) g_n
  WRITE(fo, *) g_npls1
  WRITE(fo, *) p_n
  WRITE(fo, *) p_nmin1
  WRITE(fo, *) p_nmin2
  WRITE(fo, *) r_nmin1
  WRITE(fo, *) s_n
  WRITE(fo, *) s_nmin1
  WRITE(fo, *) s_nmin2
  WRITE(fo, *) t_n_n
  WRITE(fo, *) t_nmin1_n
  WRITE(fo, *) t_nmin2_n
  WRITE(fo, *) t_npls1_n
  WRITE(fo, *) w_n
  WRITE(fo, *) w_nmin1
  WRITE(fo, *) x_nmin1
  WRITE(fo, *) b
  WRITE(fo, *) res_fast
  WRITE(fo, *) t_nmin2_nmin2
  WRITE(fo, *) t_nmin1_nmin1
  WRITE(fo, *) v_nmin1
  WRITE(fo, *) v_npls1
  WRITE(fo, *) Av_n
  WRITE(fo, *) v_n
  !
  close(fo)
  !
  WRITE(*,*) "  Restart File is written."
  !
END SUBROUTINE
!
SUBROUTINE komega_shifted_qmr_sym_init_restart(ndim0, nl0, nz0, x, z0, b0, itermax0, threshold0,&
        & Av_n, v_n, fi, filename, status)
  !
  USE komega_parameter, ONLY : lz_conv, nl, nz, iter, itermax, threshold, ndim
  USE komega_vals_shifted_qmr_sym, ONLY: alpha_n, beta_n, beta_nmin1, c_n, c_nmin1, c_nmin2, g_n, &
          & g_npls1, p_n, p_nmin1, p_nmin2, r_nmin1, s_n, s_nmin1, s_nmin2, t_n_n, t_nmin1_n, t_nmin2_n, &
          & t_npls1_n, w_n, w_nmin1, x_nmin1, b, res_fast, t_nmin2_nmin2, t_nmin1_nmin1, &
          & v_nmin1, v_npls1,  z
  !
  IMPLICIT NONE

  COMPLEX(8),INTENT(OUT) :: Av_n(nl), v_n(nl)
  INTEGER, INTENT(IN) :: fi
  CHARACTER(*), INTENT(IN) :: filename
  INTEGER,INTENT(IN) :: ndim0, nl0, nz0, itermax0
  REAL(8),INTENT(IN) :: threshold0
  COMPLEX(8),INTENT(IN) :: z0(nz0)
  COMPLEX(8),INTENT(IN) :: b0(nl0)
  COMPLEX(8),INTENT(OUT) :: x(nl0,nz0)
  INTEGER, INTENT(OUT) :: status(2)

  INTEGER :: ierr, iz
  LOGICAL :: flag_converge

  ndim = ndim0
  nl = nl0
  nz = nz0
  itermax = itermax0
  threshold = threshold0
  !
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Check Restart File  #####"
  WRITE(*,*)
  !
  open(fi, file = filename, status="old", action = 'read',iostat = ierr)
  !
  IF(ierr == 0) THEN
    !
    WRITE(*,*) "  Restart file is found."
    !
    ALLOCATE(z(nz), c_n(nz), c_nmin1(nz), c_nmin2(nz), g_n(nz), g_npls1(nz), p_n(nl, nz), p_nmin1(nl, nz),&
            & p_nmin2(nl, nz), r_nmin1(nz), s_n(nz), &
            & s_nmin1(nz), s_nmin2(nz), t_n_n(nz), t_nmin1_n(nz), t_nmin2_n(nz), t_npls1_n(nz), &
            & w_n(nl, nz), w_nmin1(nl, nz), x_nmin1(nl, nz), lz_conv(nz), b(nl), res_fast(nz, itermax), &
            & t_nmin2_nmin2(nz), t_nmin1_nmin1(nz), v_nmin1(nl), v_npls1(nl))
    z(:) = z0(:)
    READ(fi,*) iter
    READ(fi, *) alpha_n
    READ(fi, *) beta_n
    READ(fi, *) beta_nmin1
    READ(fi, *) c_n
    READ(fi, *) c_nmin1
    READ(fi, *) c_nmin2
    READ(fi, *) g_n
    READ(fi, *) g_npls1
    READ(fi, *) p_n
    READ(fi, *) p_nmin1
    READ(fi, *) p_nmin2
    READ(fi, *) r_nmin1
    READ(fi, *) s_n
    READ(fi, *) s_nmin1
    READ(fi, *) s_nmin2
    READ(fi, *) t_n_n
    READ(fi, *) t_nmin1_n
    READ(fi, *) t_nmin2_n
    READ(fi, *) t_npls1_n
    READ(fi, *) w_n
    READ(fi, *) w_nmin1
    READ(fi, *) x_nmin1
    READ(fi, *) b
    READ(fi, *) res_fast
    READ(fi, *) t_nmin2_nmin2
    READ(fi, *) t_nmin1_nmin1
    READ(fi, *) v_nmin1
    READ(fi, *) v_npls1
    READ(fi, *) Av_n
    READ(fi, *) v_n
    !
    close(fi)
    x(:,:) = x_nmin1(:,:)

    flag_converge = .TRUE.
    DO iz = 1,nz
      IF (lz_conv(iz) .EQV. .FALSE.) THEN
        flag_converge = .FALSE.
      END IF
    END DO

    IF (flag_converge .EQV. .TRUE.) THEN
      status(1) = -iter
      status(2) = 0
    ELSE IF (iter == itermax) THEN
      status(1) = iter
      status(2) = 1
    ELSE
      status(1) = iter
      status(2) = 0
    END IF
    !
  ELSE
    !
    WRITE(*,*) "  Restart file is NOT found."
    STOP
    !
  END IF
  !
END SUBROUTINE komega_shifted_qmr_sym_init_restart


SUBROUTINE komega_shifted_qmr_sym_getresidual(res)
  !
  USE komega_vals_shifted_qmr_sym, ONLY : r_nmin1
  USE komega_parameter, ONLY: nz
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(OUT) :: res(nz)
  !
  res(:) = r_nmin1(:)
  !
END SUBROUTINE komega_shifted_qmr_sym_getresidual

END MODULE komega_shifted_qmr_sym
