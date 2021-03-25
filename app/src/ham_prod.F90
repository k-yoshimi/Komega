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
MODULE ham_vals
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & ndiag,   & ! Diagonal components
  & nham,    & ! Non-zero elements of compressed Hamiltonian
  & nsite,   & ! Number of sites for the Heisenberg Chain
  & nsitep     ! Number of sites in inter-process region
  !
  REAL(8),SAVE :: &
  & Jx, & ! Heisenberg Jx
  & Jy, & ! Heisenberg Jy
  & Jz, & ! Heisenberg Jz
  & Dz    ! D-M interaction
  !
  LOGICAL,ALLOCATABLE,SAVE :: &
  & ud(:), & ! i=Up   & i+1=Down flag
  & du(:), & ! i=Down & i+1=Up   flag
  & para(:) ! i, i+1 Parallel flag
  !
  INTEGER,ALLOCATABLE,SAVE :: &
  & ham_indx(:,:), & ! row & column index of Hamiltonian
  & pair(:), &       ! pair connected by S+S- or S+S+
  & colind(:), &     ! index for column
  & rowptr(:)        ! row pointer
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & ham(:), &     ! Compressed Hamiltonian
  & dvalues(:), & ! Diagonal element of Hamiltonian
  & values(:)     ! Off-diagonal element of Hamiltonian
  !
  integer,allocatable,save :: row_ptr(:),col_ind(:),row_se(:,:)
  complex(kind(0.0D0)),allocatable,save :: ham_crs_val(:)
  real(8),allocatable,save :: ham_crs_val_r(:)
  !
#if defined(__MPI)
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & veci_buf(:)
#endif
  !
END MODULE ham_vals
!
! Module contains routines for Hamiltonian-Vector product
!
MODULE ham_prod_mod
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC ham_prod, print_ham, onthefly_init, finalize_ham
  !
CONTAINS
!
! Driver routine for the Hamiltonian-Vector Product
!
SUBROUTINE ham_prod(veci,veco,t11,t12,lcollect)
  !
  USE shiftk_vals, ONLY : ndim, inham, solver
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(IN) :: veci(ndim)
  COMPLEX(8),INTENT(OUT) :: veco(ndim)
  double precision,intent(OUT) :: t11,t12
  logical,intent(IN) :: lcollect
  !
  veco(1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
  !
  !WRITE(*,*) veci(1:10)
  !
  IF(inham == "") THEN
     CALL ham_prod_onthefly(veci,veco)
     t11 = 0.0D0
     t12 = 0.0D0
  ELSE
     IF (TRIM(solver) == 'bicg') THEN
     !CALL ham_prod_compress_csr(veci,veco,t11,t12,lcollect)
        call ham_prod_compress_crs(veci,veco,t11,t12,lcollect)
     ELSE
        call ham_prod_compress_crs_r(veci,veco,t11,t12,lcollect)
     ENDIF
  END IF
  !
END SUBROUTINE ham_prod
!
! Hamiltonian-vector product with the Compressed representation
!
SUBROUTINE ham_prod_compress(veci,veco)
  !
  USE shiftk_vals, ONLY : ndim
  USE ham_vals, ONLY : nham, ndiag, ham, ham_indx
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(IN) :: veci(ndim)
  COMPLEX(8),INTENT(OUT) :: veco(ndim)
  !
  INTEGER :: iham
  !
  DO iham = 1, ndiag
     veco(ham_indx(1,iham)) = veco(ham_indx(1,iham)) &
     &          + ham(iham) * veci(ham_indx(2,iham))
  END DO
  !
  DO iham = ndiag + 1, nham
     veco(ham_indx(1,iham)) = veco(ham_indx(1,iham)) &
     &          + ham(iham) * veci(ham_indx(2,iham))
     veco(ham_indx(2,iham)) = veco(ham_indx(2,iham)) &
     &   + CONJG(ham(iham)) * veci(ham_indx(1,iham))
  END DO
  !
END SUBROUTINE ham_prod_compress
!
! Hamiltonian-vector product with CSR format (Gkountouvas et al. (2013))
!
SUBROUTINE ham_prod_compress_csr(veci,veco,t11,t12,lcollect)
  !$ USE omp_lib, only : OMP_GET_NUM_THREADS, omp_get_wtime
  USE shiftk_vals, ONLY : ndim, myrank, nproc
  USE ham_vals, ONLY : nham, ndiag, ham, ham_indx, colind, rowptr, dvalues, values
#if defined(__MPI)
  USE mpi, only : MPI_COMM_WORLD, MPI_DOUBLE_COMPLEX, MPI_SUM
#endif
  !
  IMPLICIT NONE
  !
#if defined(__MPI)
  INTEGER :: ierr, iend
#endif
  !
  COMPLEX(8),INTENT(IN) :: veci(ndim)
  COMPLEX(8),INTENT(OUT) :: veco(ndim)
  double precision,intent(OUT) :: t11,t12
  logical,intent(IN) :: lcollect
  !
#if defined(__MPI)
  COMPLEX(8) :: vecbuf(ndim)
#else
  COMPLEX(8),ALLOCATABLE :: vecbuf(:,:)
#endif
  INTEGER, ALLOCATABLE :: is(:), ie(:)
  COMPLEX(8) :: zero, vtmp
  !DOUBLE PRECISION :: DOT_PRODUCT
  INTEGER :: i, ib, iham, jham, col, ham_indx0(2), nthreads
  !
#if defined(__MPI)
  vecbuf = CMPLX(0d0, 0d0, KIND(0d0))
  ib = ndim/nproc
  IF(myrank == nproc-1) THEN
     iend = ndim
  ELSE
     iend = (myrank+1)*ib
  ENDIF
  DO iham = myrank*ib+1, iend
     vecbuf(iham) = vecbuf(iham) + dvalues(iham) * veci(iham)
     DO jham = rowptr(iham), rowptr(iham+1)-1
        col = colind(jham)
        vecbuf(iham) = vecbuf(iham) + values(jham) * veci(col)
        vecbuf(col) = vecbuf(col) + CONJG(values(jham)) * veci(iham)
     END DO
  END DO
  CALL MPI_ALLREDUCE(vecbuf, veco, ndim, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  nthreads = 1
  !$OMP PARALLEL
  !$OMP MASTER
  !$ nthreads = OMP_GET_NUM_THREADS()
  !$OMP END MASTER
  !$OMP END PARALLEL
  ALLOCATE(vecbuf(nthreads,ndim),is(nthreads),ie(nthreads))
  zero = CMPLX(0d0, 0d0, KIND(0d0))
  vecbuf = zero
  veco = zero
  !
  ib = ndim/nthreads
  DO i = 1, nthreads
    is(i) = (i-1)*ib+1
    ie(i) = i*ib
  END DO
  ie(nthreads) = ndim
  !WRITE(*,*) nthreads
  !WRITE(*,*) ndim
  !STOP
  !$ t11 = omp_get_wtime()
  !if(lcollect) call start_collection("region11")
  !$OMP PARALLEL default(shared), private(i, iham, jham, col, vtmp), firstprivate(ndim, nthreads, veci)
  !$OMP DO
  DO i = 1, nthreads
     DO iham = is(i), ie(i)
        vtmp = dvalues(iham) * veci(iham)
        DO jham = rowptr(iham), rowptr(iham+1)-1
           col = colind(jham)
           vtmp = vtmp + values(jham) * veci(col)
           vecbuf(i,col) = vecbuf(i,col) + CONJG(values(jham)) * veci(iham)
        END DO
        vecbuf(i,iham) = vtmp
     END DO
  END DO
  !$OMP END DO
  !$OMP DO
  DO iham = 1, ndim
     vtmp = zero
     DO i = 1, nthreads
        vtmp = vtmp + vecbuf(i,iham)
     END DO
     veco(iham) = vtmp
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  !if(lcollect) call stop_collection("region11")
  !$ t12 = omp_get_wtime()
  DEALLOCATE(vecbuf,is,ie)
#endif
  !
END SUBROUTINE ham_prod_compress_csr
  !
subroutine ham_prod_compress_crs(veci,veco,t11,t12,lcollect)
  !$ use omp_lib
  use shiftk_vals,only : ndim
  use ham_vals,only : row_ptr,col_ind,ham_crs_val,row_se
  implicit none
  complex(kind(0.0D0)),intent(IN ) :: veci(ndim)
  complex(kind(0.0D0)),intent(OUT) :: veco(ndim)
  double precision,intent(OUT) :: t11,t12
  logical,intent(IN) :: lcollect
  !
  integer :: i,j,ithread,is,ie
  complex(kind(0.0D0)) :: czero
  !
#if defined(__MPI)
  write(*,*) "MPI version of subroutine ham_prod_compress_crs has not been made."
  stop
#else
  !
  !$ t11 = omp_get_wtime()
  !if(lcollect) call start_collection("region11")
  ithread = 0
  !$OMP PARALLEL default(shared), private(i,j,ithread,is,ie,czero)
  !$ ithread = omp_get_thread_num()
  is = row_se(1,ithread)
  ie = row_se(2,ithread)
  czero = cmplx(0.0D0,0.0D0,kind(0.0D0))
  !write(6,*) "is ie = ", is, ie
  do i = is, ie
     !write(6,*) "i", i
     veco(i) = czero
     !write(6,*) "start loop row_ptr"
     do j = row_ptr(i), row_ptr(i+1) - 1
        !write(6,*) "j", j
        !write(6,*) "col_ind(j)", col_ind(j)
        veco(i) = veco(i) + ham_crs_val(j) * veci(col_ind(j))
     end do
  end do
  !$OMP END PARALLEL
  !if(lcollect) call stop_collection("region11")
  !$ t12 = omp_get_wtime()
  !
#endif
end subroutine ham_prod_compress_crs
!
subroutine ham_prod_compress_crs_r(veci,veco,t11,t12,lcollect)
  !$ use omp_lib
  use shiftk_vals,only : ndim
  use ham_vals,only : row_ptr,col_ind,ham_crs_val_r,row_se
  implicit none
  complex(kind(0.0D0)),intent(IN ) :: veci(ndim)
  complex(kind(0.0D0)),intent(OUT) :: veco(ndim)
  double precision,intent(OUT) :: t11,t12
  logical,intent(IN) :: lcollect
  !
  integer :: i,j,ithread,is,ie
  complex(kind(0.0D0)) :: czero
  !
#if defined(__MPI)
  write(*,*) "MPI version of subroutine ham_prod_compress_crs_r has not been made."
  stop
#else
  !
  !$ t11 = omp_get_wtime()
  !if(lcollect) call start_collection("region11")
  ithread = 0
  !$OMP PARALLEL default(shared), private(i,j,ithread,is,ie,czero)
  !$ ithread = omp_get_thread_num()
  is = row_se(1,ithread)
  ie = row_se(2,ithread)
  czero = cmplx(0.0D0,0.0D0,kind(0.0D0))
  do i = is, ie
     veco(i) = czero
     do j = row_ptr(i), row_ptr(i+1) - 1
        veco(i) = veco(i) + ham_crs_val_r(j) * veci(col_ind(j))
     end do
  end do
  !$OMP END PARALLEL
  !if(lcollect) call stop_collection("region11")
  !$ t12 = omp_get_wtime()
  !
#endif
end subroutine ham_prod_compress_crs_r
!
! Hamiltonian-vector product with On-The-Fly Hamiltonian generation
!
SUBROUTINE ham_prod_onthefly(veci,veco)
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_COMM_WORLD, MPI_DOUBLE_COMPLEX, MPI_STATUS_SIZE
  USE shiftk_vals, ONLY : myrank
  USE ham_vals, ONLY : veci_buf
#endif
  USE shiftk_vals, ONLY : ndim, almost0
  USE ham_vals, ONLY : Jx, Jy, Jz, Dz, nsite, ud, du, para, pair, nsitep
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(IN) :: veci(ndim)
  COMPLEX(8),INTENT(OUT) :: veco(ndim)
  !
  INTEGER :: isite, isite1, mask1, mask2, mask12, spin, idim, nsite2
  COMPLEX(8) :: matrix, cmatrix
#if defined(__MPI)
  INTEGER :: origin, status(MPI_STATUS_SIZE), ierr
#endif
  !
  ! Intra process region
  !
  IF(nsitep == 0) THEN
     nsite2 = nsite
  ELSE
     nsite2 = nsite - nsitep - 1
  END IF
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nsite2, nsite, ndim, pair, ud, du, para, Jx, Jy, Jz, Dz, veci, veco) &
  !$OMP & PRIVATE(isite, isite1, idim, spin, mask1, mask2, mask12, &
  !$OMP &         matrix, cmatrix)
  !
  DO isite = 1, nsite2
     !
     isite1 = MOD(isite, nsite) + 1
     !
     mask12 = 0
     mask12 = IBSET(mask12, isite  - 1)
     mask12 = IBSET(mask12, isite1 - 1)
     !
     mask1 = 0
     mask1 = IBSET(mask1, isite  - 1)
     !
     mask2 = 0
     mask2 = IBSET(mask2, isite1 - 1)
     !
     !$OMP DO
     DO idim = 1, ndim
        ud(idim) = .FALSE.
        du(idim) = .FALSE.
        para(idim) = .FALSE.
     END DO
     !$OMP END DO NOWAIT
     !
     !$OMP DO
     DO idim = 1, ndim
        !
        spin = IAND(idim - 1, mask12)
        pair(idim) = IEOR(idim - 1, mask12) + 1
        !
        IF(spin == mask1) THEN
           ud(idim) = .TRUE.
        ELSE IF(spin == mask2) THEN
           du(idim) = .TRUE.
        ELSE
           para(idim) = .TRUE.
        END IF
        !
     END DO
     !$OMP END DO NOWAIT
     !
     ! S_{i z} S_{i+1 z}
     !
     IF(ABS(Jz) > almost0) THEN
        !
        matrix = CMPLX(0.25d0 * Jz, 0d0, KIND(0d0))
        !$OMP DO
        DO idim = 1, ndim
           IF(para(idim)) THEN
              veco(idim) = veco(idim) + matrix * veci(idim)
           ELSE
              veco(idim) = veco(idim) - matrix * veci(idim)
           END IF
        END DO
        !$OMP END DO NOWAIT
        !
     END IF        
     !
     ! S_{i}^+ S_{i+1}^- + S_{i}^- S_{i+1}^+
     !
     IF(ABS(Jx + Jy) > almost0 .OR. ABS(Dz) > almost0) THEN
        !
        matrix = CMPLX(0.25d0 * (Jx + Jy), 0.5d0 * Dz, KIND(0d0))
        cmatrix = CONJG(matrix)
        !$OMP DO
        DO idim = 1, ndim
           IF(du(idim)) THEN
              veco(idim) = veco(idim) +  matrix * veci(pair(idim))
           ELSE IF(ud(idim)) THEN
              veco(idim) = veco(idim) + cmatrix * veci(pair(idim))
           END IF
        END DO
        !$OMP END DO NOWAIT
        !
     END IF
     !
     ! S_{i}^+ S_{i+1}^+ + S_{i}^- S_{i+1}^
     !
     IF(ABS(Jx - Jy) > almost0) THEN
        !
        matrix = CMPLX(0.25d0 * (Jx - Jy), 0d0, KIND(0d0)) 
        !$OMP DO
        DO idim = 1, ndim
           IF(para(idim)) veco(idim) = veco(idim) + matrix * veci(pair(idim))
        END DO
        !$OMP END DO NOWAIT
        !
     END IF
     !
  END DO ! isite = 1, nsite
  !
  !$OMP END PARALLEL
  !
  ! isite = nsite - nsitep .OR. isite = nsite
  !
#if defined(__MPI)
  IF(nsitep > 0) THEN
     !
     DO isite = 1, 2
        !
        IF(isite == 1) THEN
           mask1 = 0
           mask1 = IBSET(mask1, nsite - nsitep - 1)
           mask2 = 0
           mask2 = IBSET(mask2, 0)
        ELSE
           mask1 = 0
           mask1 = IBSET(mask1, 0)
           mask2 = 0
           mask2 = IBSET(mask2, nsitep - 1)
        END IF
        origin = IEOR(myrank, mask2)
        !
        CALL mpi_sendrecv(veci, ndim, MPI_DOUBLE_COMPLEX, origin, 0, &
        &             veci_buf, ndim, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, status, ierr);
        !
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP & SHARED(ndim, mask1, myrank, mask2, para, pair, Jx, Jy, Jz, Dz, &
        !$OMP &        veci, veco, veci_buf, isite) &
        !$OMP & PRIVATE(idim, spin, matrix)
        !
        !$OMP DO
        DO idim = 1, ndim
           !
           spin = IAND(idim - 1, mask1)
           pair(idim) = IEOR(idim - 1, mask1) + 1
           !
           IF(spin == mask1) THEN
              para(idim) = .TRUE.
           ELSE
              para(idim) = .FALSE.
           END IF
           !
        END DO
        !$OMP END DO NOWAIT
        !
        IF(IAND(myrank, mask2) == 0) THEN
           !$OMP DO
           DO idim = 1, ndim
              para(idim) = .NOT. para(idim)
           END DO
           !$OMP END DO NOWAIT
        END IF
        !
        ! S_{i z} S_{i+1 z}
        !
        IF(ABS(Jz) > almost0) THEN
           !
           matrix = CMPLX(0.25d0 * Jz, 0d0, KIND(0d0))
           !$OMP DO
           DO idim = 1, ndim
              IF(para(idim)) THEN
                 veco(idim) = veco(idim) + matrix * veci(idim)
              ELSE
                 veco(idim) = veco(idim) - matrix * veci(idim)
              END IF
           END DO
           !$OMP END DO NOWAIT
           !
        END IF
        !
        ! S_{i}^+ S_{i+1}^- + S_{i}^- S_{i+1}^+
        !
        IF(ABS(Jx + Jy) > almost0 .OR. ABS(Dz) > almost0) THEN
           !
           IF((IAND(myrank, mask2) == mask2 .AND. isite == 1) .OR. &
           &  (IAND(myrank, mask2) == 0     .AND. isite == 2) ) THEN
              matrix = CMPLX(0.25d0 * (Jx + Jy), 0.5d0 * Dz, KIND(0d0))
           ELSE
              matrix = CMPLX(0.25d0 * (Jx + Jy), - 0.5d0 * Dz, KIND(0d0))
           END IF
           !
           !$OMP DO
           DO idim = 1, ndim
              IF(.NOT. para(idim)) THEN
                 veco(idim) = veco(idim) + matrix * veci_buf(pair(idim))
              END IF
           END DO
           !$OMP END DO NOWAIT
           !
        END IF
        !
        ! S_{i}^+ S_{i+1}^+ + S_{i}^- S_{i+1}^
        !
        IF(ABS(Jx - Jy) > almost0) THEN
           !
           matrix = CMPLX(0.25d0 * (Jx - Jy), 0d0, KIND(0d0)) 
           !$OMP DO
           DO idim = 1, ndim
              IF(para(idim)) veco(idim) = veco(idim) + matrix * veci_buf(pair(idim))
           END DO
           !$OMP END DO NOWAIT
           !
        END IF
        !
        !$OMP END PARALLEL
        !
     END DO ! isite = 1, 2
     !
  END IF ! (nsitep > 0)
  !
  ! nsite - nsitep < isite < nsite
  !
  DO isite = 1, nsitep - 1
     !
     isite1 = isite + 1
     !
     mask12 = 0
     mask12 = IBSET(mask12, isite  - 1)
     mask12 = IBSET(mask12, isite1 - 1)
     !
     mask1 = 0
     mask1 = IBSET(mask1, isite  - 1)
     !
     mask2 = 0
     mask2 = IBSET(mask2, isite1 - 1)
     !
     origin = IEOR(myrank, mask12)
     !
     CALL mpi_sendrecv(veci, ndim, MPI_DOUBLE_COMPLEX, origin, 0, &
     &             veci_buf, ndim, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, status, ierr);
     !
     spin = IAND(myrank, mask12)
     !
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(Jx, Jy, Jz, Dz, spin, mask1, mask2, mask12, veci, veco, veci_buf, ndim) &
     !$OMP & PRIVATE(idim, matrix)
     !
     ! S_{i z} S_{i+1 z}
     !
     IF(ABS(Jz) > almost0) THEN
        !
        IF(spin == mask12 .OR. spin == 0) THEN
           matrix = CMPLX(0.25d0 * Jz, 0d0, KIND(0d0))
        ELSE
           matrix = CMPLX(- 0.25d0 * Jz, 0d0, KIND(0d0))
        END IF
        !
        !$OMP DO
        DO idim = 1, ndim
           veco(idim) = veco(idim) + matrix * veci(idim)
        END DO
        !$OMP END DO NOWAIT
        !
     END IF        
     !
     ! S_{i}^+ S_{i+1}^- + S_{i}^- S_{i+1}^+
     !
     IF(ABS(Jx + Jy) > almost0 .OR. ABS(Dz) > almost0) THEN
        !
        IF(spin == mask1 .OR. spin == mask2) THEN
           !
           IF(spin == mask1) THEN
              matrix = CMPLX(0.25d0 * (Jx + Jy), - 0.5d0 * Dz, KIND(0d0))
           ELSE
              matrix = CMPLX(0.25d0 * (Jx + Jy),   0.5d0 * Dz, KIND(0d0))
           END IF
           !
           !$OMP DO
           DO idim = 1, ndim
              veco(idim) = veco(idim) + matrix * veci_buf(idim)
           END DO
           !$OMP END DO NOWAIT
           !
        END IF
        !
     END IF
     !
     ! S_{i}^+ S_{i+1}^+ + S_{i}^- S_{i+1}^
     !
     IF(ABS(Jx - Jy) > almost0) THEN
        !
        IF(spin == mask12 .OR. spin == 0) THEN
           !
           matrix = CMPLX(0.25d0 * (Jx - Jy), 0d0, KIND(0d0)) 
           !$OMP DO
           DO idim = 1, ndim
              veco(idim) = veco(idim) + matrix * veci_buf(idim)
           END DO
           !$OMP END DO NOWAIT
           !
        END IF
        !
     END IF
     !
     !$OMP END PARALLEL
     !
  END DO ! isite = 1, nsitep - 1
#endif
  !
END SUBROUTINE ham_prod_onthefly
!
! Allocate Flags and Pair index for On-The-Fly
!
SUBROUTINE onthefly_init()
  !
  USE shiftk_vals, ONLY : ndim
  USE ham_vals, ONLY : ud, du, para, pair
#if defined(__MPI)
  USE ham_vals, ONLY : veci_buf
#endif  
  !
  IMPLICIT NONE
  !
  ALLOCATE(ud(ndim), du(ndim), para(ndim), pair(ndim))
#if defined(__MPI)
  ALLOCATE(veci_buf(ndim))
#endif
  !
END SUBROUTINE onthefly_init
!
! Deallocate Hamiltonian parameters
!
SUBROUTINE finalize_ham()
  !
  USE shiftk_vals, ONLY : inham
  USE ham_vals, ONLY : ud, du, para, pair, ham, ham_indx
#if defined(__MPI)
  USE ham_vals, ONLY : veci_buf
#endif  
  !
  IMPLICIT NONE
  !
  IF(inham == "") THEN
     DEALLOCATE(ud, du, para, pair)
#if defined(__MPI)
     DEALLOCATE(veci_buf)
#endif
  ELSE
     !DEALLOCATE(ham, ham_indx)
  END IF
  !
END SUBROUTINE finalize_ham
!
! Print Hamiltonian with the Matrix Market format
!
SUBROUTINE print_ham()
  !
#if defined(__MPI)
  use mpi, ONLY : MPI_COMM_WORLD, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM
#endif
  USE shiftk_vals, ONLY : ndim, almost0, myrank, nproc
  USE ham_vals, ONLY : nham
  !
  IMPLICIT NONE
  !
  INTEGER :: idim, jdim, fo = 21, iproc, jproc
#if defined(__MPI)
  INTEGER :: ierr
#endif
  COMPLEX(8),ALLOCATABLE :: veci(:), veco(:)
  double precision :: t11,t12
  logical :: lcollect=.FALSE.
  !
  IF(nproc /= 1) RETURN
  !
  ALLOCATE(veci(ndim), veco(ndim))
  !
  OPEN(fo, file = "zvo_Ham.dat")
  IF(myrank == 0) WRITE(fo,*) "%%MatrixMarket matrix coordinate complex hermitian"
  !
  nham = 0
  DO iproc = 0, nproc - 1
     DO idim = 1, ndim
        !
        IF(myrank == iproc) THEN
           veci(1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
           veci(  idim) = CMPLX(1d0, 0d0, KIND(0d0))
        ELSE
           veci(1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
        END IF
        !
        CALL ham_prod(veci, veco, t11, t12, lcollect)
        !
        IF(myrank >= iproc) nham = nham + COUNT(ABS(veco(idim:ndim)) > almost0)
        !
     END DO ! idim = 1, ndim
  END DO ! iproc = 1, nproc
  !
#if defined(__MPI)
  call MPI_allREDUCE(MPI_IN_PLACE,nham, 1, &
  &                  MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  IF(myrank == 0) WRITE(fo,'(3i10)') ndim*nproc, ndim*nproc, nham 
  !
  DO iproc = 0, nproc - 1
     DO idim = 1, ndim
        !
        IF(myrank == iproc) THEN
           veci(1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
           veci(  idim) = CMPLX(1d0, 0d0, KIND(0d0))
        ELSE
           veci(1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
        END IF
        !
        CALL ham_prod(veci, veco, t11, t12, lcollect)
        !
        DO jproc = iproc, nproc - 1
           DO jdim = idim, ndim
              IF(ABS(veco(jdim)) > almost0 .AND. myrank == jproc) &
              &  WRITE(fo,'(2i10,2f15.8)') jdim + ndim * myrank, idim + ndim *iproc, &
              & DBLE(veco(jdim)), AIMAG(veco(jdim))
           END DO
#if defined(__MPI)
           CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
        END DO
        !
     END DO ! idim = 1, ndim
  END DO ! iproc = 1, nproc
  !
  CLOSE(fo)
  !
  DEALLOCATE(veci, veco)
  !
END SUBROUTINE print_ham
!
END MODULE ham_prod_mod
