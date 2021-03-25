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
! Input/Output ROUTINES
!
MODULE shiftk_io
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Initialization of MPI
!
SUBROUTINE shiftk_init()
  !
#if defined(__MPI)
  USE mpi, only : MPI_COMM_WORLD
#endif
  !$ USE omp_lib, only : OMP_GET_NUM_THREADS
  USE shiftk_vals, ONLY : myrank, nproc, stdout, inpunit
  !
  IMPLICIT NONE
  !
  CHARACTER(256) :: fname
  INTEGER ierr, iargc
  !
  !
  if(iargc() /= 1) then
     write(*,*) "Argument error. Usage: "
     write(*,*) "ShiftK.out namelist.def"
     stop
  end if
  !
#if defined(__MPI)
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, nproc, ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, myrank, ierr)
#else
  nproc = 1
  myrank = 0
#endif
  !
  IF(myrank == 0) THEN
     !
     CALL system("mkdir -p output")
     inpunit = 5
     stdout = 6
     !
     call getarg(1, fname)
     WRITE(*,*) "fname ", TRIM(fname)
     OPEN(inpunit, file = TRIM(fname), status = "OLD", iostat = ierr)
     !
     IF(ierr /= 0) THEN
        WRITE(*,*) "Cannot open input file ", TRIM(fname)
#if defined(__MPI)
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
        STOP
     ELSE
        WRITE(*,*) "  Open input file ", TRIM(fname)
     END IF
     !
  ELSE
     stdout = 6
     OPEN(stdout, file='/dev/null', status='unknown')
  END IF
  !
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "  Number of processes : ", nproc
  !$OMP PARALLEL
  !$OMP MASTER
  !$ WRITE(stdout,*) "  Number of threads : ", OMP_GET_NUM_THREADS()
  !$OMP END MASTER
  !$OMP END PARALLEL
  WRITE(stdout,*)
  !
  !
END SUBROUTINE shiftk_init
!
! Filname for Hamiltonian and RHS vector
!
SUBROUTINE input_filename()
  !
#if defined(__MPI)
  USE mpi, only : MPI_COMM_WORLD, MPI_CHARACTER
#endif
  USE shiftk_vals, ONLY : inham, invec, stdout, myrank, inpunit, solver
  !
  IMPLICIT NONE
  CHARACTER(20) hamtype
  !
#if defined(__MPI)
  INTEGER ierr
#endif
  NAMELIST /filename/ inham, invec, hamtype
  !
  inham = ""
  invec = ""
  hamtype = "complex"
  !
  IF(myrank == 0) READ(inpunit,filename,err=100)
  !
  IF(TRIM(hamtype) == "real") THEN
      solver = "cocg"
  ELSE IF(TRIM(hamtype) == "complex") THEN
      solver = "bicg"
  ELSE
      goto 100
  ENDIF
  !
#if defined(__MPI)
  call MPI_BCAST(inham, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(invec, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
#endif
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input FileName  ##########"
  WRITE(stdout,*)
  WRITE(stdout,*) "  Hamiltonian file : ", TRIM(inham)
  WRITE(stdout,*) "  Excited state file : ", TRIM(invec)
  !
  return
  !
100 write(*,*) "Stop in INPUT_FILENAME. reading namelist FILENAME"
  !
#if defined(__MPI)
     CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  STOP
  !
END SUBROUTINE input_filename
!
! Read parameter for the On-The-Fly Hamiltonian generation
!
SUBROUTINE input_parameter_ham()
  !
#if defined(__MPI)
  USE mpi, only : MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_PRECISION
#endif
  USE shiftk_vals, ONLY : ndim, almost0, solver, stdout, myrank, nproc, inpunit
  USE ham_vals, ONLY : Jx, Jy, Jz, Dz, nsite, nsitep
  !
  IMPLICIT NONE
  !
#if defined(__MPI)
  INTEGER ierr
#endif
  NAMELIST /ham/ Jx, Jy, Jz, Dz, nsite
  !
  Jx = 1d0
  Jy = 1d0
  Jz = 1d0
  Dz = 0d0
  nsite = 4
  !
  IF(myrank == 0) READ(inpunit,ham,err=100)
  !
#if defined(__MPI)
  call MPI_BCAST(nsite, 1, MPI_INTEGER,       0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(Jx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(Jy, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(Jz, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(Dz, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input Parameter for Hamiltonian ##########"
  WRITE(stdout,*)
  WRITE(stdout,*) "    TOTAL Number of sites : ", nsite
  !
  nsitep = NINT(LOG(DBLE(nproc)) / LOG(2d0))
  IF(2**nsitep /= nproc) THEN
     WRITE(*,*) "ERROR ! Number of processes is not 2-exponent."
#if defined(__MPI)
     CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
     STOP
  END IF
  WRITE(stdout,*) "    LOCAL Number of sites : ", nsite - nsitep
  !
  WRITE(stdout,*) "                       Jx : ", Jx
  WRITE(stdout,*) "                       Jy : ", Jy
  WRITE(stdout,*) "                       Jz : ", Jz
  WRITE(stdout,*) "                       Dz : ", Dz
  ndim = 2**(nsite - nsitep)
  WRITE(stdout,*) "  Dim. of Hamiltonian : ", ndim
  !
  ! Hermitian(BiCG) or Real-Symmetric(COCG)
  !
!  IF(ABS(Dz) > almost0) THEN
!     WRITE(stdout,*) "  BiCG mathod is used."
!!     WRITE(stdout,*) "  shifted_qmr_sym mathod is used."
!     solver = "shifted_qmr_sym"
!  ELSE
!     WRITE(stdout,*) "  COCG mathod is used."
!     solver = "COCG"
!  END IF
  !
  return
  !
100 write(*,*) "Stop in INPUT_PARAMETER for Hamiltonian. reading namelist HAM"
  !
#if defined(__MPI)
  CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  STOP
  !
END SUBROUTINE input_parameter_ham
!
! Parameter for All CG calculations
!
SUBROUTINE input_parameter_cg()
  !
#if defined(__MPI)
  USE mpi, only : MPI_COMM_WORLD, MPI_INTEGER
#endif
  USE shiftk_vals, ONLY : maxloops, threshold, ndim, stdout, myrank, inpunit, solver
  !
  IMPLICIT NONE
  !
  INTEGER :: convfactor
  CHARACTER(20) :: method
#if defined(__MPI)
  INTEGER ierr
#endif
  NAMELIST /cg/ method, maxloops, convfactor
  !WRITE(*,*) "maxloops", maxloops
  !WRITE(*,*) "convfactor", convfactor
  !
  maxloops = ndim
  convfactor = 8
  method = ""
  !
  IF(myrank == 0) READ(inpunit,cg,err=100)
  IF(TRIM(method) /= "") solver = method
  !
#if defined(__MPI)
  call MPI_BCAST(maxloops,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(convfactor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(solver, 20, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
#endif
  !
  threshold = 10d0**(-convfactor)
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input Parameter for CG Iteration ##########"
  WRITE(stdout,*)
  WRITE(stdout,*) "                  Method : ", TRIM(solver)
  WRITE(stdout,*) "  Maximum number of loop : ", maxloops
  WRITE(stdout,*) "   Convergence Threshold : ", threshold
  !
  return
  !
100 write(*,*) "Stop in INPUT_PARAMETER for CG. reading namelist CG"
  !
#if defined(__MPI)
  CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  STOP
  !
END SUBROUTINE input_parameter_cg
!
! Read Parameter for the Spectrum calculation
!
SUBROUTINE input_parameter_dyn()
  !
#if defined(__MPI)
  USE mpi, only : MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_COMPLEX, &
  &               MPI_CHARACTER, MPI_LOGICAL
#endif
  USE shiftk_vals, ONLY : nomega, outrestart, calctype, z, e_max, e_min, stdout, myrank, inpunit
  !
  IMPLICIT NONE
  !
  INTEGER :: iomega
#if defined(__MPI)
  INTEGER ierr
#endif
  COMPLEX(8) :: omegamax, omegamin
  NAMELIST /dyn/ omegamax, omegamin, nomega, calctype, outrestart
  !
  calctype = "normal"
  nomega = 10
  omegamin = CMPLX(e_min, 0.01d0*(e_max - e_min), KIND(0d0))
  omegamax = CMPLX(e_max, 0.01d0*(e_max - e_min), KIND(0d0))
  !omegamin = CMPLX(e_min, 0.d0, KIND(0d0))
  !omegamax = CMPLX(e_max, 0.d0, KIND(0d0))
  outrestart = .FALSE.
  !
  IF(myrank == 0) READ(inpunit,dyn,err=100)
  !
#if defined(__MPI)
  call MPI_BCAST(nomega,   1, MPI_INTEGER,        0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(omegamin, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(omegamax, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(calctype, 256, MPI_CHARACTER,    0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(outrestart, 1, MPI_LOGICAL,      0, MPI_COMM_WORLD, ierr)
#endif
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input Parameter for Spectrum  ##########"
  WRITE(stdout,*)
  WRITE(stdout,*) "           Max. of Omega : ", omegamax
  WRITE(stdout,*) "           Min. of Omega : ", omegamin
  WRITE(stdout,*) "           Num. of Omega : ", nomega
  WRITE(stdout,'(a,a)') "         Calculation type : ", calctype
  !
  ALLOCATE(z(nomega))
  z(1) = omegamin
  DO iomega = 2, nomega
     z(iomega) = omegamin + (omegamax - omegamin) * DBLE(iomega - 1) / DBLE(nomega - 1)
  END DO
  !
  return
  !
100 write(*,*) "Stop in INPUT_PARAMETER_DYN. reading namelist DYN"
  !
#if defined(__MPI)
  CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  STOP
  !
END SUBROUTINE input_parameter_dyn
!
! Input Hamiltonian with the Matrix Market Format
!
SUBROUTINE input_hamiltonian()
  !$ use omp_lib
  !
#if defined(__MPI)
  USE mpi, only : MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_COMPLEX, MPI_WTIME
#endif
  USE shiftk_vals, ONLY : ndim, inham, solver, almost0, nproc, stdout, myrank
  USE ham_vals, ONLY : ham, ham_indx, nham, ndiag, colind, rowptr, dvalues, values
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ndim2, iham, jham, ham_indx0(2)
#if defined(__MPI)
  INTEGER :: ierr
  DOUBLE PRECISION :: t1, t2
#endif
  double precision :: t01, t02
  REAL(8) :: ham_r, ham_i
  COMPLEX(8) :: ham0
  CHARACTER(100) :: ctmp
  INTEGER,ALLOCATABLE :: rowcount(:)
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input Hamiltonian  ##########"
  WRITE(stdout,*)
  !
#if defined(__MPI)
  t1 = MPI_WTIME()
  IF(myrank == 0) THEN
#endif
  !
  !$ t01 = omp_get_wtime()
  OPEN(fi, file = TRIM(inham))
  READ(fi, *) ctmp
  READ(fi,*) ndim, ndim2, nham
  WRITE(stdout,*) "          Dim. of Hamiltonian : ", ndim, ndim2
  WRITE(stdout,*) "  Num. of Non-Zero Components : ", nham
  !
  IF(ndim2 /= ndim) THEN
     WRITE(stdout,*) "ERROR ! Hamiltonian is not square."
     STOP
  END IF
  !
  ALLOCATE(ham(nham), ham_indx(2,nham), rowcount(ndim), rowptr(ndim+1), dvalues(ndim))
  !
  DO iham = 1, nham
     READ(fi,*) ham_indx(1,iham), ham_indx(2,iham), ham_r, ham_i
     ham(iham) = CMPLX(ham_r, ham_i, KIND(0d0))
  END DO
  !
  ! Count the number of the Diagonal components
  ! Sort: Diagonal component should be first
  !
  ndiag = 0
  DO iham = 1, nham
     !
     IF(ham_indx(1,iham) == ham_indx(2,iham)) THEN
        !
        ndiag = ndiag + 1
        !
        ham_indx0(1:2) = ham_indx(1:2,iham)
        ham_indx(1:2,iham) = ham_indx(1:2,ndiag)
        ham_indx(1:2,ndiag) = ham_indx0(1:2)
        !
        ham0 = ham(iham)
        ham(iham) = ham(ndiag)
        ham(ndiag) = ham0
        !
     END IF
     !
  END DO
  !
  ALLOCATE(values(nham-ndiag), colind(nham-ndiag))
  !
  CALL QSORT(ham_indx(1,:), ham_indx(2,:), ham, ndiag+1, nham)
!  DO iham = ndiag + 1, nham
!     !
!     DO jham = iham + 1, nham
!        IF(ham_indx(1,iham) > ham_indx(1,jham)) THEN
!           !
!           ham_indx0(1:2) = ham_indx(1:2,jham)
!           ham_indx(1:2,jham) = ham_indx(1:2,iham)
!           ham_indx(1:2,iham) = ham_indx0(1:2)
!           !
!           ham0 = ham(jham)
!           ham(jham) = ham(iham)
!           ham(iham) = ham0
!           !   
!        END IF
!     END DO
!     !
!  END DO
  !
  DO iham = 1, ndim
     rowcount(iham) = 0
     dvalues(iham) = CMPLX(0)
  END DO
  !
  DO iham = 1, ndiag
     dvalues(ham_indx(1,iham)) = ham(iham)
  END DO
  !
  DO iham = ndiag + 1, nham
     IF(ham_indx(1,iham) > ham_indx(2,iham)) THEN
        rowcount(ham_indx(1,iham)) = rowcount(ham_indx(1,iham)) + 1
        colind(iham-ndiag) = ham_indx(2,iham)
        values(iham-ndiag) = ham(iham)
     ENDIF
  END DO
  !
  rowptr(1) = 1
  DO iham = 2, ndim+1
     rowptr(iham) = rowptr(iham-1) + rowcount(iham-1)
  END DO
  !
  DEALLOCATE(rowcount)
  !
  CLOSE(fi)
#if defined(__MPI)
  END IF
  !
  call MPI_BCAST(ndim,     1,      MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ndiag,    1,      MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nham,     1,      MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  !
  IF(myrank /= 0) THEN
     ALLOCATE(ham(nham), ham_indx(2,nham), rowptr(ndim+1), & 
            & dvalues(ndim), values(nham-ndiag), colind(nham-ndiag))
  END IF
  !
  call MPI_BCAST(ham,      nham,       MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ham_indx, 2*nham,     MPI_INTEGER,        0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(rowptr,   ndim+1,     MPI_INTEGER,        0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(dvalues,  ndim,       MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(values,   nham-ndiag, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(colind,   nham-ndiag, MPI_INTEGER,        0, MPI_COMM_WORLD, ierr)
  !WRITE(*,*) myrank, "ham", ham(1:10)
  !WRITE(*,*) myrank, "rowptr", rowptr(1:10)
  !WRITE(*,*) myrank, "dvalues", dvalues(1:10)
  !WRITE(*,*) myrank, "values", values(1:10)
  !WRITE(*,*) myrank, "colind", colind(1:10)
  !STOP
  !$ t02 = omp_get_wtime()
  t2 = MPI_WTIME()
  OPEN(10, file='timer_inham.dat', status='replace')
  IF(myrank == 0) THEN
     WRITE(10,*) t2-t1
     !$ WRITE(10,*) "timer:",t02-t01
  END IF
  CLOSE(10)
#else
  !$ t02 = omp_get_wtime()
  open(10, file='timer_inham.dat', status='replace')
  !$ WRITE(10,*) "timer:",t02-t01
  close(10)
#endif
  !
  ! Hermitian(BiCG) or Real-Symmetric(COCG)
  !
!  IF(MAXVAL(ABS(AIMAG(ham(1:nham)))) > almost0) THEN
!     WRITE(stdout,*) "  shifted_qmr_sym mathod is used."
!     solver = "shifted_qmr_sym"
!  ELSE
!     WRITE(stdout,*) "  COCG mathod is used."
!     solver = "COCG"
!  END IF
  !
END SUBROUTINE input_hamiltonian
!
! Input Hamiltonian with the CRS FORMAT
!
subroutine input_hamiltonian_crs()
  !$ use omp_lib
  use shiftk_vals,only : ndim,inham,solver,almost0,stdout,solver,inpunit,myrank
  use ham_vals,only : row_ptr,col_ind,ham_crs_val,ham_crs_val_r,row_se
  implicit none
  integer :: fi=10
  integer :: i,j,n,idim,jdim,nham,ie,ns,ne,numd,num,numham,numlim,nthreads
  integer,allocatable :: irow(:),jcol(:)
  complex(kind(0.0D0)),allocatable :: vval(:)
  double precision,allocatable :: vval_r(:)
  double precision :: rval,cval,t01,t02
  logical :: lbal
  !
  write(stdout,*)
  write(stdout,*) "##########  Input Hamiltonian  ##########"
  write(stdout,*)
  !
  !$ t01 = omp_get_wtime()
  !
  open(fi,file=trim(inham))
  read(fi,*)
  read(fi,*) idim,jdim,nham
  !
  write(stdout,*) "          Dim. of Hamiltonian : ",idim,jdim
  write(stdout,*) "  Num. of Non-Zero Components : ",nham
  !
  if(idim.eq.jdim) then
     ndim=idim
  else
     write(stdout,*) "ERROR ! Hamiltonian is not square."
     stop
  end if
  !
  allocate(irow(nham*2),jcol(nham*2))
  if(trim(solver) == "bicg") then
     allocate(vval(nham*2))
  else
     allocate(vval_r(nham*2))
  endif
  numd=0
  do n=1,nham
     numd=numd+1
     read(fi,*) irow(numd),jcol(numd),rval,cval
     if(trim(solver) == "bicg") then
        vval(numd)=cmplx(rval,cval,kind(0.0D0))
     else
        vval_r(numd)=rval
     endif
     if(irow(numd).eq.jcol(numd))then
     else
        numd=numd+1
        irow(numd)=jcol(numd-1)
        jcol(numd)=irow(numd-1)
        if(trim(solver) == "bicg") then
           vval(numd)=conjg(vval(numd-1))
        else
           vval_r(numd)=vval_r(numd-1)
        endif
     end if
  end do
  close(fi)
  if(trim(solver) == "bicg") then
     call quicksort3(irow,jcol,vval,1,numd)
  else
     call quicksort3_r(irow,jcol,vval_r,1,numd)
  endif
  !
  i=irow(1)
  ns=1
  ne=0
  do n=1,numd
     if(irow(n).eq.i.and.n.lt.numd) then
        ne=ne+1
     else
        if(n.eq.numd) ne=ne+1
        if(trim(solver) == "bicg") then
           call quicksort2(jcol,vval,ns,ne)
        else
           call quicksort2_r(jcol,vval_r,ns,ne)
        endif
        do j=ns+1,ne
           if(jcol(j).eq.jcol(j-1)) then
              write(stdout,*) "ERROR! A duplicate matrix element."
              stop
           end if
        end do
        i=irow(n)
        ne=ne+1
        ns=ne
     end if
  end do
  !
  allocate(row_ptr(idim+1),col_ind(numd))
  if(trim(solver) == "bicg") then
     allocate(ham_crs_val(numd))
  else
     allocate(ham_crs_val_r(numd))
  endif
  row_ptr(:)=-1
  do n=1,numd
     col_ind(n)=jcol(n)
     if(row_ptr(irow(n)).lt.0) then
        row_ptr(irow(n))=n
     end if
     if(trim(solver) == "bicg") then
        ham_crs_val(n)=vval(n)
     else
        ham_crs_val_r(n)=vval_r(n)
     endif
  end do
  row_ptr(idim+1)=n
  deallocate(irow,jcol)
  if(trim(solver) == "bicg") then
     deallocate(vval)
  else
     deallocate(vval_r)
  endif
  !
  nthreads = 1
  !$OMP PARALLEL default(shared)
  !$OMP MASTER
  !$ nthreads = omp_get_num_threads()
  !$OMP END MASTER
  !$OMP END PARALLEL
  allocate(row_se(2,0:nthreads-1))
  lbal=.TRUE.
  if(lbal) then
     numham = row_ptr(ndim+1) - row_ptr(1)
     numlim = numham/nthreads + 1
     ne = 0
     do i = 0, nthreads-1
        row_se(1,i) =  0
        row_se(2,i) = -1
        if(ne.lt.ndim) then
           num = 0
           ns = ne + 1
           row_se(1,i) = ns
           do j = ns, ndim
              num = num + row_ptr(j+1) - row_ptr(j)
              if(num.lt.numlim) then
              else
                 exit
              end if
           end do
           if(j.lt.ndim) then
              ne = j
           else
              ne = ndim
           end if
           row_se(2,i) = ne
        end if
     end do
  else
     numlim = ndim/nthreads
     do i = 0, nthreads-1
        row_se(1,i) = i*numlim+1
        row_se(2,i) = (i+1)*numlim
     end do
     row_se(2,nthreads-1) = ndim
  end if
  !
  !$ t02=omp_get_wtime()
  !
  open(fi, file='timer_inham.dat', status='replace')
  !$ write(fi,*) "timer:",t02-t01
  close(fi)
  !
  return
  !
  ! Hermitian(BiCG) or Real-Symmetric(COCG)
  !
!  IF(MAXVAL(ABS(AIMAG(ham_crs_val(1:numd)))) > almost0) THEN
!     WRITE(stdout,*) "  shifted_qmr_sym mathod is used."
!     solver = "shifted_qmr_sym"
!  ELSE
!     WRITE(stdout,*) "  COCG mathod is used."
!     solver = "COCG"
!  END IF
100 write(*,*) "Stop in INPUT_PARAMETER for CG. reading namelist CG method"
  !
#if defined(__MPI)
  CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  STOP
  !
end subroutine input_hamiltonian_crs
!!
RECURSIVE SUBROUTINE QSORT(a, b, c, left, right)
  IMPLICIT NONE
  COMPLEX(8) :: c(*), ctmp
  INTEGER :: a(*), b(*), pivot, tmp, left, right
  INTEGER :: i, j
  IF(left < right) THEN
     i = left
     j = right
     pivot = a((i+j)/2)
     DO
        DO WHILE(a(i) < pivot)
           i = i+1
        END DO
        DO WHILE(pivot < a(j))
           j = j-1
        END DO
        IF(i >= j) EXIT
        tmp = a(i); a(i) = a(j); a(j) = tmp
        tmp = b(i); b(i) = b(j); b(j) = tmp
        ctmp = c(i); c(i) = c(j); c(j) = ctmp
        i = i+1
        j = j-1
     END DO
     CALL QSORT(a, b, c, left, i-1)
     CALL QSORT(a, b, c, j+1, right)
  END IF
END SUBROUTINE QSORT
!
! Input Right Hand Side Vector
!
SUBROUTINE input_rhs_vector()
  !
#if defined(__MPI)
  USE mpi, only : MPI_COMM_WORLD
#endif
  USE shiftk_vals, ONLY : ndim, rhs, invec, e_min, e_max, nproc, stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ndim2, idim
#if defined(__MPI)
  INTEGER ierr
#endif
  REAL(8) :: rhs_r, rhs_i
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input Right Hand Side Vector ##########"
  WRITE(stdout,*)
  !
  !IF(nproc /= 1) THEN
  !   WRITE(*,*) "ERROR ! MPI is not available in this mode."
#if defined(__MPI)
  !   CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  !   STOP
  !END IF
  !
  OPEN(fi, file = TRIM(invec))
  READ(fi,*) ndim2
  WRITE(stdout,*) "  Dim. of RHS vector : ", ndim2
  !
  IF(ndim2 /= ndim) THEN
     WRITE(stdout,*) "ERROR ! Dimension is Incorrect."
     STOP
  END IF
  !
  ALLOCATE(rhs(ndim))
  !
  DO idim = 1, ndim
     READ(fi,*) rhs_r, rhs_i
     rhs(idim) = CMPLX(rhs_r, rhs_i, KIND(0d0))
  END DO
  !
  CLOSE(fi)
  !
  e_min = 0d0
  e_max = 1d0
  !
END SUBROUTINE input_rhs_vector
!
! Input restart variables from file
!
SUBROUTINE input_restart_parameter()
  !
#if defined(__MPI)
  USE mpi, only : MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_COMPLEX
#endif
  USE shiftk_vals, ONLY : iter_old, alpha, beta, z_seed, nl, r_l_save, myrank, stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, iter, il
#if defined(__MPI)
  INTEGER ierr
#endif
  REAL(8) :: z_seed_r, z_seed_i, alpha_r, alpha_i, beta_r, beta_i
  REAL(8) :: r_l_save_r, r_l_save_i
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input Restart Parameter  ##########"
  WRITE(stdout,*)
  !
  IF(myrank == 0) THEN
     !
     OPEN(fi, file = 'output/TriDiagComp.dat')
     !
     READ(fi,*) iter_old
     WRITE(stdout,*) "  Num. of Iteration (Previous Run) : ", iter_old
     !
     ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl, iter_old))
     !
     READ(fi,*) z_seed_r, z_seed_i
     WRITE(stdout,*) "  Previous Omega_Seed : ", z_seed_r, z_seed_i
     z_seed = CMPLX(z_seed_r, z_seed_i, KIND(0d0))
     !
     ! alpha & beta for CG
     !
     DO iter = 1, iter_old
        !
        READ(fi,*) alpha_r, alpha_i, beta_r, beta_i
        alpha(iter) = CMPLX(alpha_r, alpha_i, KIND(0d0))
        beta(iter) = CMPLX(beta_r, beta_i, KIND(0d0))
        !
     END DO
     !
     ! Projected residual vector
     !
     DO iter = 1, iter_old
        !
        DO il = 1, nl
           READ(fi,*) r_l_save_r, r_l_save_i
           r_l_save(il,iter) = CMPLX(r_l_save_r, r_l_save_i, KIND(0d0))
        END DO
        !
     END DO
     !
     CLOSE(fi)
     !
  END IF ! (myrank == 0)
  !
#if defined(__MPI)
  call MPI_BCAST(iter_old, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  IF(myrank /= 0) THEN
     ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl, iter_old))
  END IF
  call MPI_BCAST(z_seed,   1,             MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(alpha,    iter_old,      MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(beta,     iter_old,      MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(r_l_save, iter_old * nl, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
#endif
  !
END SUBROUTINE input_restart_parameter
!
! Input restart variables from file
!
SUBROUTINE input_restart_vector()
  !
#if defined(__MPI)
  USE mpi, only : MPI_COMM_WORLD
#endif
  USE shiftk_vals, ONLY : ndim, v2, v12, v4, v14, solver, myrank, stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ndim2, idim
  REAL(8) :: v2_r, v2_i, v12_r, v12_i
  CHARACTER(256) :: fname, cmyrank
#if defined(__MPI)
  INTEGER ierr
#endif
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input Restart Vector  ##########"
  WRITE(stdout,*)
  !
  WRITE(cmyrank,*) myrank
  WRITE(fname,'(a,a)') 'output/ResVec.dat', TRIM(ADJUSTL(cmyrank))
  OPEN(fi, file = TRIM(fname))
  !
  READ(fi,*) ndim2
  WRITE(stdout,*) "  Dim. of Residual vector : ", ndim2
  !
  IF(ndim2 /= ndim) THEN
     WRITE(stdout,*) "ERROR ! Dimension is Incorrect."
#if defined(__MPI)
     CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
     STOP
  END IF
  !
  ! Last two residual vectors
  !
  DO idim = 1, ndim
     READ(fi,*) v2_r, v2_i, v12_r, v12_i
     v2(idim) = CMPLX(v2_r, v2_i, KIND(0d0))
     v12(idim) = CMPLX(v12_r, v12_i, KIND(0d0))
  END DO
  !
  ! Last two Shadow residual vectors (Only for BiCG)
  !
  IF(TRIM(solver) == "bicg") THEN
     DO idim = 1, ndim
        READ(fi,*) v2_r, v2_i, v12_r, v12_i
        v4(idim) = CMPLX(v2_r, v2_i, KIND(0d0))
        v14(idim) = CMPLX(v12_r, v12_i, KIND(0d0))
     END DO
  END IF
  !
  CLOSE(fi)
  !
END SUBROUTINE input_restart_vector
!
! Input restart variables from file
!
SUBROUTINE output_restart_parameter()
  !
  USE shiftk_vals, ONLY : iter_old, alpha, beta, z_seed, r_l_save, nl, myrank, stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, iter, il
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Output Restart Parameter  ##########"
  WRITE(stdout,*)
  !
  IF(myrank == 0) THEN
     !
     OPEN(fo, file = 'output/TriDiagComp.dat')
     !
     WRITE(fo,*) iter_old
     WRITE(stdout,*) "  Num. of Iteration (Current Run) : ", iter_old
     !
     WRITE(stdout,'(a,2e13.5)') "   Current Omega_Seed : ", z_seed
     WRITE(fo,'(2e25.16)') z_seed
     !
     DO iter = 1, iter_old
        !
        WRITE(fo,'(4e25.16)') alpha(iter), beta(iter)
        !
     END DO
     !
     DO iter = 1, iter_old
        !
        DO il = 1, nl
           WRITE(fo,'(2e25.16)') r_l_save(il,iter)
        END DO
        !
     END DO
     !
     CLOSE(fo)
     !
  END IF !(myrank == 0)
  !
END SUBROUTINE output_restart_parameter
!
! Input restart variables from file
!
SUBROUTINE output_restart_vector()
  !
  USE shiftk_vals, ONLY : ndim, v2, v12, v4, v14, solver, myrank, stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, idim
  CHARACTER(256) :: fname, cmyrank
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Output Restart Vector  ##########"
  WRITE(stdout,*)
  !
  WRITE(cmyrank,*) myrank
  WRITE(fname,'(a,a)') 'output/ResVec.dat', TRIM(ADJUSTL(cmyrank))
  OPEN(fo, file = TRIM(fname))
  !
  WRITE(fo,*) ndim
  WRITE(stdout,*) "  Dim. of Residual vector : ", ndim
  !
  DO idim = 1, ndim
     WRITE(fo,'(4e25.16)') v2(idim), v12(idim)
  END DO
  !
  IF(TRIM(solver) == "bicg") THEN
     DO idim = 1, ndim
        WRITE(fo,'(4e25.16)') v4(idim), v14(idim)
     END DO
  END IF
  !
  CLOSE(fo)
  !
END SUBROUTINE output_restart_vector
!
! Output Resulting Spectrum
!
SUBROUTINE output_result()
  !
  USE shiftk_vals, ONLY : x_l, z, nomega, myrank
  !
  IMPLICIT NONE
  !
  INTEGER :: iz, fo = 20
  !
  IF(myrank == 0) THEN
     !
     OPEN(fo, file = "output/dynamicalG.dat")
     !
     DO iz = 1, nomega
        !
        write(fo,'(4e13.5)') DBLE(z(iz)), AIMAG(z(iz)), DBLE(x_l(1,iz)), AIMAG(x_l(1,iz))
        !
     END DO
     !
     CLOSE(fo)
     !
  END IF ! (myrank == 0)
  !
END SUBROUTINE output_result
!
! Output Resulting Spectrum, True residual
!
SUBROUTINE output_result_debug()
  !
  USE shiftk_vals, ONLY : v2, ndim, x_l, rhs, z, nomega, myrank, stdout
  USE ham_prod_mod, ONLY : ham_prod
  USE lobpcg_mod, ONLY : zdotcMPI
  !
  IMPLICIT NONE
  !
  INTEGER :: iz, fo = 20
  CHARACTER(100) :: cndim, form
  COMPLEX(8) :: Gii
  REAL(8) :: res
  double precision :: t11,t12
  logical :: lcollect=.FALSE.
  !
  WRITE(cndim,*) ndim * 2
  WRITE(form,'(a,a,a)') "(", TRIM(ADJUSTL(cndim)), "e13.5)"
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "#####  Check Results  #####"
  WRITE(stdout,*)
  !
  WRITE(stdout,*) "  Residual Vector"
  !
  DO iz = 1, nomega
     !
     CALL ham_prod(x_l(1:ndim,iz), v2, t11, t12, lcollect)
     v2(1:ndim) = z(iz) * x_l(1:ndim,iz) - v2(1:ndim) - rhs(1:ndim)
     !
     res = SQRT(DBLE(zdotcMPI(ndim, v2(1:ndim), v2(1:ndim))))
     !write(*,form) v2(1:ndim)
     IF (myrank == 0) write(*,'(a,i5,a,2e13.5,a,e13.5)') &
     &        "DEBUG (", iz, "), omega = ", z(iz), ", Res. = ", res
     !
  END DO
  !
  IF (myrank == 0) OPEN(fo, file = "output/dynamicalG.dat")
  !
  DO iz = 1, nomega
     !
     Gii = zdotcMPI(ndim, rhs,x_l(1:ndim,iz))
     IF (myrank == 0) write(fo,'(4e13.5)') DBLE(z(iz)), AIMAG(z(iz)), DBLE(Gii), AIMAG(Gii)
     !
  END DO
  !
  IF (myrank == 0) CLOSE(fo)
  !
END SUBROUTINE output_result_debug
!
recursive subroutine quicksort3(ii,jj,cc,is,ie)
  implicit none
  integer,intent(IN) :: is,ie
  integer,intent(INOUT) :: ii(:),jj(:)
  complex(kind(0.0D0)),intent(INOUT) :: cc(:)
  !
  integer :: i,j,h,ti
  complex(kind(0.0D0)) :: tc
  !
  h=ii((is+ie)/2)
  i=is
  j=ie
  do
     do while (ii(i).lt.h)
        i=i+1
     end do
     do while (h.lt.ii(j))
        j=j-1
     end do
     if(i.ge.j) exit
     ti=ii(i)
     ii(i)=ii(j)
     ii(j)=ti
     ti=jj(i)
     jj(i)=jj(j)
     jj(j)=ti
     tc=cc(i)
     cc(i)=cc(j)
     cc(j)=tc
     i=i+1
     j=j-1
  end do
  !
  if(is.lt.i-1) call quicksort3(ii,jj,cc,is,i-1)
  if(j+1.lt.ie) call quicksort3(ii,jj,cc,j+1,ie)
  !
end subroutine quicksort3
!
recursive subroutine quicksort3_r(ii,jj,rr,is,ie)
  implicit none
  integer,intent(IN) :: is,ie
  integer,intent(INOUT) :: ii(:),jj(:)
  double precision,intent(INOUT) :: rr(:)
  !
  integer :: i,j,h,ti
  double precision :: tr
  !
  h=ii((is+ie)/2)
  i=is
  j=ie
  do
     do while (ii(i).lt.h)
        i=i+1
     end do
     do while (h.lt.ii(j))
        j=j-1
     end do
     if(i.ge.j) exit
     ti=ii(i)
     ii(i)=ii(j)
     ii(j)=ti
     ti=jj(i)
     jj(i)=jj(j)
     jj(j)=ti
     tr=rr(i)
     rr(i)=rr(j)
     rr(j)=tr
     i=i+1
     j=j-1
  end do
  !
  if(is.lt.i-1) call quicksort3_r(ii,jj,rr,is,i-1)
  if(j+1.lt.ie) call quicksort3_r(ii,jj,rr,j+1,ie)
  !
end subroutine quicksort3_r
!
recursive subroutine quicksort2(ii,cc,is,ie)
  implicit none
  integer,intent(IN) :: is,ie
  integer,intent(INOUT) :: ii(:)
  complex(kind(0.0D0)),intent(INOUT) :: cc(:)
  !
  integer :: i,j,h,ti
  complex(kind(0.0D0)) :: tc
  !
  h=ii((is+ie)/2)
  i=is
  j=ie
  do
     do while (ii(i).lt.h)
        i=i+1
     end do
     do while (h.lt.ii(j))
        j=j-1
     end do
     if(i.ge.j) exit
     ti=ii(i)
     ii(i)=ii(j)
     ii(j)=ti
     tc=cc(i)
     cc(i)=cc(j)
     cc(j)=tc
     i=i+1
     j=j-1
  end do
  !
  if(is.lt.i-1) call quicksort2(ii,cc,is,i-1)
  if(j+1.lt.ie) call quicksort2(ii,cc,j+1,ie)
  !
end subroutine quicksort2
!
recursive subroutine quicksort2_r(ii,rr,is,ie)
  implicit none
  integer,intent(IN) :: is,ie
  integer,intent(INOUT) :: ii(:)
  double precision,intent(INOUT) :: rr(:)
  !
  integer :: i,j,h,ti
  double precision :: tr
  !
  h=ii((is+ie)/2)
  i=is
  j=ie
  do
     do while (ii(i).lt.h)
        i=i+1
     end do
     do while (h.lt.ii(j))
        j=j-1
     end do
     if(i.ge.j) exit
     ti=ii(i)
     ii(i)=ii(j)
     ii(j)=ti
     tr=rr(i)
     rr(i)=rr(j)
     rr(j)=tr
     i=i+1
     j=j-1
  end do
  !
  if(is.lt.i-1) call quicksort2_r(ii,rr,is,i-1)
  if(j+1.lt.ie) call quicksort2_r(ii,rr,j+1,ie)
  !
end subroutine quicksort2_r
!
END MODULE shiftk_io
