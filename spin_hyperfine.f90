!*********************************************
!To compile program:
!
! gfortran spin_hyperfine_module.F90 -lblas -fopenmp -llapack -o spin_hyperfine.exe
!
! Program can also be compiled without -lblas, -llapack and/or -fopenmp for debugging
!
!Tested with LAPACK-3.4.2
!*********************************************


module decoherence_stuff
  implicit none
  
  !Simulation parameters
  integer,          parameter :: N = 1 				!Number of Nuclei in simulation
  integer,          parameter :: t_max = 10000 			!Maximum number of timesteps
  integer,          parameter :: dp=selected_real_kind(15,300)		!precision of numerical values
  real(kind=dp),    parameter :: interaction_constant = -0.04_dp 	!Interaction constant for hyperfine interaction
  complex(kind=dp), parameter :: TDSE_constant = (0.0_dp, 1.0E-5_dp)	!TDSE wavefunction constant
  integer         , parameter :: QN = 2**(N+1)				!'Quantum' N, defined as 2^(N+1) which is the number of possible states
  
  
  !Simulation Variables
  complex(kind=dp), save, dimension(1:QN)      :: C            !Complex numbers corresponding to basis states
  real(kind=dp),    save, dimension(1:QN,1:QN) :: H            !Hamiltonian matrix for system
  real(kind=dp),    save, dimension(1:QN,1:QN) :: basis        !matrix of basis vectors (diagonalised)
  real(kind=dp),    save, dimension(1:QN)      :: e_spin_up_wf !Wavefunction corresponding to superpositon of electron spin up states
  
  !Error reporting
  integer,          save :: ierr !Error integer
  
  
contains

  !*********************************************
  !
  !Performs Kronecker product of two matrices A and B producing matrix P
  !n and m correspond to the dimensions of A and x and y correspond to the dimensions of B
  !dimensions of P are assumed to be n*x and m*y
  !
  !*********************************************
  subroutine kronecker_product(A, B, P, n, m, x, y)
    implicit none
    integer :: i, j, k, l
    integer :: n, m, x, y
    complex(kind=dp), dimension(n, m) :: A
    complex(kind=dp), dimension(x, y) :: B
    complex(kind=dp), dimension(x*n, y*m) :: P
    
    do i = 1, n
      do j = 1, m
      
	do k = 1, x
	  do l = 1, y
	    
	    P((i-1)*x + k, (j-1)*y + l) = A(i, j) * B(k, l)

	  end do
	end do

      end do
    end do
    
    return
  
  end subroutine kronecker_product
  
  !*********************************************
  !Generates a hamiltonian for the hyperfine interaction from pauli matrices
  !*********************************************
  subroutine hyperfine_hamiltonian
    implicit none
    integer :: i, j, k
    complex(kind=dp), dimension(3,2,2) :: pauli
    real(kind=dp), dimension(2,2) :: identity
    integer, dimension(1:N,1:N) :: diag_matrix
    complex(kind=dp), dimension(:,:), Allocatable :: A
    complex(kind=dp), dimension(2,2) :: B
    complex(kind=dp), dimension(:,:), Allocatable :: P
    
    pauli(1, 1, 2) = 1.0_dp		!Pauli x
    pauli(1, 1, 1) = 0.0_dp		!Pauli x
    pauli(1, 2, 2) = 0.0_dp		!Pauli x
    pauli(1, 2, 1) = 1.0_dp		!Pauli x
    
    pauli(2, 1, 1) = 1.0_dp		!Pauli z
    pauli(2, 1, 2) = 0.0_dp		!Pauli z
    pauli(2, 2, 1) = 0.0_dp		!Pauli z
    pauli(2, 2, 2) = -1.0_dp		!Pauli z
    
    pauli(3, 1, 2) = (0.0_dp, -1.0_dp) !Pauli y
    pauli(3, 1, 1) = 0.0_dp		!Pauli y
    pauli(3, 2, 2) = 0.0_dp		!Pauli y
    pauli(3, 2, 1) = (0.0_dp, 1.0_dp)  !Pauli y
    
    identity(1, 1) = 1.0_dp		!Identity matrix
    identity(2, 2) = 1.0_dp		!Identity matrix
    identity(1, 2) = 0.0_dp		!Identity matrix
    identity(2, 1) = 0.0_dp		!Identity matrix
    
    !Creates a diagonal matrix for permuations of spin matrices
    do i = 1, N
      diag_matrix(i, i) = 1
    end do
    
    !k signifies each of the cartesian dimensions
    !x = 1
    !z = 2
    !y = 3
    do k = 1, 3
    
      !Loops over number of possible combinations of identity and pauli matrices
      do i = 1, N
	
	!Allocates A input matrix for kronecker_product first iteration
	Allocate(A(1:2, 1:2))
	if (ierr/=0) stop 'Error in allocating matrix initial A'
	
	!Allocate matrix result for kronecker product so it set as first pauli matrix
	!and then can be deallocated in loop
	Allocate(P(1:2, 1:2))
	if (ierr/=0) stop 'Error in allocating matrix initial P'
	
	!Sets up the first A input matrix as the k dimension pauli matrix
	P = pauli(k, 1:2, 1:2)
	
	!Loops over each matrix to be multiplied for this iteration
	!Each j is a possible interaction pair
	do j = 1, N
	  
	  !Sets previous product as next input
	  A = P
	  
	  !Deallocates product matrix so it can be reshaped
	  Deallocate(P)
	  if (ierr/=0) stop 'Error in deallocating matrix P in loop'
	  
	  !Finds the next matrix product based on diagonlised matrix
	  if(diag_matrix(i, j) == 1) then
	    !Sets the B input matrix as the kth dimensional pauli matrix
	    B(1:2, 1:2) = pauli(k, 1:2, 1:2)
	  else
	    !Sets the B input matrix as the identity matrix (non-interaction)
	    B(1:2, 1:2) = identity(1:2, 1:2)
	  end if
	  
	  !Allocates P matrix to the size of size of the next k_product matrix
	  Allocate(P(1:2**(j+1), 1:2**(j+1)))
	  
	  !Calls the tensor product for A and B producing P
	  !Dimensions of A are determined by place in j loop
	  call kronecker_product(A, B, P, 2**j, 2**j, 2, 2)
	  
	  !Deallocates input matrix A to prepare for next loop
	  Deallocate(A)
	  if (ierr/=0) stop 'Error in deallocating matrix A in loop'
	  
	  !Allocates input matrix A in next size to prepare for next loop
	  Allocate(A(1:2**(j+1), 1:2**(j+1)))
	  if (ierr/=0) stop 'Error in allocating matrix A in loop'
	  
	end do
	
	!Deallocates A for this interaction
	Deallocate(A)
	if (ierr/=0) stop 'Error in deallocating matrix A at end of loop'
	
	!Adds P to hamiltonian in matrix at each loop
	!Allows for additional hamiltonian elements to be added later
	H(1:QN,1:QN) = H(1:QN,1:QN) + P(1:QN, 1:QN)
	
	Deallocate(P)
	if (ierr/=0) stop 'Error in deallocating matrix P at end'
	
      end do
    end do

    !Writes Hamiltonian to file for debugging
    open(400, file="hamiltonian.dat")
    do i = 1, QN
      write(400, *) (INT(H(i, j)), j = 1, QN)
    end do
    close(400)
    
    !Multiplies all hamiltonian elements by hyperfine interaction constant
    !H = interaction_constant * H
    
    return
  
  end subroutine hyperfine_hamiltonian
  
  !*********************************************
  !
  !*********************************************
  subroutine magnetic_hamiltonian
  !placeholder
  end subroutine magnetic_hamiltonian
  
  !*********************************************
  ! Uses LAPACK to find eigenstates of hamiltonian H using DGEEV functions
  ! DSYEVD finds the eigenvectors and eigenvalues of symmetric matrices
  !
  ! LAPACK reference:
  !@BOOK{laug,
  !    AUTHOR = {Anderson, E. and Bai, Z. and Bischof, C. and
  !              Blackford, S. and Demmel, J. and Dongarra, J. and
  !              Du Croz, J. and Greenbaum, A. and Hammarling, S. and
  !              McKenney, A. and Sorensen, D.},
  !    TITLE = {{LAPACK} Users' Guide},
  !    EDITION = {Third},
  !    PUBLISHER = {Society for Industrial and Applied Mathematics},
  !    YEAR = {1999},
  !    ADDRESS = {Philadelphia, PA},
  !    ISBN = {0-89871-447-8 (paperback)} }
  !*********************************************
  subroutine find_eigenstates
    implicit none
    
    integer :: i, j							!Loop integers
    integer, parameter :: LIWORK = 64

    
    real(kind=dp), dimension(1:QN) :: WR			!Work (real)
    real(kind=dp), dimension(1:QN) :: WI			!Work (imaginary)
    real(kind=dp), dimension(1:QN, 1:QN) :: A
    real(kind=dp), dimension(1:QN,1:QN) :: VL
    real(kind=dp), dimension(1:QN,1:QN) :: VR
    integer, dimension(1:LIWORK) :: WORK
    
    
    !Because A is overwritten by function
    A = H
    

    call DGEEV( 'N', 'V', QN, A, QN, WR, WI, VL, 1, VR, QN, WORK, LIWORK, ierr)


    if (ierr/=0) then
      print *, 'Error in solving eigenstates of hamiltonian'
      if (ierr<0) then
	print *, 'value', ierr, 'had illegal value'
      else
	print *, 'dgeev failed to converge'
      end if
      stop
    end if

    open(600, file='eigenstates.dat', iostat=ierr)
    if (ierr==0) then
      do j=1, QN
	if(WI(j)==0.0_dp) then
	  write(600,*) 'Eigenvalue', j, '=', WR(j)
	else
	  write(600,*) 'Eigenvalue', j, '=', WR(j), '+', WI(i), 'i'
	end if
	
	write(600,*)
	write(600,*) 'Eigenvector', j,  '='
	write(600,*)
	
	if(WI(j)==0.0_dp) then
	  do i=1, QN
	    write(600, *) VR(i, j)
	  end do
	elseif(WI(j)>=0.0_dp) then
	  do i=1, QN
	    write(600,*) VR(i, j), '+', VR(i, j+1), 'i'
	  end do
	else
	  do i=1, QN
	    write(600,*) VR(i, j-1), '+', VR(i, j), 'i'
	  end do
	 end if
	 
	write(600,*)
	
      end do
    else
      stop 'Error in opening file eigenstates.dat'
    end if
    
    close(600, iostat=ierr)
    if (ierr/=0) stop 'Error in closing file eigenstates.dat'
    
  end subroutine find_eigenstates
  
  !*********************************************
  !Generates basis functions
  !
  !
  !
  !
  !*********************************************
  subroutine init_basis
    implicit none
    integer :: i, j          !loop variables
    
    !Sets all values in basis matrix to 0 with double precision  
    basis(:, :) = 0.0_dp
    
    !Sets up basis functions
    do i = 1, QN
      basis(i,i) = 1.0_dp
    end do
    
    !sets up complex values and 
    e_spin_up_wf(:) = 0.0_dp
    C(:) = 0.0_dp
    
    do i = 1, (2**N)
      e_spin_up_wf(i) = 1.0_dp / (SQRT(REAL(2**N, kind=dp)))
      C(i) = 1.0_dp / (SQRT(REAL(2**N, kind=dp)))
    end do
    
    
  end subroutine init_basis
  
  !*********************************************
  !Finds "fidelity" of system in current state
  !
  !
  !
  !
  !
  !*********************************************
  subroutine fidelity(Fid)
    implicit none
    
    
    !Subroutine variables
    real(kind=dp) :: Fid           !Fidelity variable for subroutine
    real(kind=dp) :: DDOT          !DDOT for blas library
    
!If BLAS is enabled when compiling (-lblas), this code segment
!will be compiled into final program using BLAS functions
#ifdef lblas
    
      Fid = (abs(DDOT(QN, e_spin_up_wf,1, C, 1)))**2
    
!If BLAS is not enabled when compiling, this section of intrinsic functions
!will be used instead.
#else
    
      Fid = (abs(DOT_PRODUCT(e_spin_up_wf, C)))**2
    
!ends the preprocessor if statement
#endif
    
    !returns value of Fid from subroutine to main program
    return
 
  end subroutine fidelity
  
  !*********************************************
  !Integrates time dependant schroedinger equation with simple numeric differnetiation
  !
  !
  !
  !
  !
  !*********************************************
  subroutine intergrate_TDSE
    implicit none
    
    !Library header for OpenMP
    !Required for single memory multiple processor operations
    !Requires compiler flag -fopenmp in gnu compilers to enable
    include 'omp_lib.h'
    
    !subroutine variables
    integer :: i, j                                     !Loop integers
    real(kind=dp) :: DDOT                               !BLAS library double precision dot product (required)

    real(kind=dp) :: Hij                                !Dot product of Hi and transpose of basis(j)
    complex(kind=dp) :: sum_C_H                         !Sum of Complex values multipled by Hij
    real(kind=dp), dimension(1:QN):: Hi           !Vector resulting from matrix vector multiply of H and basis(i)
    
    
    !iterate j from 1 to QN
    do j=1,QN
      
      !Set summation to zero to initial addition in next loop
      sum_C_H = 0.0_dp
      
      !$OMP parallel do private(i, Hij, Hi) shared(C, H, j) reduction(+:sum_C_H)
      !iterate i from 1 to QN
      !With OpenMP enabled this loop is performed on each thread seperately
      !Values of sum_C_H are added together from each thread at end of loop
      do i=1,QN
	
!If BLAS is enabled when compiling (-lblas), this code segment
!will be compiled into final program using BLAS functions
#ifdef lblas
	
	  !call for double precision matrix vector multiply from BLAS
	  !Multiplies Hamiltonian matrix with basis(i) to form vector Hi
	  call DGEMV( 'N', QN, QN, 1.0_dp, H(1:QN, 1:QN), QN, &
	  basis(i,1:QN), 1, 0.0_dp, Hi(1:QN), 1) 
      
	  !call for double precision dot product from BLAS
	  !Performs dot product of transpose of basis(j) with Hi to form value Hij
	  Hij = DDOT(QN,basis(j, 1:QN),1,Hi(1:QN),1)
	
!If BLAS is not enabled when compiling, this section of intrinsic functions
!will be used instead.
#else
	
	  !call for intrinsic matrix vector multiply
	  !Multiplies Hamiltonian matrix with basis(i) to form vector Hi
	  Hi = MATMUL(H, basis(i, 1:QN))
      
	  !call for intrinsic dot product
	  !Performs dot product of transpose of basis(j) with Hi to form value Hij
	  Hij = DOT_PRODUCT(basis(j, 1:QN), Hi)
	
!Ends preprocessor if block
#endif
	

	!Adds product of C(i) and Hij to summation
	!This is processed individually on each thread then summed at end
	!if openMP is enabled (-fopenmp)
	sum_C_H = sum_C_H + C(i) * Hij      

      end do
      !$OMP end parallel do
      !ends the parallel do loop for openMP
    
      !Find C(j) by multiplying TDSE_constant with sum_C_H
      C(j) = C(j) - TDSE_constant * sum_C_H
      
    end do

    
  
  end subroutine intergrate_TDSE

  
end module decoherence_stuff

!*********************************************
!
!Main program
!
!
!
!
!*********************************************
program hyperfine_spin_interaction
  use decoherence_stuff
  implicit none
  integer :: t
  integer :: i, j
  real(kind=dp) :: Fid
  real(kind=dp) :: sum_C_H
  real(kind=dp), dimension(1:QN):: Hi
  real(kind=dp) :: Hij
  
  !calls the hyperfine interaction hamiltonian subroutine
  !adds a hyperfine interaction hamiltonian for N nuclei and 1 electron
  !to H
  call hyperfine_hamiltonian

!If LAPACK is enabled when compiling (-llapack), this code segment
!will be compiled into final program using LAPACK functions 
!#ifdef llapack

  !calls the eigenstate solution subroutine which finds the
  !eigenvectors and eigenvalues of the hamiltonian
  !and writes them to file eigenstates.dat
  call find_eigenstates
  
!#endif
  
  !calls the subroutine to initiate basis states
  call init_basis
 
  
    !Open output files for use in time loop
    !Opens file for c_vs_time.dat and assigns it value 100 and reports error if failure
    open(100, file="c_vs_time.dat",iostat=ierr)
    if (ierr/=0) stop 'Error in opening file c_vs_time.dat'
    
    !Opens file for sum_c_vs_time.dat and assigns it value 200 and reports error if failure
    open(200, file="sum_c_vs_time.dat",iostat=ierr)
    if (ierr/=0) stop 'Error in opening file sum_c_vs_time.dat'
    
    !Opens file for fidelity_vs_time.dat and assigns it value 300 and reports error if failure
    open(300, file="fidelity_vs_time.dat",iostat=ierr)
    if (ierr/=0) stop 'Error in opening file fidelity_vs_time.dat'
  
  !Loops over time from t= 1 to tmax which is set in module header
  do t = 1, t_max
  
    !Calls the subroutine to integrate the TDSE
    call intergrate_TDSE
    
    !Calls the fidelity function which finds the fidelity
    !of the system compared to e_spin_up_wf which
    !is a superpostion of the states which contain a spin up electron
    call fidelity(Fid)


    !write to output file
    !writes probability of each state to file in sequence as well as timesteps
    !writes to c_vs_time.dat
    write(100, *) t, (abs(C(i))**2.0_dp, i = 1, QN)
    
    !writes sum of probabilities and time step to file
    !writes to sum_c_vs_time.dat
    write(200, *) t, sum(abs(C)**2.0_dp)
    
    !writes fidelity of system compared to stated value with timestep
    !writes to fidelity_vs_time.dat
    write(300, *) t, Fid
    
  
  end do
  
  !Closes files
  !Closes file c_vs_time.dat and reports error if failure
  close(100,iostat=ierr)
  if (ierr/=0) stop 'Error in closing file c_vs_time.dat'
  
  !Closes file sum_c_vs_time.dat and reports error if failure
  close(200,iostat=ierr)
  if (ierr/=0) stop 'Error in closing file sum_c_vs_time.dat'
  
  !Closes file fidelity_vs_time.dat and reports error if failure
  close(300,iostat=ierr)
  if (ierr/=0) stop 'Error in opening file fidelity_vs_time.dat'
  
end program hyperfine_spin_interaction
  
  
