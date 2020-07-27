program lapack_solve
      !Purpose: To solve a set of linear equations: -H.d = qE
      !Format is Ax = B where A = -H, x = d and B = qE
      use iso_c_binding
      implicit none

      !Generating hessian and born matrices for extraction from .npy files
      integer, parameter :: DPR = selected_real_kind(p=15), nout = 14
      
      character(len=*), parameter :: hess_data = 'hess_matrix', &
              born_data = 'born_matrix',e_vector = 'e_vector', &
              n_file = 'n_array',sol_file = 'disp.sols.cart', UPLO = 'L'
      
      logical :: exist
      
      !Defining variables required for calling dysv
      
      !Integers
      integer :: N, i,j
      integer :: INFO, NRHS, LDA, LDB, LWORK
      
      !Arrays
      double precision, allocatable, dimension(:,:) :: h_tensor, b_tensor
      !real(DPR), allocatable, dimension(:,:) :: h_tensor, b_tensor
      integer, allocatable, dimension(:) :: IPIV
      integer, dimension(1) :: nlist
      double precision, allocatable, dimension(:,:) :: A, q
      double precision, allocatable, dimension(:) :: E, B, B_COPY,DIFF, WORK
      
      open(10, file = n_file, status='old',access = 'stream',&
              form = 'unformatted')
      read(10) nlist
      close(10, status = 'delete')
      N = nlist(1)
      print "(a27,i5)", "Dimension of the matrices: ", N
      NRHS = 1
      LDA = N
      LDB = N
      LWORK = 3*N
      allocate (h_tensor(N,N), b_tensor(N,N), A(LDA,N), q(N,N), E(N),&
              B(N), B_COPY(N), WORK(LWORK), IPIV(N))

      !Writing Hessian and Born tensors from binary files generated from disp_solve.py
      open(11, file=hess_data, status = 'old', access = 'stream',&
              form = 'unformatted')
      read(11) h_tensor
      close(11,status = 'delete')

      open(12, file=born_data, status = 'old', access = 'stream',&
              form = 'unformatted')
      read(12) b_tensor
      close(12,status = 'delete')
      !Same for electric field vector

      open(13,file=e_vector, status = 'old', access= 'stream',&
              form = 'unformatted')
      read(13) E
      close(13,status = 'delete')
      
      !Generating the A and q matrices
      !Have to transpose q, but not A since it is symmetric
      do i = 1,N
        do j = 1,N
                A(i,j) = dble(-1.0*h_tensor(i,j))
                q(i,j) = dble(b_tensor(j,i))
        end do
      end do

      !Generating B matrix
      B = matmul(q,E)
      B_COPY = matmul(q,E)
      !write(*,*) B
      
      !do i = 1,2
      !  do j = 1,N
      !          print *, A(i,j)
      !          if(j==N) then
      !                  print *, "##############################"
      !          end if
      !  end do
      !end do
      !print *, "before: ", B

      call dsysv(UPLO,N,NRHS,A,LDA,IPIV,B,LDB,WORK,LWORK,INFO)

      !print *, "after: ", B

      !Difference error matrix
      DIFF = abs(matmul(A,B) - B_COPY)
      print *, "Difference matrix:"
      print *, DIFF

      inquire(file=sol_file,exist=exist)
      if(exist) then
              open(nout,file=sol_file,status='old')
      else
              open(nout,file=sol_file,status='new')
      end if
      if(info==0) then
              !Print the solution
              write (nout,*) 'DSYSV -Hd = qE Cartesian Results'
              write (nout,*) 'Solutions'
              write (nout,100) B(1:N)

              !Print factorization details

              write (nout, *)
              flush (nout)

              !Print pivot indices
              write (nout,*)
              write (nout, *) 'Pivot Indices'
              write (nout,110) ipiv(1:n)
              print '(a58)', "Cartesian solution has been written to disp.sols.cart file"
      else
              write (nout, 120) 'The diagonal block ', info, ' of D is zero'
              print *, "See disp.sol.cart file for outcome"
      end if
      close(nout)
      100 format (7(e18.10,2X))
      110 format (1x, 7i11)
      120 format (1x, A, I3, A)

end program lapack_solve
