program lapack_solve
      !Purpose: To solve a set of linear equations: -H.d = qE
      !Format is Ax = B where A = -H, x = d and B = qE
      use iso_c_binding
      implicit none

      !Generating hessian and born matrices for extraction from .npy files
      integer, parameter :: DPR = selected_real_kind(p=15)
      character(len=*), parameter :: hess_data = 'hess_matrix', &
              born_data = 'born_matrix', UPLO = 'L'
      
      !Important variables : matrix dimension = 3*No. atoms in unit cell
      !i,j for do loops
      integer :: N, i,j
      real(DPR), allocatable, dimension(:,:) :: h_tensor, b_tensor
      
      !Defining variables required for calling dysv
      
      !Integers
      integer, parameter :: nout = 12
      integer :: INFO, NRHS, LDA, LDB, LWORK
      
      !Arrays
      integer, allocatable, dimension(:) :: IPIV
      double precision, allocatable, dimension(:,:) :: A, q
      double precision, allocatable, dimension(:) :: E, B, WORK
      !Defining electric field parameters to be able to generate B
      real, parameter :: Ex = 3.0, Ey = 2.0, Ez = 1.0
      
      print *, 'Dimension of Hessian and Born tensor: '
      read *, N
      NRHS = 1
      LDA = N
      LDB = N
      LWORK = 3*N
      allocate (h_tensor(N,N), b_tensor(N,N), A(LDA,N), q(N,N), E(N), B(N), WORK(LWORK), IPIV(N))
      write (nout,*) 'DSYSV -Hd = qE Results'

      !Writing Hessian and Born tensors from .npy files from disp_solve.py
      open(10, file=hess_data, status = 'old', access = 'stream',&
              form = 'unformatted')
      read(10) h_tensor
      close(10)
      open(11, file=born_data, status = 'old', access = 'stream',&
              form = 'unformatted')
      read(11) b_tensor
      close(11)
      
      !print *, "Hessian: "
      !print *, h_tensor
      !print *, "Born tensor: "
      !print *, b_tensor
      
      !Generating the A and q matrices
      do i = 1,N
        do j = 1,N
                A(i,j) = dble(-1.0*h_tensor(i,j))
                q(i,j) = dble(b_tensor(j,i))
        end do
      end do

      !print *, A
      !Generating electric field matrix
      do i = 1,N
        if(mod(i,3) == 1) then
                E(i) = Ex
        else if(mod(i,3) == 2) then
                E(i) = Ey
        else
                E(i) = Ez
        end if
      end do

            !print *, q
      !print *, "######################"

      !Generating B matrix
      B = matmul(q,E)
      
      write(*,*) B
      
      !do i = 1,2
      !  do j = 1,N
      !          print *, A(i,j)
      !          if(j==N) then
      !                  print *, "##############################"
      !          end if
      !  end do
      !end do
      
      call dsysv(UPLO,N,NRHS,A,LDA,IPIV,B,LDB,WORK,LWORK,INFO)

      if(info==0) then
              !Print the solution
              
              write (12,*) 'Solutions'
              write (nout,100) B(1:N)

              !Print factorization details

              write (nout, *)
              flush (nout)

              !Print pivot indices
              write (nout,*)
              write (nout, *) 'Pivot indices'
              write (nout,110) ipiv(1:n)
      else
              write (nout, 120) 'The diagonal block ', info, ' of D is zero'
      end if
      100 format (7(e18.10,2X))
      110 format (1x, 7i11)
      120 format (1x, A, I3, A)

! call as:
!gfortran lapack_solve.f03 -L/usr/lib -llapack -L/usr/lib -lblas
end program lapack_solve
