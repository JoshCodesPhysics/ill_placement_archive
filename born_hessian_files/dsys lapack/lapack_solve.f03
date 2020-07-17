program lapack_solve
      !Purpose: To solve a set of linear equations: -H.d = qE
      !Format is Ax = B where A = -H, x = d and B = qE
      use iso_c_binding
      implicit none

      !Generating hessian and born matrices for extraction from .npy files
      integer, parameter :: DPR = selected_real_kind(p=15)
      character(len=*), parameter :: hess_data = 'hess_matrix', &
              born_data = 'born_matrix'
      
      !Important variable : matrix dimension = 3*No. atoms in unit cell
      integer, parameter :: N = 192
      real(DPR), dimension(N,N) :: h_tensor, b_tensor
      
      !Defining variables required for calling dysv
      
      !Integers
      integer, parameter :: NRHS = 1, LDA = N, LDB = N
      !Not sure what to assign these
      integer :: LWORK, INFO
      
      !Arrays
      integer, dimension(N) :: IPIV
      real, dimension(N) :: E
      double precision, dimension(N,N) :: A, B, q
      !Allocatable because I don't know how to define LWORK which
      !determines the length of WORK
      double precision, allocatable, dimension(:) :: WORK
      !Defining electric field parameters to be able to generate B
      real, parameter :: Ex = 3.0, Ey = 2.0, Ez = 1.0
      integer :: i,j
      
      !I'm not sure if I even need these to be honest:
      
      !Local scalars
      !logical :: lquery
      !integer :: lwkopt
      !External functions
      !logical :: LSAME
      !external :: lsame
      !External subroutines
      !external :: xerbla, dsytrf, dsytrs, dsytrs2
      !Intrinsic functions
      !intrinsic :: max
      
      !Writing Hessian and Born tensors from .npy files from disp_solve.py
      open(10, file=hess_data, status = 'old', access = 'stream',&
              form = 'unformatted')
      read(10) h_tensor
      close(10)
      open(11, file=born_data, status = 'old', access = 'stream',&
              form = 'unformatted')
      read(11) b_tensor
      close(11)
      
      !Generating the A matrix
      do i = 1,N
        do j = 1,N
                A(i,j) = dble(-1.0*h_tensor(i,j))
                q(i,j) = dble(b_tensor(i,j))
        end do
      end do

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

      !These have incompatible ranks for some reason
      print *, "Rank of born tensor: ", rank(q), " Rank of E vector: ", rank(E), " Cannot multiply :("
      !Trying to generate B matrix unsuccessfully:
      !B = matmul(q,E)
       
      !call dsysv(UPLO,N,NRHS,A,LDA,B,LDB,LWORK)
end program lapack_solve
