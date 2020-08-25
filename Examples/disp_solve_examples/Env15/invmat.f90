!!$***********************************************************************
!!$
subroutine det(d,M,n)
!!$
!!$***********************************************************************
!!$
!!$          Routine de calcul du determinant d'une matrice
!!$          (avec les permutations)
!!$
!!$***********************************************************************
!!$
!!$----- Déclarations  
!!$
       implicit none
       real*8, intent(out) :: d
       integer, intent(in) :: n
       real*8, dimension(n,n) :: M
       integer, dimension(n) :: p, saut
       integer :: f, i, j, k, sgn
       real*8 :: x
!!$
!!$--- nombre de permutations :
       f = 1
       do i=1, n
          f = f * i
       end do
!!$--- initialisation
       saut(:) = 0
       p(:) = 0
       sgn = 1
       d = 0.
!!$--- boucle sur les permutations
       do i=1, f
          !--- generation de la permutation
          p(:) = (/ (j,j=1,n) /)
          sgn = 1
          do j=1, n-1
             if ( saut(j) > 0 ) then
                sgn = -sgn
             end if
             if ( saut(j) /= 0 ) then
                k = p(j)
                p(j) = p(j+saut(j))
                p(j+saut(j)) = k
             end if
          end do
          !--- calcul du terme correspondant
          x =1.
          do j=1, n
             x = x *M(p(j),j)
          end do
          d = d + x*sgn
          !--- incrementation de la ref de la permutation
          do j=1,n-1
             if ( saut(n-j) < j ) then
                saut(n-j) = saut(n-j) + 1
                do k=1, j-1
                   saut(n-k) = 0
                end do
                exit
             end if
          end do
!!$--- fin de la boucle sur les permutations
       end do
!!$
end subroutine det
!!$
!!$***********************************************************************
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$***********************************************************************
!!$
subroutine inv_mat(IM,M,n)
!!$
!!$***********************************************************************
!!$
!!$          Routine d'inversion d'une matrice
!!$          (avec les la comatrices, pas optimisée ...)
!!$
!!$***********************************************************************
!!$
!!$----- Déclarations  
!!$
       implicit none
       integer, intent(in) :: n
       real*8, dimension(n,n), intent(in) :: M
       real*8, dimension(n,n), intent(out) :: IM
       real*8, dimension(n,n) :: CM
       real*8, dimension(n-1,n-1) :: PM
       real*8 :: dM
       integer :: i, j, ic, jc
!!$       
!!$--- det :
       call det(dM,M,n)
!!$--- comatrice :
       do j= 1, n
          do i= 1, n
             ! calc comatrice (lent...)
             do jc=1, j-1
                PM(1:i-1,jc)= M(1:i-1,jc)
                PM(i:n-1,jc)= M(i+1:n,jc)
             end do
             do jc=j+1, n
                PM(1:i-1,jc-1)= M(1:i-1,jc)
                PM(i:n-1,jc-1)= M(i+1:n,jc)
             end do
             ! le det
             call det(CM(i,j),PM,n-1)
             if ( btest(i+j,0) )  CM(i,j)= -CM(i,j)
          end do
       end do
!!$--- matrice inverse :
       IM= transpose(CM)/dM
!!$
end subroutine inv_mat
!!$
!!$***********************************************************************
!!$
