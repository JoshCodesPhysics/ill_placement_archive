!***********************************************************************
!
!   Routine de calcul des fonctions mu de renormalisation des charges
!                 Alain Gellé - Octobre 2004
!
!
subroutine  fctmu(nc,p,pos,lambda)
       !
       !  
       !***********************************************************************
       !***********************************************************************
       !
       !----- Déclarations  
       implicit none

       integer, intent(in) :: p     ! nbr de moments annulés = p-1
       integer, intent(in) :: nc    ! nbr de charges de la maille élémentaire

       ! positions fractionnaires des atomes dans la maille élémentaire
       real*8, dimension(3,nc), intent(in) :: pos

       ! variables de renormalisation des charges (n° de l'atome, coordonnée, cellule)
       real*8, dimension(nc,3,p), intent(out) :: lambda

       ! correspond à (voir thèse) :  lambda0(i) = lambda_i(0)
       real*8, dimension(:), allocatable :: lambda0

       ! correpsond à (voir thèse) :  lambda1(i) = lambda_1(i-1)
       real*8, dimension(:), allocatable :: lambda1


       integer :: i, j, k, l, ic, id, ip
       real*8 :: x,sgn,a
       !
       !***********************************************************************  
       !
       !----- Initialisations 
       !
       write(6,*) 
!       write(6,*) 'DEBUG fctmu'
!       do i=1,nc
!          write(6,'(3(2x,f10.5))') pos(:,i)
!       end do

       allocate (lambda0(p), lambda1(p))

       !----- Calcul de lambda_1
       ! lambda1(i) = lambda_1(i-1) = (i-1)^(p-1) / (p-1)!

       lambda1(:) = 1.0d0
       lambda1(1) = 0.0d0
       do i= 2, p
          do j= 1, p-1
             lambda1(i) = lambda1(i)*dble(i-1)/dble(j)
          end do
       end do


       !----- calcul de lambda_0
       ! (notations these)
       !   lambda0(i) = lambda_i(0) = sum(j=0,i-1) (-1)^j binom(p,j) lambda_1(i-j-1)
       ! (ici)
       !   lambda0(i) = sum(j=0,i-2) (-1)^j binom(p,j) lambda1(i-j)
       !   car pour j=i-1, lambda1(1) = 0

       lambda0(:) = 0.0d0
       do i= 1, p
          do j= 0, i-2
             sgn= 1.d0
             if (btest(j,0)) sgn = -1.d0
             lambda0(i) = lambda0(i) + sgn*binom(p,j)*lambda1(i-j)
          end do
       end do


       !----- calcul de lambda
       ! lambda_i(x) = sum(l=1,p) lambda_j(0) prod(j=1,p,\=i) (x+j-l)/(j-i) ??
       ! lambda_i(x) = sum(l=1,p) lambda_l(0) prod(j=1,p,\=i) (x+j-l)/(j-i) 

       lambda(:,:,:) = 0.0d0
       do ip = 1, p
          !--- calcul pour chaque atome
          do ic= 1, nc
             !--- calcul pour chaque dimensions
             do id = 1, 3
                x = pos(id,ic)
                !--- calcul lambda_i(x)
                do l=2, p
                   a = lambda0(l)
                   do j=1, p
                      if (j/=ip) then
                         a= a*(x+dble(j-l))/dble(j-ip)
                      end if
                   end do
                   lambda(ic,id,ip)= lambda(ic,id,ip) + a
                end do
             end do
          end do
       end do

contains

       function binom(a,b)
              implicit none
              real*8 :: binom
              integer, intent(in) :: a, b
              integer :: ii
              binom = 1.0d0
              do ii= 1, b
                 binom = binom * dble(a+ii-b) / dble(ii)
              end do
       end function binom

end subroutine fctmu

