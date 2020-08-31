!***********************************************************************
!
!                Routine de mise en place d'une symétrie d'axe 3 
!                                autour du centre
!                 Marie-Bernadette Lepetit  - Octobre 2015
!
!
  subroutine axe_sym3(natom,atom_env,az_env,pos_env,q_env,&
       lim2,&
       natom_axis,atom_axe3,az_axe3,pos_axe3,q_axe3)
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!----- Déclarations  
!
    implicit none
! Variables globales
    real*8, parameter :: pi=3.14159265358979d0 
    real*8 :: lim2
    logical, parameter :: debug=.false.

! nombre d'atomes en entrée 
    integer :: natom
! atomes, positions fractionnaires, num atomique, charges en entrée
    character*3, dimension(natom) :: atom_env ! nom des atomes
    integer, dimension(natom) :: az_env   ! numeros atomiques 
    real*8, dimension(3,natom) :: pos_env !positions fractionnaires
    real*8, dimension(natom) :: q_env ! charge renormalisée

! atomes, etc... intermédiaires
    character*3, dimension(natom,2:3) :: atom_tmp ! nom des atomes
    integer,     dimension(natom,2:3) :: az_tmp   ! numeros atomiques 
    real*8,    dimension(3,natom,2:3) :: pos_tmp !positions fractionnaires
    real*8,      dimension(natom,2:3) :: q_tmp ! charge renormalisée
    integer,     dimension(3*natom)   :: deg_axis3

! atomes, etc... en sortie
    integer :: natom_axis
    character*3, dimension(3*natom) :: atom_axe3! nom des atomes
    integer, dimension(3*natom)     :: az_axe3   ! numeros atomiques 
    real*8, dimension(3,3*natom)    :: pos_axe3 !positions fractionnaires
    real*8, dimension(3*natom)      :: q_axe3 ! charge renormalisée

! variables muettes
    integer :: i, j, iatom, natom_tmp, ideg
    real*8 :: cc, ss, dist, deg

!    write (6,*) ">>>>>>> entree axe_sym3 <<<<<<<"
!-----------------------------------------------------------------------
!
!----- theta = 2 Pi/3
!
    atom_tmp(1:natom, 2) = atom_env(1:natom)
    az_tmp(1:natom, 2)   = az_env(1:natom)
    q_tmp(1:natom, 2) = q_env(1:natom)

    pos_tmp(3, 1:natom, 2) = pos_env(3, 1:natom)
    do i =1, natom
       pos_tmp(1,i,2) =  - pos_env(2,i)
       pos_tmp(2,i,2) =    pos_env(1,i) - pos_env(2,i)
    end do

!-----------------------------------------------------------------------
!
!----- theta = -2 Pi/3 = 4 Pi/3
!
    atom_tmp(1:natom, 3) = atom_env(1:natom)
    az_tmp(1:natom, 3)   = az_env(1:natom)
    q_tmp(1:natom, 3)    = q_env(1:natom)

    pos_tmp(3, 1:natom, 3) = pos_env(3, 1:natom)
    do i =1, natom
       pos_tmp(1,i,3) =  - pos_env(1,i) + pos_env(2,i)
       pos_tmp(2,i,3) =  - pos_env(1,i) 
    end do

!-----------------------------------------------------------------------
!
!----- Identification des atomes indépendants sur les trois ensembles
!
    
    atom_axe3(1:natom)     = atom_env(1:natom)
    az_axe3(1:natom)       = az_env(1:natom)
    pos_axe3(1:3, 1:natom) = pos_env(1:3, 1:natom)
    q_axe3(1:natom)        = q_env(1:natom)
    
    deg_axis3(1:natom)         = 1
    deg_axis3(1+natom:3*natom) = 0

    ! comparaisons pour les atomes tournes de 2 pi/3
    natom_axis = natom
    do j = 1,natom  ! tmp
       do i = 1,natom  ! ref
          if (abs(pos_tmp(1,j,2) - pos_axe3(1,i)) .lt. lim2 .and. &
              abs(pos_tmp(2,j,2) - pos_axe3(2,i)) .lt. lim2 .and. &
              abs(pos_tmp(3,j,2) - pos_axe3(3,i)) .lt. lim2 ) then 
              if (atom_tmp(j,2).ne.atom_axe3(i)) then
                write(6,*)
                write(6,*)
                write(6,'( " Atome ", a3,  3(F8.6,2x), " est a la meme position que atome ",&
                     a3, 3(F8.6,2x) )') &
                     atom_tmp(j,2),pos_tmp(1:3,j,2), atom_axe3(i), pos_axe3(1:3,i)
                write(6,*)
                stop "******* ERREUR *********"
             end if
             ! identique à atome 1
             deg_axis3(i) = deg_axis3(i) +1
             q_axe3(i)    = q_axe3(i) + q_tmp(j,2)
             if (debug) then 
                write(6,*) atom_axe3(i), pos_axe3(1:3,i),q_axe3(i) 
                write(6,*) atom_tmp(j,2), pos_tmp(1:3,j,2),q_tmp(j,2)
             end if
             goto 10
          else 
             iatom = j
          end if
       end do
       ! nouvel atome
       natom_axis               = natom_axis + 1
       atom_axe3(natom_axis)    = atom_tmp(j,2)
       az_axe3(natom_axis)      = az_tmp(j,2)
       pos_axe3(1:3,natom_axis) = pos_tmp(1:3,j,2)
       q_axe3(natom_axis)       = q_tmp(j,2)
       deg_axis3(natom_axis)    = 4 ! premier atom 2
!       write(6,*) atom_axe3(natom_axis),pos_axe3(1:3,natom_axis),q_axe3(natom_axis), deg_axis3(natom_axis)
10     continue
    end do
    !if (debug) &
         write(6,*) "Nombre d'atomes ajoute apres rotation de 2Pi/3 :", natom_axis-natom 

    ! comparaisons pour les atomes tournes de 4 pi/3
    natom_tmp = natom_axis 
    do j = 1,natom  ! tmp
       do i = 1,natom_tmp  ! ref
          if (abs(pos_tmp(1,j,3) - pos_axe3(1,i)) .lt. lim2 .and. &
              abs(pos_tmp(2,j,3) - pos_axe3(2,i)) .lt. lim2 .and. &
              abs(pos_tmp(3,j,3) - pos_axe3(3,i)) .lt. lim2 ) then 
             if (atom_tmp(j,3).ne.atom_axe3(i)) then
                write(6,*)
                write(6,*)
                write(6,'( " Atome ",A3,3(F8.6,2x)," est a la meme position que atome ",A3,3(F8.6,2x) )') &
                     atom_tmp(j,3), pos_tmp(1:3,j,3), atom_axe3(i), pos_axe3(1:3,i)
                write(6,*)
                stop "******* ERREUR *********"
             end if
             ! identique à atome 1 ou 2
             deg_axis3(i) = deg_axis3(i) +1
             q_axe3(i)    = q_axe3(i) + q_tmp(j,3)
             if (debug) then 
                write(6,*) atom_axe3(i),  pos_axe3(1:3,i),q_axe3(i) 
                write(6,*) atom_tmp(j,3), pos_tmp(1:3,j,3),q_tmp(j,3)
             end if
             goto 20
          else 
             iatom = i
          end if
       end do
       ! nouvel atome
       natom_axis               = natom_axis + 1
       atom_axe3(natom_axis)    = atom_tmp(j,3)
       az_axe3(natom_axis)      = az_tmp(j,3)
       pos_axe3(1:3,natom_axis) = pos_tmp(1:3,j,3)
       q_axe3(natom_axis)       = q_tmp(j,3)
       deg_axis3(natom_axis)    = 7 ! premier atom 3
!       write(6,*) atom_axe3(natom_axis),pos_axe3(1:3,natom_axis),q_axe3(natom_axis), deg_axis3(natom_axis)
20     continue
    end do

    !if (debug) &
    write(6,*) "Nombre d'atomes ajoute apres rotation de 4Pi/3 :",&
         natom_axis-natom_tmp, natom_axis-natom 

!-----------------------------------------------------------------------
!
!----- moyenne des charges
!
    do i=1, natom_axis
       if (deg_axis3(i).eq.0) stop  "******* ERREUR deg *****"
       ideg = mod(deg_axis3(i)-1,3)+1
       deg  = dble(ideg)
       q_axe3(i) = q_axe3(i) / deg 
!      write(6,*) i, deg
    end do

    return
  end subroutine axe_sym3


!!$***********************************************************************
!!$***********************************************************************
!***********************************************************************
!
!                Routine de test d'une symétrie d'axe 3 
!                           autour du centre
!                 Marie-Bernadette Lepetit  - Octobre 2015
!
!
  subroutine test_sym3(natom,atom_env,coord_env,q_env,lim2,lim4)
! 
!***********************************************************************
!***********************************************************************
!
!----- Déclarations  
!
    implicit none
! Variables globales
    real*8, parameter :: pi=3.14159265358979d0 
    real*8 :: lim2, lim4

! nombre d'atomes en entrée 
    integer :: natom
! atomes, positions fractionnaires, num atomique, charges en entrée
    character*3, dimension(natom) :: atom_env ! nom des atomes
!    integer, dimension(natom) :: az_env   ! numeros atomiques 
!    real*8, dimension(3,natom) :: pos_env !positions fractionnaires
    real*8, dimension(3,natom) :: coord_env !positions cartésiennes 
    real*8, dimension(natom) :: q_env ! charge renormalisée

! atomes, etc... intermédiaires
    character*3, dimension(natom)     :: atom_tmp ! nom des atomes
!    integer,     dimension(natom,2:3) :: az_tmp   ! numeros atomiques 
!    real*8,    dimension(3,natom,2:3) :: pos_tmp !positions fractionnaire
    real*8,      dimension(3,natom)   :: coord_tmp !positions cartésiennes 
    real*8,      dimension(natom)     :: q_tmp ! charge renormalisée


! variables muettes
    integer :: i, j, iatom
    real*8 :: cc, ss

!***********************************************************************
    cc = cos(2.d0*pi/3.d0)
    ss = sin(2.d0*pi/3.d0)

    write(6,*)
    write(6,*) '---------------------------------------------------'
    write(6,*) '       Entre dans test de symmetrie ordre 3        '
    write(6,*) 
    write(6,*) ' Les atomes suivants sont symetriques mais de charge differentes '

!-----------------------------------------------------------------------
!
!----- Rotation de theta = 2 Pi/3
!
    atom_tmp(1:natom)    = atom_env(1:natom) 
    q_tmp(1:natom)       = q_env(1:natom) 
    coord_tmp(3,1:natom) =  coord_env(3,1:natom) 
    do i=1,natom
       coord_tmp(1,i) =  cc*coord_env(1,i) + ss*coord_env(2,i)  
       coord_tmp(2,i) = -ss*coord_env(1,i) + cc*coord_env(2,i)  
    end do

!    write(6,*)
!    do i=1,natom
!       write(6,*) atom_env(i), coord_env(1:3,i), q_env(i) 
!       write(6,*) atom_tmp(i), coord_tmp(1:3,i), q_tmp(i) 
!       write(6,*)
!    end do
!-----------------------------------------------------------------------
!
!----- Test de charge 
!

    iatom = 0
    do j = 1,natom  ! tmp
       do i = 1,natom  ! ref
          if (abs(coord_tmp(1,j) - coord_env(1,i)) .lt. lim2 .and. &
              abs(coord_tmp(2,j) - coord_env(2,i)) .lt. lim2 .and. &
              abs(coord_tmp(3,j) - coord_env(3,i)) .lt. lim2 ) then 
             ! atomes identiques 
             if (atom_tmp(j).ne.atom_env(i)) then
                write(6,*)
                write(6,*)
                write(6,'( " Atome ", a3,  3(F8.6,2x), " est a la meme position que atome ",&
                     a3, 3(F8.6,2x) )') &
                     atom_tmp(j),coord_tmp(1:3,j), atom_env(i), coord_env(1:3,i)
                write(6,*)
                stop "******* ERREUR *********"
             end if
             ! test de la charge
             if (abs(q_tmp(j) - q_env(i)) .gt. lim4) then 
                ! charges differentes
                iatom = iatom +1
                if (iatom .le. 10) &
                     write(6,9001) atom_env(i), coord_env(1:3,i), q_env(i), &
                     atom_tmp(j), coord_tmp(1:3,j), q_tmp(j), &
                     q_env(i)-q_tmp(j)
             else
                goto 10
             end if
          end if
       end do
10     continue
    end do

    write(6,*)
    if (iatom.gt.10) write(6,*) ' .......' 
    write(6,*)
    write(6,9002) '+',iatom

!-----------------------------------------------------------------------
!
!----- Rotation de theta = -2 Pi/3
!
    atom_tmp(1:natom)    = atom_env(1:natom) 
    q_tmp(1:natom)       = q_env(1:natom) 
    coord_tmp(3,1:natom) =  coord_env(3,1:natom) 
    do i=1,natom
       coord_tmp(1,i) =  cc*coord_env(1,i) - ss*coord_env(2,i)  
       coord_tmp(2,i) =  ss*coord_env(1,i) + cc*coord_env(2,i)  
    end do

!    write(6,*)
!    do i=1,natom
!       write(6,*) atom_env(i), coord_env(1:3,i), q_env(i) 
!       write(6,*) atom_tmp(i), coord_tmp(1:3,i), q_tmp(i) 
!       write(6,*)
!    end do
!-----------------------------------------------------------------------
!
!----- Test de charge 
!
    iatom = 0
    do j = 1,natom  ! tmp
       do i = 1,natom  ! ref
          if (abs(coord_tmp(1,j) - coord_env(1,i)) .lt. lim2 .and. &
              abs(coord_tmp(2,j) - coord_env(2,i)) .lt. lim2 .and. &
              abs(coord_tmp(3,j) - coord_env(3,i)) .lt. lim2 ) then 
             ! atomes identiques 
             if (atom_tmp(j).ne.atom_env(i)) then
                write(6,*)
                write(6,*)
                write(6,'( " Atome ", a3,  3(F8.6,2x), " est a la meme position que atome ",&
                     a3, 3(F8.6,2x) )') &
                     atom_tmp(j),coord_tmp(1:3,j), atom_env(i), coord_env(1:3,i)
                write(6,*)
                stop "******* ERREUR *********"
             end if
             ! test de la charge
             if (abs(q_tmp(j) - q_env(i)) .gt. lim4) then 
                ! charges differentes
                iatom = iatom +1
                if (iatom .le. 10) &
                     write(6,9001) atom_env(i), coord_env(1:3,i), q_env(i), &
                     atom_tmp(j), coord_tmp(1:3,j), q_tmp(j), &
                     q_env(i)-q_tmp(j)
             else
                goto 20
             end if
          end if
       end do
20     continue
    end do

    write(6,*)
    if (iatom.gt.10) write(6,*) ' .......' 
    write(6,*)
    write(6,9002) '-',iatom
 
9001 format (1x,a3,4(2x,f10.4),5x,'et',5x,a3,4(2x,f10.4),5x,'Delta Q =',5x,D12.5)
9002 format (1x,a1,"2 Pi/3",2x,i10, 2x, 'atomes sont symetriques et de charges differentes')
    return
  end subroutine test_sym3
