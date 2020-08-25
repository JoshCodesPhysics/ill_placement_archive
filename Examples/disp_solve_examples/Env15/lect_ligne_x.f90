! version 3
!

!***********************************************************************
!
!   Routines de de lecture dans un ligne 
!                 Alain Gellé -  2007
!
!
!
!
!***********************************************************************

!***********************************************************************
!
!      donees partagees ligne a lire et position sur la ligne
!
!***********************************************************************

module ligne_data
       integer, parameter :: lmaxl=200
       character(len=lmaxl), save :: ligne
       integer, save :: l_p1, l_p2
end module ligne_data

!***********************************************************************
!
!      interfaces
!
!***********************************************************************

module interface_lect_ligne
       implicit none
       interface
          !
          !--- nouvelle ligne=a
          subroutine ligne_nouvelle(a)
                 use ligne_data
                 character(len=*), intent(in) :: a
          end subroutine ligne_nouvelle
          !
          !--- se mettre au debut
          subroutine ligne_debut()
                 use ligne_data
          end subroutine ligne_debut
          !
          !--- paquet suivant
          subroutine ligne_suivant()
                 use ligne_data
          end subroutine ligne_suivant
          !
          !--- nbr de paquet 
          subroutine ligne_n(n)
                 use ligne_data
          end subroutine ligne_n
          !
          !--- rechercher un mot clef
          subroutine ligne_clef(clei,OK)
                 use ligne_data
                 logical, intent(out) :: OK
                 character(len=*), intent(in) :: clei
          end subroutine ligne_clef
          !
          !--- lecture entier
          subroutine ligne_int(vi,OK)
                 use ligne_data
                 integer, intent(out) :: vi
                 logical, intent(out) ::  OK
          end subroutine ligne_int
          !
          !--- idem sinon stop
          subroutine ligne_int_stop(vi,mess)
                 use ligne_data
                 integer, intent(out) :: vi
                 character(len=*), intent(in) ::  mess
          end subroutine ligne_int_stop
          !
          !--- lecture tableau d'entiers
          subroutine ligne_tint(ti,n,OK)
                 use ligne_data
                 integer, intent(in) :: n
                 integer, dimension(n), intent(out) :: ti
                 logical, intent(out) ::  OK
          end subroutine ligne_tint
          !
          !--- idem sinon stop
          subroutine ligne_tint_stop(ti,n,mess)
                 use ligne_data
                 integer, intent(in) :: n
                 integer, dimension(n), intent(out) :: ti
                 character(len=*), intent(in) ::  mess
          end subroutine ligne_tint_stop
          !
          !--- lecture d'un reel 8
          subroutine ligne_r8_stop(vr,mess)
                 use ligne_data
                 real*8, intent(out) :: vr
                 character (len=*), intent(in) ::  mess
          end subroutine ligne_r8_stop
          !
          !--- lecture d'un reel 8
          subroutine ligne_r8(vr,OK)
                 use ligne_data
                 real*8, intent(out) :: vr
                 logical, intent(out) ::  OK
          end subroutine ligne_r8
          !
          !--- lecture d'un tableau de reel 8
          subroutine ligne_tr8_stop(tvr,n,mess)
                 use ligne_data
                 integer, intent(in) :: n
                 real*8, dimension(n), intent(out) :: tvr
                 character (len=*), intent(in) ::  mess
          end subroutine ligne_tr8_stop
          !
          !--- lecture d'une chaine
          subroutine ligne_ch(a,OK)
                 use ligne_data
                 character (len=*), intent(out) :: a
                 logical, intent(out) ::  OK
          end subroutine ligne_ch
          !
          !--- lecture d'une chaine
          subroutine ligne_ch_stop(a,mess)
                 use ligne_data
                 character (len=*), intent(out) :: a
                 character (len=*) ::  mess
          end subroutine ligne_ch_stop
          !
          !--- trouver une chaine
          subroutine ligne_trouver(a,OK)
                 use ligne_data
                 character (len=*), intent(in) :: a
                 logical, intent(out) ::  OK
          end subroutine ligne_trouver
          !
          !--- remplace les caract non numerique par des espaces
          subroutine ligne_numonly()
                 use ligne_data
          end subroutine ligne_numonly
       end interface
end module interface_lect_ligne


!***********************************************************************
!
!      nouvelle ligne
!
!***********************************************************************

subroutine ligne_nouvelle(a)
       use ligne_data
       implicit none
       character(len=*), intent(in) :: a

       ligne=""
       if ( len_trim(a)<=lmaxl ) then
          ligne=trim(a)
       else
          write(6,"('**** erreur (ligne) : nouvelle ligne trop grande ****')")
          write(6,"(a)") trim(a)
          ligne=a(1:lmaxl)
       end if
       l_p1= 0
       l_p2= 0
end subroutine ligne_nouvelle

!***********************************************************************
!
!      position au debut de la ligne
!
!***********************************************************************

subroutine ligne_debut()
       use ligne_data
       implicit none

       l_p1= 0
       l_p2= 0
end subroutine ligne_debut


!***********************************************************************
!
!      position d'un nouveau paquet de characters
!
!***********************************************************************

subroutine ligne_suivant()
       use ligne_data
       implicit none
       
       do l_p1= l_p2+1, lmaxl
          if ( ligne(l_p1:l_p1)/=" " ) exit
       end do
       do l_p2= l_p1+1, lmaxl
          if ( ligne(l_p2:l_p2)==" " ) exit
       end do
       if ( l_p2>lmaxl )  l_p2=lmaxl
end subroutine ligne_suivant

!***********************************************************************
!
!      comptage du nombre de paquet
!
!***********************************************************************

subroutine ligne_n(n)
       use ligne_data
       implicit none
       integer, intent(out) :: n
       integer :: l_p1s, l_p2s
       
       n= 0
       l_p1s= l_p1
       l_p2s= l_p2
       do
          call ligne_suivant()
          if ( l_p1>lmaxl )  exit
          n= n+1
       end do
       l_p1= l_p1s
       l_p2= l_p2s
end subroutine ligne_n


!***********************************************************************
!
!      recherche d'un mot clef dans une ligne
!
!***********************************************************************

subroutine ligne_clef(clei,OK)
       use ligne_data
       implicit none
       logical, intent(out) :: OK
       character(len=*), intent(in) :: clei
       character(len=len(clei)) :: cle
       character(len=lmaxl) ::lignem
       integer :: l
       
       OK= .false.
       cle= clei
       call lett_min(cle)
       l= len_trim(cle)
       lignem = adjustl(ligne)
       call lett_min(lignem)
       if ( cle(:l)==lignem(:l) )  OK= .true.
end subroutine ligne_clef

!***********************************************************************
!
!      lecture d'un entier
!
!***********************************************************************

subroutine ligne_int(vi,OK)
       use ligne_data
       implicit none
       integer, intent(out) :: vi
       logical, intent(out) ::  OK
       integer :: ierr, i
       
       OK= .false.
       call ligne_suivant()
       if ( l_p1>lmaxl ) return
       read(ligne(l_p1:l_p2),*,iostat=ierr) vi
       if ( ierr==0 )  OK= .true.
end subroutine ligne_int

!***********************************************************************
!
!      lecture d'un entier
!      avec arret si echec  (message "mess")
!
!***********************************************************************

subroutine ligne_int_stop(vi,mess)
       use ligne_data
       implicit none
       integer, intent(out) :: vi
       character(len=*), intent(in) ::  mess
       integer :: ierr, i
       
       call ligne_suivant()
       if ( l_p1>lmaxl ) then
          write(6,"('**** erreur (ligne) : lecture d''un entier, ligne vide ****')")
          write(6,"('**** ',a,' ****')") mess
          stop
       end if
       read(ligne(l_p1:l_p2),*,iostat=ierr) vi
       if ( ierr/=0 ) then
          write(6,"('**** erreur (ligne) : lecture d''un entier ****')")
          write(6,"('**** ',a,' ****')") mess
          stop
       end if
end subroutine ligne_int_stop

!***********************************************************************
!
!      lecture d'un tableau d'entier
!
!***********************************************************************

subroutine ligne_tint(ti,n,OK)
       use ligne_data
       implicit none
       integer, intent(in) :: n
       integer, dimension(n), intent(out) :: ti
       logical, intent(out) ::  OK
       integer :: ierr, i
       
       do i= 1, n
          call ligne_int(ti(i),OK)
          if ( .not. OK )  return
       end do
end subroutine ligne_tint

!***********************************************************************
!
!      lecture d'un tableau d'entier
!      avec arret si echec  (message "mess")
!
!***********************************************************************

subroutine ligne_tint_stop(ti,n,mess)
       use ligne_data
       implicit none
       integer, intent(in) :: n
       integer, dimension(n), intent(out) :: ti
       character(len=*), intent(in) ::  mess
       integer :: ierr, i
       character(len=8) :: bla
       
       do i= 1, n
          write(bla,"(i8)") i
          bla= adjustl(bla)
          call ligne_int_stop(ti(i),mess//' - indice '//trim(bla))
       end do
end subroutine ligne_tint_stop

!***********************************************************************
!
!      lecture d'un reel 8
!      avec arret si echec  (message "mess")
!
!***********************************************************************

subroutine ligne_r8_stop(vr,mess)
       use ligne_data
       implicit none
       real*8, intent(out) :: vr
       character (len=*), intent(in) ::  mess
       integer :: ierr, i

       call ligne_suivant()
       if ( l_p1>lmaxl ) then
          write(6,"('**** erreur (ligne) : lecture d''un entier, ligne vide ****')")
          write(6,"('**** ',a,' ****')") mess
          stop
       end if
       read(ligne(l_p1:l_p2),*,iostat=ierr) vr
       if ( ierr/=0 ) then
          write(6,"('**** erreur (ligne) : lecture d''un entier ****')")
          write(6,"('**** ',a,' ****')") mess
          stop
       end if
!!$       write(6,"(a,' : ',ES12.5)") mess, vr
end subroutine ligne_r8_stop

!***********************************************************************
!
!      lecture d'un reel 8
!      
!
!***********************************************************************

subroutine ligne_r8(vr,OK)
       use ligne_data
       implicit none
       real*8, intent(out) :: vr
       logical, intent(out) ::  OK
       integer :: ierr, i

       OK= .false.
       call ligne_suivant()
       if ( l_p1>lmaxl ) return
       read(ligne(l_p1:l_p2),*,iostat=ierr) vr
       if ( ierr/=0 ) return
       OK= .true.
end subroutine ligne_r8

!***********************************************************************
!
!      lecture d'un tableau de reel 8
!      avec arret si echec  (message "mess")
!
!***********************************************************************

subroutine ligne_tr8_stop(tvr,n,mess)
       use ligne_data
       implicit none
       integer, intent(in) :: n
       real*8, dimension(n), intent(out) :: tvr
       character (len=*), intent(in) ::  mess
       integer :: ierr, i

       do i=1, n
          call ligne_r8_stop(tvr(i),mess)
       end do
end subroutine ligne_tr8_stop

!***********************************************************************
!
!      lecture d'une ch. de caract.,
!      
!
!***********************************************************************

subroutine ligne_ch(a,OK)
       use ligne_data
       implicit none
       character (len=*), intent(out) :: a
       logical, intent(out) :: OK
       integer :: ierr
       
       a= ""
       OK= .false.
       call ligne_suivant()
       if ( l_p1>lmaxl ) return
       read(ligne(l_p1:l_p2),*,iostat=ierr) a
       if ( ierr/=0 ) then
          a= ""
          return
       end if
       a= adjustl(a)
       OK= .true.
       
end subroutine ligne_ch


!***********************************************************************
!
!      lecture d'une ch. de caract.,
!      avec arret si echec  (message "mess")
!
!***********************************************************************

subroutine ligne_ch_stop(a,mess)
       use ligne_data
       implicit none
       character (len=*), intent(out) :: a
       character (len=*) ::  mess
       integer :: ierr
       
       call ligne_suivant()
       if ( l_p1>lmaxl ) then
          write(6,"('**** erreur (ligne) : lecture d''une ch. de caract., ligne vide ****')")
          write(6,"('**** ',a,' ****')") mess
          stop
       end if
       read(ligne(l_p1:l_p2),*,iostat=ierr) a
       if ( ierr/=0 ) then
          write(6,"('**** erreur (ligne) : lecture d''une ch. de caract. ****')")
          write(6,"('**** ',a,' ****')") mess
          stop
       end if
       a= adjustl(a)
!!$       write(6,"(x,a,' : ',a)") mess, a
end subroutine ligne_ch_stop

!***********************************************************************
!
!      trouver une chaine
!      
!
!***********************************************************************

subroutine ligne_trouver(a,OK)
       use ligne_data
       implicit none
       character (len=*), intent(in) :: a
       logical, intent(out) ::  OK
       integer :: ierr, i
       
       OK= .false.
       call ligne_suivant()
       if ( l_p1>lmaxl ) then
          return
       end if
       i= index(ligne(l_p1:),a)
       if ( i>0 ) then
          OK= .true.
          l_p1= i
          do l_p2= l_p1+1, lmaxl
             if ( ligne(l_p2:l_p2)==" " ) exit
          end do
          if ( l_p2>lmaxl )  l_p2=lmaxl
       end if
end subroutine ligne_trouver


!***********************************************************************
!
!      remplacer les caract non numerique par des espaces
!      
!
!***********************************************************************

subroutine ligne_numonly()
       use ligne_data
       implicit none
       integer :: i, l
       
       l=len_trim(ligne)
       do i= 1, l
          if ( index('+-0123456789',ligne(i:i))==0 ) then
             ligne(i:i)=' '
          end if
       end do
end subroutine ligne_numonly
