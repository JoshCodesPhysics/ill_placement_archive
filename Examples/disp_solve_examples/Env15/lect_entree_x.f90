! version 3
!

!***********************************************************************
!
!   Routines de de lecture des entrées  -  format :  clef  /  valeur    
!                 Alain Gellé - Octobre 2006
!
!
!
!
!***********************************************************************

!***********************************************************************
!
!      donees partagees, unite et nom du fichier tmp
!
!***********************************************************************

module entree_data
       ! caracteristique fichier tmp
       integer, parameter :: lmaxfe=50
       integer, save :: ufe
       character(len=lmaxfe), save :: fe
       logical, save :: okfe= .false.
       ! longueur des mots clef
       integer, parameter :: lmaxclef=50
       ! clef fin d'entree
       character(len=*), parameter :: cle_fin='end', cle_fin0='&end'
end module entree_data

!***********************************************************************
!
!      interfaces
!
!***********************************************************************

module interface_lect_entree
       implicit none
       interface
          !
          !--- creation du fichier tmp
          subroutine  entree_init(n,nf)
                 use entree_data
                 use ligne_data
                 integer, intent(in) :: n
                 character(len=*), intent(in) :: nf
          end subroutine entree_init
          !
          !--- idem avec mot clef debut/fin
          subroutine entree_init_clef(n,nf,cle_deb)
                 use entree_data
                 use ligne_data
                 integer, intent(in) :: n
                 character(len=*), intent(in) :: nf, cle_deb
          end subroutine entree_init_clef
          !
          !--- fermer les entrees
          subroutine  entree_close()
                 use entree_data
          end subroutine entree_close
          !
          !--- chercher une clef
          subroutine entree_clef(clei,OK)
                 use entree_data
                 use ligne_data
                 logical, intent(out) :: OK
                 character(len=*), intent(in) :: clei
          end subroutine entree_clef
          !
          !--- idem sinon stop
          subroutine entree_clef_stop(cle,mess)
                 use entree_data
                 character (len=*), intent(in) :: cle, mess
          end subroutine entree_clef_stop
          !
          !--- envoie une ligne de donnees vers lect_ligne
          subroutine entree_ligne()
                 use entree_data
                 use ligne_data
          end subroutine entree_ligne
          !
          !--- idem sinon stop
          subroutine entree_ligne_stop()
                 use entree_data
                 use ligne_data
          end subroutine entree_ligne_stop
          !
          !--- lecture entier
          subroutine entree_int_stop(vi,mess)
                 use entree_data
                 integer, intent(out) :: vi
                 character (len=*), intent(in) ::  mess
          end subroutine entree_int_stop
          !
          !--- lecture d'un tableau
          subroutine entree_tint_stop(tvi,d,mess)
                 use entree_data
                 integer, intent(in) :: d
                 integer, dimension(d), intent(out) :: tvi
                 character (len=*), intent(in) ::  mess
          end subroutine entree_tint_stop
          !
          !--- lecture reel 8
          subroutine entree_r8_stop(vi,mess)
                 use entree_data
                 real*8, intent(out) :: vi
                 character (len=*), intent(in) ::  mess
          end subroutine entree_r8_stop
          !
          !--- lecture tableau reel 8
          subroutine entree_tr8_stop(tvi,d,mess)
                 use entree_data
                 integer, intent(in) :: d
                 real*8, dimension(d), intent(out) :: tvi
                 character (len=*), intent(in) ::  mess
          end subroutine entree_tr8_stop
          !
          !--- lecture chaine
          subroutine entree_ch_stop(a,mess)
                 use entree_data
                 character (len=*), intent(out) :: a
                 character (len=*), intent(in) ::  mess
          end subroutine entree_ch_stop
          !
          !--- remplace une chaine si on trouve la clef
          subroutine entree_repl_ch_clef(clei,a,mess)
                 use entree_data
                 character(len=*), intent(in) :: clei
                  character (len=*), intent(out) :: a
                 character (len=*), intent(in) ::  mess
          end subroutine entree_repl_ch_clef
          !
          !---
       end interface
end module interface_lect_entree

!***********************************************************************
!
!      recopie des entrees dans un fichier
!
!      fichier  "nom",   unité "n"
!
!***********************************************************************

subroutine  entree_init(n,nf)
       use entree_data
       use ligne_data
       implicit none
       integer, intent(in) :: n
       character(len=*), intent(in) :: nf
!!$       character(len=lmaxl) :: ligne
       integer :: ierr

       write(6,"(' lecture des entrees')")
       if ( okfe ) then
          write(6,"('**** erreur (entree) : un fichier d''entree est deja ouvert ****')")
          write(6,"('**** ancien :',a,' ,nouveau ',a,' ****')") trim(fe), trim(nf)
          stop
       end if
       if ( len_trim(nf) > lmaxfe ) then
          write(6,"('**** erreur (entree) : nom de fichier d''entree trop long ****')")
          write(6,"('**** ',i4,' caract. > ',i4,' ****')") len_trim(nf), lmaxfe
          stop
       end if
       open (unit=n,file=trim(nf),form="formatted",status="replace")
       do
          read(5,"(a)",iostat=ierr) ligne
          if ( ierr /= 0 ) exit
          write(n,"(a)") trim(ligne)
       end do
       rewind (n)
       ufe= n
       fe= trim(nf)
       okfe= .true.
end subroutine entree_init

!***********************************************************************
!
!      recopie des entrees dans un fichier, avec mot clef donnees
!
!      fichier  "nom",   unité "n"
!
!***********************************************************************

subroutine entree_init_clef(n,nf,cle_deb)
       use entree_data
       use ligne_data
       implicit none
       integer, intent(in) :: n
       character(len=*), intent(in) :: nf, cle_deb
       logical :: OK, OKD, OKF
       integer :: ierr

       write(6,"(' lecture des entrees')")
       if ( okfe ) then
          write(6,"('**** erreur (entree) : un fichier d''entree est deja ouvert ****')")
          write(6,"('**** ancien :',a,' ,nouveau ',a,' ****')") trim(fe), trim(nf)
          stop
       end if
       if ( len_trim(nf) > lmaxfe ) then
          write(6,"('**** erreur (entree) : nom de fichier d''entree trop long ****')")
          write(6,"('**** ',i4,' caract. > ',i4,' ****')") len_trim(nf), lmaxfe
          stop
       end if
       open (unit=n,file=trim(nf),form="formatted",status="replace")
       OKD= .false.
       OK= .false.
       do
          read(5,"(a)",iostat=ierr) ligne
          if ( ierr /= 0 ) then
             if (OKD)  write(6,"('**** pas de mot clef de fin : ',a,' ****')") trim(cle_fin)
             exit
          end if
          if (.not. OKD) then
             call ligne_clef(cle_deb,OK)
             if (OK) then
                write(6,"(' nom des entree : ',a)") trim(ligne)
                rewind(n)
                OKD= .true.
                OK= .false.
             end if
          end if
          write(n,"(a)") trim(ligne)
          call ligne_clef(cle_fin,OKF)
          if (OKF)  then
             if (OKD) exit
             write(6,"('**** mot clef de fin sans celui de debut (',a,' sans ',a,') ****')") &
                  trim(cle_fin), trim(cle_deb)
             write(6,"('**** mot clef ignore ****')")
          end if
          call ligne_clef(cle_fin0,OKF)
          if (OKF)  then
             if (OKD) exit
             write(6,"('**** mot clef de fin sans celui de debut (',a,' sans ',a,') ****')") &
                  trim(cle_fin), trim(cle_deb)
             write(6,"('**** mot clef ignore ****')")
          end if
       end do
       rewind (n)
       ufe= n
       fe= trim(nf)
       okfe= .true.
end subroutine entree_init_clef

!***********************************************************************
!
!      fermeture et effacement du fichier
!
!***********************************************************************

subroutine  entree_close()
       use entree_data
       implicit none
       
       close(ufe,status='delete')
end subroutine entree_close


!***********************************************************************
!
!      recherche d'un mot clef
!
!      dans unite "n",
!      independant minuscule/majuscules
!      resultat dans OK,   pointe sur ligne suivante
!
!***********************************************************************

subroutine entree_clef(clei,OK)
       use entree_data
       use ligne_data
       implicit none
       logical, intent(out) :: OK
       character(len=*), intent(in) :: clei
       character(len=len(clei)) :: cle
!!$       character(len=lmaxl) :: ligne
       integer :: ierr

       cle=adjustl(clei)
       call lett_min(cle)
       rewind (ufe)
       OK = .false.
       do
          ligne= ""
          read (ufe,"(a)",iostat=ierr ) ligne
          if ( ierr < 0 ) return
          call ligne_clef(cle,OK)
          if (OK)  return
       end do
end subroutine entree_clef
       
!***********************************************************************
!
!      recherche d'un mot clef
!
!      arret du prog en cas d'echec  (avec message "mess")
!
!***********************************************************************

subroutine entree_clef_stop(cle,mess)
       use entree_data
       implicit none
       character (len=*), intent(in) :: cle, mess
       logical :: OK
       
       call entree_clef(cle,OK)
       
       if (.not.OK) then
          write(6,"('**** erreur (entree) : mot clef non trouve ****')")
          write(6,"('**** ',a,' ****')") mess
          stop
       end if
end subroutine entree_clef_stop
       
!***********************************************************************
!
!      mise en memoire d'une ligne
!
!***********************************************************************

subroutine entree_ligne()
       use entree_data
       use ligne_data
       implicit none
       integer :: ierr

       read(ufe,"(a)",iostat=ierr) ligne
       if ( ierr/=0 ) then
          write(6,"('**** erreur (entree) : lecture d''une ligne ****')")
       end if
       call ligne_debut()
end subroutine entree_ligne

!***********************************************************************
!
!      mise en memoire d'une ligne, avec arret
!
!***********************************************************************

subroutine entree_ligne_stop()
       use entree_data
       use ligne_data
       implicit none
       integer :: ierr
       
       read(ufe,"(a)",iostat=ierr) ligne
       if ( ierr/=0 ) then
          write(6,"('**** erreur (entree) : lecture d''une ligne ****')")
          stop
       end if
       call ligne_debut()
end subroutine entree_ligne_stop

!***********************************************************************
!
!      lecture d'un entier, 
!      avec arret si echec  (message "mess")
!
!***********************************************************************

subroutine entree_int_stop(vi,mess)
       use entree_data
       implicit none
       integer, intent(out) :: vi
       character (len=*), intent(in) ::  mess
       integer :: ierr

       read(ufe,*,iostat=ierr) vi
       if ( ierr/=0 ) then
          write(6,"('**** erreur (entree) : lecture d''un entier ****')")
          write(6,"('**** ',a,' ****')") mess
          stop
       end if
       if ( vi < 9999 ) then
          write(6,"(x,a,' : ',i4)") mess, vi
       else 
          write(6,"(x,a,' : ',i12)") mess, vi
       end if
end subroutine entree_int_stop


!***********************************************************************
!
!      lecture d'un tableau d'entier, 
!      avec arret si echec  (message "mess")
!
!***********************************************************************

subroutine entree_tint_stop(tvi,d,mess)
       use entree_data
       implicit none
       integer, intent(in) :: d
       integer, dimension(d), intent(out) :: tvi
       character (len=*), intent(in) ::  mess
       integer :: ierr, l

       read(ufe,*,iostat=ierr) tvi
       if ( ierr/=0 ) then
          write(6,"('**** erreur (entree) : lecture d''un entier ****')")
          write(6,"('**** ',a,' ****')") mess
          stop
       end if
       l= len(mess)
       if ( maxval(tvi) < 9999 ) then
          if ( l+5*d <= 80 )  then
             write(6,"(x,a,' :',100(x,i4))") mess, tvi
          else
             write(6,"(x,a,' :',100(/4x,15(x,i4)))") mess, tvi
          end if
       else 
          if ( l+13*d <= 80 )  then
             write(6,"(x,a,' :',100(x,i12))") mess, tvi
          else
             write(6,"(x,a,' :',100(/4x,5(x,i12)))") mess, tvi
          end if
       end if
end subroutine entree_tint_stop

!***********************************************************************
!
!      lecture d'un reel 8, 
!      avec arret si echec  (message "mess")
!
!***********************************************************************

subroutine entree_r8_stop(vi,mess)
       use entree_data
       implicit none
       real*8, intent(out) :: vi
       character (len=*), intent(in) ::  mess
       integer :: ierr

       read(ufe,*,iostat=ierr) vi
       if ( ierr/=0 ) then
          write(6,"('**** erreur (entree) : lecture d''un reel ****')")
          write(6,"('**** ',a,' ****')") mess
          stop
       end if
       write(6,"(x,a,' : ',ES9.2)") mess, vi
end subroutine entree_r8_stop


!***********************************************************************
!
!      lecture d'un tableau de reels 8, 
!      avec arret si echec  (message "mess")
!
!***********************************************************************

subroutine entree_tr8_stop(tvi,d,mess)
       use entree_data
       implicit none
       integer, intent(in) :: d
       real*8, dimension(d), intent(out) :: tvi
       character (len=*), intent(in) ::  mess
       integer :: ierr, l

       read(ufe,*,iostat=ierr) tvi
       if ( ierr/=0 ) then
          write(6,"('**** erreur (entree) : lecture d''un tableau de reel ****')")
          write(6,"('**** ',a,' ****')") mess
          stop
       end if
       l= len(mess)
       if ( l+d*10<80 ) then
          write(6,"(x,a,' :',100(x,ES9.2))") mess, tvi
       else
          write(6,"(x,a,' :',100(/4x,8(x,ES9.2)))") mess, tvi
       end if
end subroutine entree_tr8_stop

!***********************************************************************
!
!      lecture d'une ch. de caract.,
!      avec arret si echec  (message "mess")
!
!***********************************************************************

subroutine entree_ch_stop(a,mess)
       use entree_data
       implicit none
       character (len=*), intent(out) :: a
       character (len=*), intent(in) ::  mess
       integer :: ierr
       
       read(ufe,"(a)",iostat=ierr) a
       if ( ierr/=0 .or. trim(a)=="" ) then
          write(6,"('**** erreur (entree) : lecture d''une chaine de caracter ****')")
          write(6,"('**** ',a,' ****')") mess
          stop
       end if
       a=adjustl(a)
       write(6,"(x,a,' : ',a)") mess, trim(a)
end subroutine entree_ch_stop

      
!***********************************************************************
!
!      lecture d'une ch. de caract.,
!      avec arret si echec  (message "mess")
!
!***********************************************************************

subroutine entree_repl_ch_clef(clei,a,mess)
       use entree_data
       implicit none
       character(len=*), intent(in) :: clei
       character (len=*), intent(out) :: a
       character (len=*), intent(in) ::  mess
       integer :: ierr
       logical :: OK
       
       call entree_clef(clei,OK)
       if (OK) then
          read(ufe,"(a)",iostat=ierr) a
          if ( ierr/=0 .or. trim(a)=="" ) then
             write(6,"('**** erreur (entree) : lecture d''une chaine de caracter ****')")
             write(6,"('**** ',a,' ****')") mess
             stop
          end if
          a=adjustl(a)
          write(6,"(x,a,' : ',a)") mess, trim(a)
       else
          write(6,"(x,a,' : ',a)") mess, trim(a)
       end if
end subroutine entree_repl_ch_clef

       

 






