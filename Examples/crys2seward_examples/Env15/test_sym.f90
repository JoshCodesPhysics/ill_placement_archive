  Program test_sym
    use ligne_data
    use interface_lect_ligne
    use String_Utility

    implicit none
    real*8, parameter ::  lim1=1.d-4, lim2=1.d-7, lim3=1.d-6

!!$ liste at
       integer, parameter :: nnat=103
       character(len=2), dimension(nnat), parameter :: lat=(/&
            &'h ','he',&
            &'li','be','b ','c ','n ','o ','f ','ne',&
            &'na','mg','al','si','p ','s ','cl','ar',&
            &'k ','ca','sc','ti','v ','cr','mn','fe',&
            'co','ni','cu','zn','ga','ge','as','se','br','kr',&
            &'rb','sr','y ','zr','nb','mo','tc','ru',&
            'rh','pd','ag','cd','in','sn','sb','te','i ','xe',&
            &'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb',&
            &'lu','hf','ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po','at','rn',&
            &'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','es','fm','md','no','lr'/)

!!$ Variable de sortie
       character(len=6) :: sortie="molcas" 

!!$ déclaration des fichiers
       character(len=100) :: nfi_in='coord.in', nfi_sys='sys.coord', nfi_ps='pseudo.coord'
       character(len=100) :: nfi_env='env.coord',  nfi_bid='tmp', prefix=''
       character(len=100) :: nfi_sew='sew0.coord', nfi_xyz='orca.xyz',  nfi_pc='orca.pc'
       character(len=100) :: nfi_fra='frag.molden', nfi_fra2='pseud.molden'
       character(len=100) :: nfi_out='out'
       integer, parameter :: fi_in=1, fi_sys=2, fi_ps=3, fi_env=4, fi_bid=7, &
            fi_sew=8, fi_xyz=9, fi_pc=10, &
            fi_fra=50, fi_fra2=51, fi_out=20

!!$ Operations de symetrie
       character*2, dimension(8) :: oper_sym
       Integer, dimension(3,8) :: opsym, opgpe
       Integer :: nopsym, nopgpe

!!$ Définitions de tous les atomes 
       integer, parameter :: maxatom=100000   ! nbr max d'atomes
       integer :: natom            ! nbre total d'atomes
       Integer ::  natomirr,nirrsys,nirrps,nirrenv
       Integer ::  natomsys,natomps,natomenv
       character*3, dimension(maxatom) :: atom_env    ! nom 
       real*8, dimension(3,maxatom)    :: coord_env   ! positions réelles 
       logical, dimension(maxatom)     :: trouv_env
       real*8, dimension(maxatom)      :: q_env       ! charge 
       character*3, dimension(:), allocatable :: atom_test    ! nom 
       real*8, dimension(:,:), allocatable    :: coord_test   ! positions réelles 
       real*8, dimension(:), allocatable      :: q_test
       logical, dimension(:), allocatable     :: trouv_test


!!$ variables muettes
       integer :: i,j,k, ierr, iatom, isym, ntmp,n
       logical :: OK, succes=.true.
       character(len=40) :: cle
       real*8, dimension(3)    :: coord_tmp
       real*8 :: dist, distmin, q_tmp


!!$
!!$
!!$***********************************************************************  
!!$
       namelist/testenvin/prefix, sortie, nfi_out
!!$
!!$
!!$***********************************************************************  
!!$  Initialisations
       read(5,testenvin)
       if ( trim(prefix)/='' ) then
          prefix=adjustl(prefix)
          if (nfi_in=='        ') nfi_in=trim(prefix)//'.cell'
          nfi_sys=trim(prefix)//'.sys'
          nfi_ps=trim(prefix)//'.psd'
          nfi_env=trim(prefix)//'.env'
          nfi_bid=trim(prefix)//'.tmp'
          nfi_sew=trim(prefix)//'.sew0'
          nfi_xyz=trim(prefix)//'.xyz'
          nfi_pc=trim(prefix)//'.pc'
          nfi_fra=trim(prefix)//'.sys.molden'
          nfi_fra2=trim(prefix)//'.psd.molden'
       end if
       sortie=trim(sortie)
       sortie=StrLowCase(sortie)

       if (sortie.ne.'molcas') then
          write(6,'(//," Pas des symetrie dans ORCA ",//)')
          stop
       else
          open(fi_sew,file=nfi_sew,form="FORMATTED",status="OLD")
       end if

!!$
!!$
!!$***********************************************************************  
!!$  operations de symetrie

       oper_sym(:)='00'
       write(6,9001) " oper_sym : ", oper_sym
       cle="SYMMETRY"
       OK = .false.
       do
          ligne= ""
          read (fi_sew,"(a)",iostat=ierr ) ligne
          if ( ierr < 0 ) stop "erreur lecture sew0"
          call ligne_clef(cle,OK)
          if (OK)  exit
       end do

       write(6,*) " Trouve SYMMETRY"
       read(fi_sew,*) oper_sym
       write(6,'(" oper_sym : ",8(2x,a2))') oper_sym
       
       opsym(:,:) = 0 
       nopsym     = 0
       do i = 1,8
          if (oper_sym(i) .eq. 'x') then 
             opsym(1:3,i) =(/-1,1,1/)
             nopsym = nopsym + 1
          else if (oper_sym(i) .eq. 'y') then 
             opsym(1:3,i) =(/1,-1,1/)
             nopsym = nopsym + 1
          else if (oper_sym(i) .eq. 'z') then 
             opsym(1:3,i) =(/1,1,-1/)
             nopsym = nopsym + 1
          else if (oper_sym(i) .eq. 'xy') then 
             opsym(1:3,i) =(/-1,-1,1/)
             nopsym = nopsym + 1
          else if (oper_sym(i) .eq. 'xz') then 
             opsym(1:3,i) =(/-1,1,-1/)
             nopsym = nopsym + 1
          else if (oper_sym(i) .eq. 'yz') then 
             opsym(1:3,i) =(/1,-1,-1/)
             nopsym = nopsym + 1
          else if (oper_sym(i) .eq. 'xyz') then 
             opsym(1:3,i) =(/-1,-1,-1/)
             nopsym = nopsym + 1
          else 
             exit
          end if
       end do

       write(6,'("Nbre op. sym : ", i4)') nopsym
       write(6,'(" operations de symetrie : ",8(2x,a2))') (oper_sym(i), i=1,nopsym)
       write (6,*)

      ! reconstruire toutes les operations du gpe et pas seulement les irreductibles
      nopgpe =  nopsym
      opgpe(1:3,1:nopsym) = opsym(1:3,1:nopsym)
      do i = 1,nopsym
         do j = 1, nopsym
            opgpe(:,nopgpe+1) = opsym(:,i)*opsym(:,j)
            do k = 1,nopgpe
               if ( opgpe(1,nopgpe+1).eq.opgpe(1,k) .and.  &
                    opgpe(2,nopgpe+1).eq.opgpe(2,k) .and. &
                    opgpe(3,nopgpe+1).eq.opgpe(3,k) ) goto 111
            end do
            nopgpe = nopgpe + 1
111          continue
         end do
      end do

      write (6,*) "Ordre du groupe=",nopgpe
      do i=1,nopgpe
         write (6,*) opgpe(:,i)
      end do
      write (6,*)
      write(6,*)" ----------------------------------------------------"

      nopsym = nopgpe
      opsym  = opgpe

!!$
!!$
!!$***********************************************************************  
!!$ Fragment

       ! lecture dans sew0
       backspace(fi_sew)
       backspace(fi_sew)
       cle='Fragment :'
       OK = .false.
       do
          ligne= ""
          read (fi_sew,"(a)",iostat=ierr ) ligne
          write(6,*) ligne
          if ( ierr < 0 ) stop "erreur lecture sew0"
          call ligne_clef(cle,OK)
          if (OK)  then 
             write(6,*) " Trouve 'Fragment :' "
             exit
          end if
       end do
       
       iatom = 0
       atom_env(:)="  "
       cle='Pseudos :'
       do
          ligne= ""
          read (fi_sew,"(a)",iostat=ierr ) ligne
          if ( ierr < 0 ) stop "erreur lecture sew0"
          call ligne_clef(cle,OK)
          if (OK) exit
          iatom = iatom + 1
          read(ligne(3:4),*,iostat=ierr) atom_env(iatom)(1:2)
          read(ligne(11:27),*,iostat=ierr) coord_env(1,iatom)
          read(ligne(32:48),*,iostat=ierr) coord_env(2,iatom)
          read(ligne(53:79),*,iostat=ierr) coord_env(3,iatom)
       end do
       nirrsys = iatom
       write(6,'("Nbre d''atomes irreductibles dans le fragment: ", i4)') nirrsys
       do i = 1, nirrsys
          write (6,"(2x,a2,3(4x,f17.12))",iostat=ierr) &
                     atom_env(i)(1:2),  coord_env(:,i)
       end do
       write (6,*)


       ! generation de tous les atomes
       natomsys  = nirrsys
       do iatom = 1, nirrsys
          do isym = 1,nopsym
             coord_tmp(:) = opsym(:,isym) * coord_env(:,iatom)
             ! l'atome est son propres sym. je ne le rajoute pas
             if (       abs(coord_tmp(1) - coord_env(1,iatom)) .lt. lim2 &
                  .and. abs(coord_tmp(2) - coord_env(2,iatom)) .lt. lim2 &
                  .and. abs(coord_tmp(3) - coord_env(3,iatom)) .lt. lim2 ) goto 10
             ! ce symetrique a deja ete compte je ne le rajoute pas
             do i = 1,isym-1
                if (       abs(coord_tmp(1) - opsym(1,i)*coord_env(1,iatom)) .lt. lim2 &
                     .and. abs(coord_tmp(2) - opsym(2,i)*coord_env(2,iatom)) .lt. lim2 & 
                     .and. abs(coord_tmp(3) - opsym(3,i)*coord_env(3,iatom)) .lt. lim2 ) goto 10
             end do
             ! c'est un nouveau symetrique
             natomsys = natomsys + 1
             atom_env(natomsys)(1:2) =  atom_env(iatom)(1:2)
             coord_env(:,natomsys) = coord_tmp(:)
!             write(6,*) " nouveau symetrique "
!             write (6,"(2x,a2,3(4x,f17.12))",iostat=ierr) &
!               atom_env(natomsys)(1:2),  coord_env(:,natomsys)
10           continue
          end do
       end do
       write(6,'("Nbre total d''atomes du fragment generes par symetrie: ", i4)') natomsys
       do i = 1, natomsys
          write (6,"(2x,a2,3(4x,f17.12))",iostat=ierr) &
               atom_env(i)(1:2),  coord_env(:,i)
       end do
       write (6,*)

       ! lecture dans sys
       write(6,*) " Lecture dans sys"
       open(fi_sys,file=nfi_sys,form="FORMATTED",status="OLD")
       allocate (atom_test(natomsys), coord_test(3,natomsys))
       read(fi_sys,*)
       do  iatom = 1,natomsys
          ligne= ""
          read (fi_sys,"(a)",iostat=ierr ) ligne
          write(6,*) ligne(1:50)
          if ( ierr < 0 ) stop "erreur lecture sys"
          read(ligne(2:3),*,iostat=ierr)   atom_test(iatom)(1:2)
          read(ligne(6:15),*,iostat=ierr)  coord_test(1,iatom)
          read(ligne(17:26),*,iostat=ierr) coord_test(2,iatom)
          read(ligne(28:37),*,iostat=ierr) coord_test(3,iatom)
       end do
       read(fi_sys,*)
       read(fi_sys,*)
       read (fi_sys,"(a)",iostat=ierr ) ligne
       if (ligne(1:13).ne."  En symetrie") then 
          write(6,*)  "Erreur nbre d''atomes  du systeme"
          write(6,*)  ligne(1:50)
          write(6,*)
          stop 
       end if
       
       ! comparaison avec les atomes dans sys
       ntmp = 0
       dist = 0.d0
       trouv_env(:) = .false.
       do iatom = 1, natomsys
          distmin=1000.d0
          do i = 1,natomsys
             dist = (coord_env(1,iatom)-coord_test(1,i))*(coord_env(1,iatom)-coord_test(1,i)) &
                  + (coord_env(2,iatom)-coord_test(2,i))*(coord_env(2,iatom)-coord_test(2,i)) &
                  + (coord_env(3,iatom)-coord_test(3,i))*(coord_env(3,iatom)-coord_test(3,i))  
             dist = sqrt(dist)
             if (dist.lt.distmin) distmin = dist
             if (atom_env(iatom)(1:2).eq.atom_test(i)(1:2) .and. dist.lt.lim1) then 
                trouv_env(iatom) = .true.
                ntmp = ntmp+1
                goto 11
             end if
          end do
          write(6,*) "atome non trouve", distmin
          write(6,9006) atom_env(iatom)(1:2), coord_env(:,iatom)
11        continue
       end do
       if (ntmp.ne.natomsys) then
          write(6,*) " Tous les atomes n''ont pas ete trouves"
          write(6,*) " Nbre d'atomes trouves: ",ntmp 
       end if
       write(6,*) " Atomes non trouves "
       j = 0
       do i=1,natomsys
          if (.not.trouv_env(i)) then 
             write(9006,*) atom_env(i), coord_env(:,i)
             j = j+1
             succes=.false.
          end if
       end do
       if (j.eq.0)  write(6,*) " Tous les atomes du systeme ont ete trouves "
       
       close(fi_sys)
       write(6,*)
       write(6,*)" ----------------------------------------------------"

!!$
!!$
!!$***********************************************************************  
!!$ Pseudos

       if (OK)  write(6,*) " Trouve 'Pseudos :' "

       ! lecture dans sew0
       iatom = 0
       cle='XFIEld'
       iatom = natomsys
       do
          ligne= ""
          read (fi_sew,"(a)",iostat=ierr ) ligne
          if ( ierr < 0 ) stop "erreur lecture sew0"
          call ligne_clef(cle,OK)
          if (OK) exit
          iatom = iatom + 1
          read(ligne(3:4),*,iostat=ierr) atom_env(iatom)(1:2)
          read(ligne(11:27),*,iostat=ierr) coord_env(1,iatom)
          read(ligne(32:48),*,iostat=ierr) coord_env(2,iatom)
          read(ligne(53:79),*,iostat=ierr) coord_env(3,iatom)
       end do
       natom    = iatom
       nirrps   = natom - natomsys
       write(6,'("Nbre d''atomes irreductibles dans les pseudos : ", i4)') nirrps
       do i = natomsys+1, natomsys+nirrps
          write (6,"(2x,a2,3(4x,f17.12))",iostat=ierr) &
                     atom_env(i)(1:2),  coord_env(:,i)
       end do
       write (6,*)

       ! generation de tous les atomes
       natomps  = nirrps
       do iatom = natomsys+1, natomsys+nirrps
          do isym = 1,nopsym
             coord_tmp(:) = opsym(:,isym) * coord_env(:,iatom)
             if (       abs(coord_tmp(1) - coord_env(1,iatom)) .lt. lim2 &
                  .and. abs(coord_tmp(2) - coord_env(2,iatom)) .lt. lim2 &
                  .and. abs(coord_tmp(3) - coord_env(3,iatom)) .lt. lim2 ) goto 20
             do i = 1,isym-1
                if (       abs(coord_tmp(1) - opsym(1,i)*coord_env(1,iatom)) .lt. lim2 &
                     .and. abs(coord_tmp(2) - opsym(2,i)*coord_env(2,iatom)) .lt. lim2 & 
                     .and. abs(coord_tmp(3) - opsym(3,i)*coord_env(3,iatom)) .lt. lim2 ) goto 20
             end do
             natomps = natomps + 1
             natom   = natom   + 1
             atom_env(natom)(1:2) =  atom_env(iatom)(1:2)
             coord_env(:,natom) = coord_tmp(:)
20           continue
          end do
       end do
       write(6,'("Nbre total d''atomes dans les pseudos generes par symetrie: ", i4)') natomps
       do i = natomsys+1, natomsys + natomps
          write (6,"(2x,a2,3(4x,f17.12))",iostat=ierr) &
               atom_env(i)(1:2),  coord_env(:,i)
       end do
       write (6,*)

       ! lecture dans pseudos
       write(6,*) " Lecture dans pseudos"
       open(fi_ps,file=nfi_ps,form="FORMATTED",status="OLD")
       deallocate(atom_test,coord_test)
       allocate (atom_test(natomps), coord_test(3,natomps))
       read(fi_ps,*)
       do  iatom = 1,natomps
          ligne= ""
          read (fi_ps,"(a)",iostat=ierr ) ligne
          write(6,*) ligne(1:50)
          if ( ierr < 0 ) stop "erreur lecture pseudos"
          read(ligne(2:3),*,iostat=ierr)   atom_test(iatom)(1:2)
          read(ligne(6:15),*,iostat=ierr)  coord_test(1,iatom)
          read(ligne(17:26),*,iostat=ierr) coord_test(2,iatom)
          read(ligne(28:37),*,iostat=ierr) coord_test(3,iatom)
       end do
       read(fi_ps,*)
       read(fi_ps,*)
       read (fi_ps,"(a)",iostat=ierr ) ligne
       if (ligne(1:13).ne."  En symetrie") then 
          write(6,*)  "Erreur nbre d''atomes de pseudos"
          write(6,*)  ligne(1:50)
          write(6,*)
          stop 
       end if
       
       ! comparaison avec les atomes dans ps
       ntmp = 0
       dist = 0.d0
       do iatom = natomsys+1, natomsys+natomps
          distmin=1000.d0
          do i = 1,natomps
             dist = (coord_env(1,iatom)-coord_test(1,i))*(coord_env(1,iatom)-coord_test(1,i)) &
                  + (coord_env(2,iatom)-coord_test(2,i))*(coord_env(2,iatom)-coord_test(2,i)) &
                  + (coord_env(3,iatom)-coord_test(3,i))*(coord_env(3,iatom)-coord_test(3,i))  
             dist = sqrt(dist)
             if (dist.lt.distmin) distmin = dist
             if (atom_env(iatom)(1:2).eq.atom_test(i)(1:2) .and. dist.lt.lim1) then 
                trouv_env(iatom) = .true.
                ntmp = ntmp+1
                goto 21
             end if
          end do
          write(6,*) "atome non trouve", distmin
          write(6,9006) atom_env(iatom)(1:2), coord_env(:,iatom)
21        continue
       end do
       if (ntmp.ne.natomps) then
          write(6,*) " Tous les atomes n''ont pas ete trouves"
          write(6,*) " Nbre d'atomes trouves: ",ntmp 
       end if
       write(6,*) " Atomes non trouves "
       j=0
       do i=natomsys+1, natomsys+natomps
          if (.not.trouv_env(i)) then 
             write(9006,*) atom_env(i), coord_env(:,i)
             j = j+1
             succes=.false.
          end if
       end do
       if (j.eq.0)  write(6,*) " Tous les atomes de pseudos ont ete trouves "
       write(6,*)
       write(6,*)" ----------------------------------------------------"
       close(fi_ps)
       write(6,*)


!!$ 
!!$
!!$***********************************************************************  
!!$ Env
       if (OK)  write(6,*) " Trouve 'XFIEld' "
       read (fi_sew,*) nirrenv

       ! lecture dans sew0
       iatom = natom
       do
          ligne= ""
          read (fi_sew,"(a)",iostat=ierr ) ligne
          if ( ierr < 0 ) exit
          iatom = iatom + 1
          read(ligne(3:19),*,iostat=ierr)  coord_env(1,iatom)
          read(ligne(24:40),*,iostat=ierr) coord_env(2,iatom)
          read(ligne(45:61),*,iostat=ierr) coord_env(3,iatom)
          read(ligne(64:79),*,iostat=ierr) q_env(iatom)
       end do
       natom    = iatom
       if (nirrenv.ne.natom-natomsys-natomps) stop " Erreur nbre d''atomes irreductibles dans env"
       write(6,'("Nbre d''atomes irreductibles dans l''environnement : ", i4)') nirrenv

       do i =  natomsys+natomps+1,  natomsys+natomps+5
          write (6,"(2x,4(4x,f10.5))",iostat=ierr)   coord_env(:,i), q_env(i)
       end do
       write (6,*) " ..... "
       do i =  natom-5,  natom
          write (6,"(2x,4(4x,f10.5))",iostat=ierr)   coord_env(:,i), q_env(i)
       end do

       ! generation de tous les atomes
       natomenv  = nirrenv
       do iatom = natomsys+natomps+1, natomsys+natomps + nirrenv
          do isym = 1,nopsym
             coord_tmp(:) = opsym(:,isym) * coord_env(:,iatom)
             if (       abs(coord_tmp(1) - coord_env(1,iatom)) .lt. lim2 &
                  .and. abs(coord_tmp(2) - coord_env(2,iatom)) .lt. lim2 &
                  .and. abs(coord_tmp(3) - coord_env(3,iatom)) .lt. lim2 ) goto 30
             do i = 1,isym-1
                if (       abs(coord_tmp(1) - opsym(1,i)*coord_env(1,iatom)) .lt. lim2 &
                     .and. abs(coord_tmp(2) - opsym(2,i)*coord_env(2,iatom)) .lt. lim2 & 
                     .and. abs(coord_tmp(3) - opsym(3,i)*coord_env(3,iatom)) .lt. lim2 ) goto 30
             end do
             natomenv = natomenv + 1
             natom    = natom   + 1
             coord_env(:,natom) = coord_tmp(:)
             q_env(natom)       = q_env(iatom)
30           continue
          end do
       end do
       write(6,'("Nbre d''atomes total dans l''environnement : ", i9)') natomenv
       write(6,'("Nbre d''atomes total                       : ", i9)') natom
       write (6,*)

       ! Lecture du nombre total d''atomes dans out
       write(6,*) " Lecture du nombre total d''atomes dans out"
       open(fi_out,file=nfi_out,form="FORMATTED",status="OLD")
       cle="Nbre total d'atomes"
       OK = .false.
       do
          ligne= ""
          read (fi_out,"(a)",iostat=ierr ) ligne
          if ( ierr < 0 ) stop "erreur lecture out"
          call ligne_clef(cle,OK)
          if (OK)  exit
       end do
       backspace(fi_out)
       read(fi_out,"(22x,i13)") ntmp
       write(6,*) " Nbre total d''atomes lu dans out : ", ntmp
       if (ntmp.ne.natom) write(6,*) " >>>> ERREUR dans la symetrie <<<< "
       write(6,*)

       ! lecture dans env 
       write(6,*) " Lecture dans env"
       open(fi_env,file=nfi_env,form="FORMATTED",status="OLD")
       deallocate(atom_test,coord_test)
       n = ntmp - (natomsys+natomps)
       allocate (coord_test(3,n), q_test(n), trouv_test(n))
       read(fi_env,*)
       do  iatom = 1,n
          ligne= ""
          read (fi_env,"(a)",iostat=ierr ) ligne
!          write(6,*) ligne(1:50)
          if ( ierr < 0 ) stop "erreur lecture env"
          read(ligne(4:18),*,iostat=ierr)  coord_test(1,iatom)
          read(ligne(22:36),*,iostat=ierr) coord_test(2,iatom)
          read(ligne(40:54),*,iostat=ierr) coord_test(3,iatom)
          read(ligne(58:72),*,iostat=ierr) q_test(iatom)
       end do
       read(fi_env,*)
       read(fi_env,*)
       read (fi_env,"(a)",iostat=ierr ) ligne
       OK=.true.
       if (ligne(1:13).ne."  En symetrie") then
          OK=.false.
          write(6,*)  "Erreur nbre d''atomes de l''environnement"
          write(6,*)  ligne(1:50)
          write(6,*)
          stop 
       end if

       ! comparaison avec les atomes dans env
       ntmp = 0
       dist = 0.d0
       trouv_test(:)=.false.
       do iatom = natomsys+natomps+1, natom
          distmin=1000.d0
          do i = 1,n
             dist = (coord_env(1,iatom)-coord_test(1,i))*(coord_env(1,iatom)-coord_test(1,i)) &
                  + (coord_env(2,iatom)-coord_test(2,i))*(coord_env(2,iatom)-coord_test(2,i)) &
                  + (coord_env(3,iatom)-coord_test(3,i))*(coord_env(3,iatom)-coord_test(3,i))  
             dist = sqrt(dist)
             if (dist.lt.distmin) distmin = dist
             if (dist.lt.lim1.and. abs(q_env(iatom)-q_test(i)).lt.lim3) then 
                if (trouv_env(iatom)) write(6,*) " charge generee par symetrie ", &
                     iatom-natomsys+natomps," déjà trouvee"
                if (trouv_test(i))  then 
                   write(6,*) " charge de env deja generee par symetrie"
                   write(6,9007) coord_test(:,i), q_test(i)
                end if
                trouv_env(iatom) = .true.
                trouv_test(i)    = .true.
                ntmp = ntmp+1
                goto 31
             end if
          end do
          write(6,*) "atome non trouve", distmin
          write(6,9007) coord_env(:,iatom), q_env(iatom)
31        continue
       end do

       if (ntmp.ne.natomenv .or. ntmp.ne.n) then
          write(6,*) " >>>> Tous les atomes n''ont pas ete trouves <<<<"
          write(6,*) " Nbre d'atomes trouves: ",ntmp 
          succes=.false.
       end if
       write(6,*)
       write(6,*) " Atomes generes par symetrie et non trouves dans env"
       i=0
       do iatom = natomsys+natomps+1, natom
          if (.not.trouv_env(iatom)) then 
             write(6,9007) coord_env(:,iatom), q_env(iatom)
             i = i+1
             succes=.false.
          end if
       end do
       write(6,*) " Nbre ", i
       write(6,*)
       i=0
       write(6,*) " Atomes trouves dans env et non generes par symetrie"
       do iatom = 1,n
          if (.not.trouv_test(iatom)) then
             write(6,9007)  coord_test(:,iatom), q_test(iatom)
             i=i+1
             succes=.false.
          end if
       end do
       write(6,*) " Nbre ", i
       write(6,*)
       

       close(fi_env)
       write(6,*)



!!$ 
!!$
!!$***********************************************************************  
!!$  Fin

       if (succes) then 
          write(6,*) " La symetrie semble OK"
       else
          write(6,*) ">>>> pb de gestion de la symetrie <<<"
       end if
       write (6,*) 
       write (6,*) " fin"



9001   format(8(2x,a2))
9005   format (1x,a3,1x,6(f10.6,1x),1x,'Angstrom',2x,10(f10.6,x))
9006   format (1x,a3,1x,3(f10.6,1x))
9007   format (1x,4(f10.6,1x))
  end Program test_sym
