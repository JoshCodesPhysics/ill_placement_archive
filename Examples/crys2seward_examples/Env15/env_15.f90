!!$***********************************************************************
!!$
!!$                Programme de calcul des environnements  
!!$        M�thode d'annulation des moments dipolaires d'Alain Gell�
!!$           Marie-Bernadette Lepetit, Alain Gell� - Octobre 2004
!!$
!!$
!!$
!!$ Donn�es pour namelist :
!!$
!!$    prefix : pr�fixe pour les fichiers de sorties ( prefix.sys, prefix.psd, prefix.env, prefix.tmp ) 
!!$             ( sans pr�fixe ce sont : sys.coord, pseudo.coord, env.coord )
!!$    nfi_in : fichiers avec les atomes de la cellule �l�mentaire (d�faut : coord.in ) 
!!$    sortie : molcas (default) g�n�re le fichier     prefix.sew0
!!$           : orca,            g�n�re les fichiers   prefix.xyz (fragment+pseudo)
!!$                                                 et prefix.pc  (env)
!!$
!!$    SINGLE : true  : pour fixer le nombre de moment a annuler
!!$             false : le nombre de moment augmente jusqu'a convergence du pot  (d�faut)
!!$    prec   : precision pour la convergence du potentiel  (defaut : 10-1 meV )
!!$             ( le potentiel est calcul� au centre et sur les atomes magn�tiques)
!!$    maxatom : nombre maximal d'atomes dans l'envirr.
!!$              (d�faut= 100 000)
!!$
!!$    ndip   : nbre de moments dipolaires � annuler (d�faut = 1)    
!!$             (inutile si SINGLE=.false.)
!!$    mindip/maxdip : valeur min et max pour la boucle sur ndip
!!$                    (d�faut : mindip=1, maxdip=10)
!!$    adipm   : nombre de moment dipolaire affich�s, en plus des ndip annul�s
!!$              (defaut : adipm=0)
!!$
!!$    nch  : nbre d'atomes dans la cellule �l�mentaire  (d�faut = 1)
!!$    a,b,c, alpha,beta,gamma : param�tres de la maille �l�mentaire 
!!$           (d�fault = 1,1,1,pi/2,pi/2,pi/2)
!!$
!!$    a2,b2,c2 :  axes de references pour creer la base orthonormee finale
!!$                X= a2  ;  Y= b2-(X.b2)X  ;  Z= c2-(X.c2)X-(Y.c2)Y 
!!$                +normalisation
!!$                ils sont exprim�s dans la base a,b,c
!!$                (d�faut a2=(1,0,0), b2=(0,1,0), c2=(0,0,1) )
!!$
!!$    ncel : nbre de cellules par dimension sans ajustement de charge 
!!$           (en sus de la cellule centrale, d�fault = 0) 
!!$    pos0 : sp�cifications de l'atome central (d�fault : pos = O,O,O)
!!$
!!$    nmag : nbre d'atomes magn�tiques  (d�fault = 0) 
!!$    atommag, posmag : sp�cifications des atomes magn�tiques
!!$
!!$    ydebug : debug (d�fault = F)
!!$    iprint : niveau d'impression de sortie (default = 0)
!!$    ysym : T (default) calcul des op�rations de sym�trie en d2h maximum
!!$    lim1 : 1.d-6  par d�fault, seuil d'�galit� de charges
!!$    lim2 : 1.d-5 par d�fault, limite d'erreur sur la position 
!!$           pour que deux atomes soient consid�r�s comme sym�triques 
!!$    lim3 : 1.d-9  par d�fault, seuil d'�limination des atomes de charge nulle
!!$    lim4 : 1.d-6  par d�fault, seuil de nullit� de position pour un atome sur un axe de sym�trie
!!$    lim5 : 1.d-9  par d�fault, seuil de nullit� de la distance pour le pot �lec. 
!!$  
!!$***********************************************************************
!!$***********************************************************************
!!$
!!$ 10/10/06  version _1 : qq corrections
!!$ 11/10/06  version _2 : entrees par mots clef au choix
!!$                        boucle sur ndip -> precision sur pot
!!$ 17/10/06  version _3 : rotation possible des axes de la base d'arriv�
!!$                        symetrie : on ne prend que les plqns de sym (pas les inversions)
!!$ 23/11/06  version _4 : conv de la diff de pot aussi sur les at. du syst.
!!$
!!$ 07/12/06  version _5 : dist pour les pseudos : utilise les atomes du systeme
!!$
!!$ 11/11/07  version _8 : supprimer des sym pour sebastien
!!$
!!$ 27/05/08  version _10 : atomes magn en A ou bohr (coord. reelles)
!!$                         et on peut donner une nouvelle base 
!!$                         avec axes en coord. reelles A ou B (ancienne base)
!!$
!!$ 27/06/12 version _12 : version namelist corrig�e + coord_rot_roti corrig�
!!$----- D�clarations  
!!$
!!$ 10/03/14 version _14 : on n'enl�ve pas les charges nulles dans sys et pseud
!!$                        adaptation des sorties aux grands syst�mes
!!$                        sortie orca
!!$
!!$ 1/10/15 version _15  : on traite la sym�trie d'odre 3 en faisant une moyenne de 3 environnements 
!!$                      : on corrige le choix des atomes de sym�trie qui ne prennaient 
!!$                        pas les atomes sur la diagonale d'un cube
!!$
Program Env 
       use interface_lect_ligne
       use interface_lect_entree
       use String_Utility
       implicit none
!!$ type de lecture d'entree : namelist / mot-clef
       logical, parameter :: entreeTypeClef=.true.
!!$ Variables globales
       logical :: ydebug=.false. , ysym=.true.
       real*8,  parameter :: pi=3.14159265358979d0, ua=0.52917720859d0, &
            mev=27211.39628d0, zero=0.0d0
       integer, parameter :: maxnmag=50   !nbre max d'atomes magn�tiques
       integer :: iprint=0 
       real*8 ::  lim1=1.d-6, lim2=1.d-5, lim3=1.d-18, lim4=1.d-6, lim5=1.d-9

!!$ Variable de sortie
       character(len=6) :: sortie="molcas" 

!!$ d�claration des fichiers
       character(len=100) :: nfi_in='coord.in', nfi_sys='sys.coord', nfi_ps='pseudo.coord'
       character(len=100) :: nfi_env='env.coord',  nfi_bid='tmp', prefix=''
       character(len=100) :: nfi_sew='sew0.coord', nfi_xyz='orca.xyz',  nfi_pc='orca.pc'
       character(len=100) :: nfi_fra='frag.molden', nfi_fra2='pseud.molden'
       integer, parameter   :: fi_in=1, fi_sys=2, fi_ps=3, fi_env=4, fi_bid=7, &
            fi_sew=8, fi_xyz=9, fi_pc=10, &
            fi_fra=50, fi_fra2=51
       character(len=50), parameter :: fich_tmp="fich_env_entree.tmp"

!!$ liste at
       integer, parameter :: nnat=103
       character(len=2), dimension(nnat), parameter :: lat=(/&
            &'h ','he',&
            &'li','be','b ','c ','n ','o ','f ','ne',&
            &'na','mg','al','si','p ','s ','cl','ar',&
            &'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu','zn','ga','ge','as','se','br','kr',&
            &'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd','in','sn','sb','te','i ','xe',&
            &'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb',&
            &'lu','hf','ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po','at','rn',&
            &'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','es','fm','md','no','lr'/)

!!$ D�finition des atomes de la cellule �l�mentaire 
       real*8, dimension(:,:), allocatable :: pos     !positions fractionnaires 
       real*8, dimension(:,:), allocatable :: coord   !positions r�elles
       real*8, dimension(:), allocatable :: q         !charge 
       character*3, dimension(:), allocatable :: atom !nom 
       integer, dimension(:), allocatable :: az   !num at
       real*8, dimension(3) :: coord_1, pos_1

!!$ D�finitions de tous les atomes 
       integer :: maxatom=100000   ! nbr max d'atomes
       integer :: ncel=0           ! nbr de cellules exactes mis � part la centrale
       integer :: natom            ! nbre total d'atomes
       real*8, dimension(:,:), allocatable :: pos_env     ! positions fractionnaires 
       real*8, dimension(:,:), allocatable :: coord_env   ! positions r�elles 
       real*8, dimension(:),   allocatable :: q_env       ! charge 
       character*3, dimension(:), allocatable :: atom_env ! nom 
       integer, dimension(:), allocatable :: az_env   ! num at
       integer :: natenv  ! nbr de charge dans env
       integer :: natsys  ! nbr d'at. dans le syst. pour la conv.
       integer :: natfra, natfra2  ! nbr at dans le
       real*8, dimension(:,:), allocatable :: coord_sys
       character*3, dimension(:), allocatable :: atom_sys ! nom 

!!$ variables de renormalisation des charges (n� de l'atome, coordonn�e, cellule)
       real*8, dimension(:,:,:), allocatable :: lambda
       real*8, dimension(:,:), allocatable :: mu
       integer, dimension(3) :: pmin, pmax

!!$ variables pour le moment dip
       logical :: SINGLE=.false.      ! une seule val. de ndip.
       real*8 :: prec=1.E-1           ! precision pour la diff de pot
       real*8 :: prec2=1.E-4           ! precision pour le pot
       integer :: ndip=1              ! nbr de moments annul�s (p-1 = ndip)
       integer :: adip=1,  adipm=0  ! pour les moments affich�s : de 1 a ndip+adipm
       integer :: mindip=1, maxdip=10 ! valeures min et max de ndip

!!$ variables de d�finition du cristal
       integer :: nch=1                 ! nbr de charges de la maille �l�mentaire
       real*8 :: a=1.0, b=1.0, c=1.0, & ! maille du cristal   
            alpha=90.0, beta=90.0, gamma=90.0
       character*3 :: atom0='x'         ! atome central
!!$ base
       real*8, dimension(3) :: a2=(/1.D0,0.D0,0.D0/), b2=(/0.D0,1.D0,0.D0/), c2=(/0.D0,0.D0,1.D0/)
       real*8, dimension(3,3) :: R, RI
       real*8, dimension(3) :: pos0=(/ 0.0,0.0,0.0 /)! coord. frac. atome central
       logical :: NVBASE=.false.
       real*8, dimension(3,3) :: nvB
       character, dimension(3) :: tpnvb
!!$ frag
       integer :: nmag=0                ! nbre d'atomes magn�tiques         
       character*3, dimension(maxnmag)   :: atommag ! atomes magn�tiques
       real*8,      dimension(3,maxnmag) :: posmag,coordmag  ! coord. frac. atomes mag.
       character, dimension(maxnmag)   :: tpmag ! type de position (rien : fractionnaire, A ou B(ohr) )

!!$ variables de partition syst�me, pseudos, environnement 
       real*8 :: dsys, dpseud
       real*8, dimension(:), allocatable :: dist
       real*8 :: dist_min, dist_1, dist_min_sys

!!$ variables de traitement de la sym�trie
       Integer :: maxsym, nopsym, natomirr,nirrsys,nirrps,nirrenv
       Integer, dimension(:), allocatable :: sym, iw, typat
       Integer, dimension(:,:), allocatable :: isym
       Integer, dimension(3,8) :: opsym
       Real*8, dimension(3) :: coordsym
       logical, dimension(3) :: psym, asym
       logical :: ptsym, oksym, yaxe3=.false.
       !--- pour en virer
       logical :: DELSYM=.false.
       integer :: ndelsym=0
       Integer, dimension(3,8) :: opdelsym


!!$ Potentiel electrostatique
       real*8, dimension(:,:), allocatable :: pot 
       real*8, dimension(:,:), allocatable :: dpot 
       real*8, dimension(:,:), allocatable :: vpot, dvpot 
       real*8 :: minpot, minvpot, maxsysdpot, maxsysdvpot

!!$ variables muettes
       character*5 :: bla
       character*8 :: blabla
       character*3 :: aat, aat1
       integer :: i,j,k,l,m,iatom,idim,p,nn,is,ii,jj,kk,jmin,jmax,ierr, ifra
       real*8  :: mpol, tmp, tpol(3), mpolm, vmax, vmin
       logical :: OK, OK2
       
       integer :: natom_axis3
       character*3, dimension(:), allocatable :: atom_axe3! nom des atomes
       integer, dimension(:), allocatable :: az_axe3   ! numeros atomiques 
       real*8, dimension(:,:), allocatable :: pos_axe3 !positions fractionnaires
       real*8, dimension(:), allocatable :: q_axe3 ! charge renormalis�e

!!$
!!$
!!$***********************************************************************  
!!$
       namelist/envin/prefix, nfi_in, sortie, &
            SINGLE,ndip,mindip,maxdip,adipm, prec, &
            ncel,nch,maxatom, &
            a,b,c, alpha,beta,gamma, &
            NVBASE, a2,b2,c2, nvb, tpnvb,&
            atom0,pos0, nmag,atommag,posmag, &
            dsys, dpseud, &
            ydebug, iprint, &
            ysym, delsym, ndelsym, opdelsym, yaxe3, &
            lim1, lim2, lim3, lim4, lim5 
!!$
!!$
!!$***********************************************************************  
!!$
       write(6,"(/x,70('=')/)")
       write(6,*)'                Programme Env '
       write(6,*) 
       write(6,*)' Ecrit par Alain Gell� et Marie-Bernadette Lepetit'
       write(6,*)' D''apr�s la m�thode d''Alain Gell�' 
       write(6,*) 
       write(6,*)' Merci de citer les r�ferences suivantes~: '  
       write(6,*)' Alain gell�, th�se de doctorat, Toulouse, d�cembre 2004'
       write(6,*)'  (article en pr�paration)' 
       write(6,*) 
       write(6,*)' Vous �tes autoris� � utiliser ce programme '
       write(6,*)' et de le diffuser en l''�tat sous r�serve d''en informer '
       write(6,*)' MB Lepetit ou A Gell�'
       write(6,*) 
       write(6,*)' Vous �tes autoris� � modifier/ am�liorer ce programme '
       write(6,*)' sous r�serve de transmettre vos modifications � '
       write(6,*)' MB Lepetit ou A Gell� de mani�re � ce que '
       write(6,*)' vos am�liorations profitent � tous. '
       write(6,*)
       write(6,*)' Pour �tre sur la liste de diffusion des nouvelles versions '
       write(6,*)' envoyez un mail � MB Lepetit'
       write(6,*)
       write(6,*)' MB Lepetit '
       write(6,*)' Institut N�el, CNRS UPR~2940'                       
       write(6,*)' 25 rue des Martyrs, BP 166, B�t. D'
       write(6,*)' 38042 Grenoble cedex 9'
       write(6,*)' FRANCE'
       write(6,*)' Courriel: Marie-Bernadette.Lepetit@Grenoble.CNRS.fr'
       write(6,*)
       write(6,*)' Alain Gell� '
       write(6,*)' Institut de Physique de Rennes '
       write(6,*)' B�t 11 A'
       write(6,*)' 263 Avenue du General Leclerc CS 74205'
       write(6,*)' 35042 Rennes Cedex'
       write(6,*)' FRANCE'
       write(6,*)' Courriel : alain.gelle@univ-rennes1.fr'
       write(6,*)
       !
       !***********************************************************************  
       !
       !----- Lecture des donnees
       !
       write(6,"(/x,70('='))")
       call entree_init(1,fich_tmp)
       call entree_clef("clef",OK2)
       !
       !----- lecture de la namelist
       !       
       if ( .not. OK2 ) then
          write(6,"(x,'entrees par namelist ')")
          blabla='coord.in'
          nfi_in='        '
          rewind(1)
          read(1,envin)
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
          if (sortie.eq.'orca') ysym=.false.
          if ( SINGLE ) then
             mindip= ndip
             maxdip= ndip
             write(6,"(' Mode single : ')") 
             write(6,"('   ndip=',I3)") ndip
          else 
             write(6,"(' Mode optim : ')") 
             write(6,"('   Prec =', f8.4,'   Mindip=',I3,'   Maxdip=',I3)") prec, mindip, maxdip
          end if
          write(6,"(' Parametres de maille : ',3(f8.4,2x))") a, b, c
          write(6,"(' Angles maille        : ',3(f8.2,2x))") alpha, beta, gamma
          if (NVBASE) then
             do i = 1,3
                tpnvb(i)= adjustl(tpnvb(i))
                call lett_min(tpnvb(i))
             end do
          end if
          tpmag=" "
          write(6,"(' Nbr de sym a supprimer =',i3)") ndelsym
          do i = 1, ndelsym
             write(6,"('   Operation de sym a supprimer =',3(i3,2x))") opdelsym(:,i)
          end do
          if (ysym.eqv..false.)  write(6,"(' Pas de symetries')")
          write(6,"(' Limites =',5(ES8.1,2x))") lim1, lim2, lim3, lim4, lim5
          write(6,"(' Niveau impression =',i3)") iprint
          
          !
          !----- lecture par mots clef
          !       
       else
          write(6,"(x,'entrees par mots clef ')")
          !--- version debug
          call entree_clef("ydebug",ydebug)
          if ( ydebug )  write(6,"(' **** version debugage')")
          !--- nom de fichiers
          ! prefix
          call entree_clef("pref",OK)
          if (OK) then
             call entree_ch_stop(prefix,"prefix")
             nfi_in=trim(prefix)//'.cell'
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
          !--- Programme a interfacer en sortie
          call entree_clef("sortie",OK)
          if (OK) then
             call entree_ch_stop(sortie,"sortie")
             sortie=StrLowCase(trim(sortie))
             if (sortie.eq.'orca') ysym=.false.
          end if
          !--- fichier de cellules
          call entree_clef("fcell",OK)
          if (OK) call entree_ch_stop(nfi_in,"fichier cell")
          !--- conditions de convergence
          ! single shoot !
          call entree_clef("single",SINGLE)
          if (SINGLE) then
             call entree_int_stop(ndip,"ndip apres single")
             mindip= ndip
             maxdip= ndip
          end if
          ! max atomes
          call entree_clef("maxatom",OK)
          if (OK) read(1,*) maxatom
          ! precision
          call entree_clef("prec",OK)
          if (OK) then
             if ( SINGLE ) then
                write(6,"(/' **** precision inutile en mode single ****')") 
             else
                call entree_r8_stop(prec,"precision")
             end if
          end if
          !--- distance sys
          call entree_clef("dsys",OK)
          if (OK) call  entree_r8_stop(dsys,"dist. du systeme")
          ! distace pseudos
          call entree_clef("dpseud",OK)
          if (OK) call  entree_r8_stop(dpseud,"dist. des pseudos")
          ! mindip : maxdip
          call entree_clef("mindip",OK)
          if (OK) then
             if ( SINGLE ) then
                write(6,"(/' **** mindip inutile en mode single ****')") 
             else
                call entree_int_stop(mindip,"nbr min de moments dip annules")
             end if
          end if
          call entree_clef("maxdip",OK)
          if (OK) then
             if ( SINGLE ) then
                write(6,"(/' **** maxdip inutile en mode single ****')") 
             else
                call entree_int_stop(maxdip,"nbr max de moments dip annules")
             end if
          end if

          !--- donnees sur le crystal
          ! nch
          call entree_clef("nch",OK)
          if (OK) call entree_int_stop(nch,"nbr atomes dans cell")
          if (.not. OK) write(6,"(' nbr atomes dans cell : ',i4)") nch
          ! a, b, c
          call entree_clef("dim",OK)
          if (OK) then
             call entree_ligne_stop()
             call ligne_r8_stop(a,"cell dimension 1")
             call ligne_r8_stop(b,"cell dimension 2")
             call ligne_r8_stop(c,"cell dimension 3")
          end if
          write(6,"(' dimensions cellule : ',3(f8.4,2x))") a, b, c
          ! angles
          call entree_clef("angle",OK)
          if (OK) then
             call entree_ligne_stop()
             call ligne_r8_stop(alpha,"cell angle 1")
             call ligne_r8_stop(beta,"cell angle 2")
             call ligne_r8_stop(gamma,"cell angle 3")
          end if
          write(6,"(' angles cellule : ',3(f8.2,2x))") alpha, beta, gamma
          ! axes de ref
          call entree_clef("base",OK)
          if (OK)  then
             call entree_tr8_stop(a2,3,"base sortie vecteur 1")
             call entree_tr8_stop(b2,3,"base sortie vecteur 2")
             call entree_tr8_stop(c2,3,"base sortie vecteur 3")
          end if
          ! nouvelle base
          call entree_clef("nvbase",NVBASE)
          if (NVBASE) then
             do i= 1, 3
                write(bla,"(i1)") i
                call entree_ligne_stop()
                call ligne_tr8_stop(nvb(:,i),3,"nouvelle base sortie vecteur "//trim(bla))
                call ligne_ch(tpnvb(i),OK)
                call lett_min(tpnvb(i))
             end do
          end if
          ! centre 
          call entree_clef("centre",OK)
          if (OK) call entree_tr8_stop(pos0,3,"centre")
          if (.not. OK) write(6,"(x,a,' :',100(x,ES9.2))") "centre", pos0
          if (.not. OK) write(6,"(' centre : ',3(f8.4,2x))") pos0
          ! ncel
          call entree_clef("ncel",OK)
          if (OK)  call entree_int_stop(ncel,"nbr cell intactes")
          if (.not. OK) write(6,"(x,a,' : ',i4)") "nbr cell intactes", ncel
          ! fragment
          call entree_clef("fragment",OK)
          if (OK)  then
             read(1,*) nmag
!!$          les tabl. suivant ne sont plus dynamique a cause de la nameliste
!!$          allocate ( atommag(nmag) )
!!$          allocate ( posmag(3,maxmag), coordmag(3,maxmag) )
             do i=1, nmag
                call entree_ligne_stop()
                call ligne_ch_stop(atommag(i),"nom atome magn")
                call ligne_tr8_stop(posmag(:,i),3,"pos atome magn")
                call ligne_ch(tpmag(i),OK)
                call lett_min(tpmag(i))
             end do
          else
             write(6,"(/' **** pas de systeme ****')")
          end if
          ! supprimer des sym
          call entree_clef("delsym",DELSYM)
          if (DELSYM) then
             call entree_int_stop(ndelsym,"nbr de sym a supprimer")
             do i= 1, ndelsym
                call entree_tint_stop(opdelsym(:,i),3,"operation de sym a supprimer")
             end do
          end if
          ! pas de sym
          call entree_clef("nosym",OK)
          ysym= .not. OK
          if (OK)  write(6,"(' pas de symetries')")
          ! Axe de symetrie d'ordre 3
          call entree_clef("yaxe3",OK)
          yaxe3=.true.
          !--- donnees techniques
          ! limites
          call entree_clef("limites",OK)
          if (OK) then
             write(6,"(/' ****  changement des limites : ****')")
             write(6,"('   defaut  :',5(ES8.1,2x)/)") lim1, lim2, lim3, lim4, lim5
             call entree_ligne_stop()
             call ligne_r8_stop(lim1,"limite 1")
             call ligne_r8_stop(lim2,"limite 2")
             call ligne_r8_stop(lim3,"limite 3")
             call ligne_r8_stop(lim4,"limite 4")
             call ligne_r8_stop(lim5,"limite 5")
             write(6,"('   entrees :',5(ES8.1,2x)/)") lim1, lim2, lim3, lim4, lim5
          end if
          ! iprint
          call entree_clef("iprint",OK)
          if (OK)  then
             call entree_int_stop(iprint,"niveau impression :")
          end if
          ! iprint
          call entree_clef("affdip",OK)
          if (OK)  then
             call entree_int_stop(adipm,"niveau impression moment dip.")
          end if
          !
          !----- fin
          !
       end if
       call entree_close()

       !
       !----- sortie des entrees
       !

       write(6,"(/x,70('='))")
       write(6,"(x,'=== fichiers '/)")
       write(6,"(' cellule : ',a)") nfi_in
       write(6,"(' systeme : ',a/' pseudos : ',a/' env     : ',a )") nfi_sys, nfi_ps, nfi_env
       write(6,"(' temp    : ',a)") nfi_bid
       write(6,"(/x,70('='))")
       write(6,"(x,'=== methode ')")
       if ( SINGLE ) then
          write(6,"(/' methode choisie :   single')")
          write(6,"(' nbr de moment a annuler fixe a : ',i4)") maxdip
       else
          write(6,"(/' methode choisie :   convergence du potentiel')")
          write(6,"(' precision sur le pot. : ',ES9.1)") prec
          write(6,"(' val min/max pour nbr de moments a annuler : ',i4,2x,i4)") mindip, maxdip
       end if

       !
       !***********************************************************************  
       !
       !----- passage en ua .... euh non ... finallement non ...
       !
       write(6,"(/x,70('='))")
       write(6,"(x,'=== dimensions de la cellule')")
       write(6,"(/' en Ang. : ',3(2x,f12.8))") a, b, c
       write(6,"(' en Bohr : ',3(2x,f12.8))") a/ua, b/ua, c/ua

       !
       !***********************************************************************  
       !
       !----- calcul de la rotation pour obtenir les coord reelles
       !
       write(6,"(/x,70('='))")
       write(6,"(x,'=== rotation de axes')")
       write(6,"(/x,'(a,b,c)    : axes du cristal')")
       write(6,"(x,'(a2,b2,c2) : axes de sortie (non orthonormes)')")
       write(6,"(x,'(X,Y,Z)    : axes de sortie (orthonormes)')")
       call coord_rot_roti(R,RI,a2,b2,c2,a,b,c,alpha,beta,gamma,iprint)

       !
       !***********************************************************************  
       !
       !----- position du centre en position fractionnaire
       !
       write(6,"(/x,70('='))")
       write(6,"(x,'=== position fractionnaire du centre')")
       write(6,"(x,'    Attention si le centre n''est pas au centre de gravite ')")
       write(6,"(x,'      des atomes magnetiques Pb avec la symetrie')")
       write(6,*)

       write(6,"(x,'Centre : ',3(3x,f9.6))") pos0(:)

       !
       !***********************************************************************  
       !
       !----- position des atomes magnetiques en position fractionnaire
       !
       write(6,"(/x,70('='))")
       write(6,"(x,'=== position fractionnaire des atomes magnetiques')")
       write(6,*)

       do i= 1, nmag
          if ( tpmag(i)=="b" ) posmag(:,i)= posmag(:,i)*ua
          if ( tpmag(i)=="a" .or. tpmag(i)=="b" ) then
             pos_1= matmul(RI,posmag(:,i))
             posmag(:,i)= pos_1 + pos0
          end if
          write(6,"(x,i3,3(3x,f9.6))") i, posmag(:,i)
       end do

       !
       !***********************************************************************  
       !
       !----- on utilise eventuellement une autre base
       !

       if (NVBASE) then
          write(6,"(/x,70('='))")
          write(6,"(x,'=== nouveaux axes de sortie')")
          write(6,*)
          !--- axe en pos fractionnaire
          do i= 1, nmag
             if ( tpnvb(i)=="b" )   nvb(:,i)= nvb(:,i)*ua
             if ( tpnvb(i)=="b" .or. tpnvb(i)=="a" ) then
                pos_1= matmul(RI,nvb(:,i))
                nvb(:,i)= pos_1 + pos0
             end if
          end do
          a2= nvb(:,1)
          b2= nvb(:,2)
          c2= nvb(:,3)
          !--- on fait la rotation
          call coord_rot_roti(R,RI,a2,b2,c2,a,b,c,alpha,beta,gamma,iprint)
       else 
          write(6,"(/x,70('='))")
          write(6,"(x,'=== Axes de sortie identiques � axes orthonormalises des axes cristallo')")
          write(6,*)
       end if
       !
       !***********************************************************************  
       !
       !----- Lecture de la maille
       !
       write(6,"(/x,70('='))")
       write(6,"(x,'=== lecture des atomes de la maille')")

       allocate (atom(nch),pos(3,nch),q(nch),coord(3,nch),az(nch))  
       az= 0

       open (unit=fi_in,file=nfi_in,form='formatted',status='old')

       vmin= 0.D0
       vmax= 0.D0
       do i = 1,nch
          read(fi_in,*,iostat=ierr) atom(i), pos(1,i), pos(2,i), pos(3,i), q(i) 
          if (ierr/=0) then
             write(6,"('**** erreur (main) : lecture du fichier de la cellule')")
             write(6,"('**** lecture (nom-position(3)-charge), ligne ',i4,' ****')") i
             stop
          end if
          vmax= max(vmax,maxval(pos(:,i)))
          vmin= min(vmin,minval(pos(:,i)))
          aat=atom(i)
          call lett_min(aat)
          do j= 1, nnat
             if ( index(aat,lat(j))>0 ) then
                az(i)= j
                exit
             end if
          end do
          if (az(i)==0) then
             aat1= aat
             do
                l= len_trim(aat1)
                if ( l<=1 ) exit
                aat1(l:l)=' '
                do j= 1, nnat
                   if ( index(aat1,lat(j))>0 ) then
                      az(i)= j
                      exit
                   end if
                end do
             end do
             if (az(i)==0) then 
                write(6,"('**** atome ',i4,' de la maille (',a3,&
                     &') n''est pas dans le tabl. period. ****')") i, aat
             else
                write(6,"('**** atome ',i4,' (',a3,&
                     &') est suppose etre un ''',a3,''' ****')") i, aat, aat1
             end if
          end if
          if (ydebug) write(6,9001) atom(i), pos(1,i), pos(2,i), pos(3,i), q(i)
       end do

       if (vmin<0.D0 .or. vmax>1.D0) then
          write(6,"('**** attention : ****')")
          write(6,"('**** lecture des atomes de la maille ****')")
          write(6,"('**** des positions sont <0 ou >1  ****')")
          write(6,"('**** les positions sont elles factionnaires ? ****')")
       end if
       close (fi_in)

       !--- si debug
       if (ydebug) then 
          write(6,*) 
          write(6,*) 'ydebug = ', ydebug
          write(6,*) ' Coordonnees des atomes tels que lus',&
               ' dans coord.in'
           write(6,*) 'Atome    Positions fractionnaires      Coordonnees reelles'   
          coord= matmul(R,pos)
          do i=1,nch
             write(6,9007) atom(i), pos(:,i), coord(:,i)
          end do
          write(6,*)
          coordmag= matmul(R,posmag)
          write(6,*) ' Coordonnees reelles des centres magnetiques'
          write(6,*) 'Atome    Positions fractionnaires      Coordonnees reelles'   
          do i=1,nmag 
             write(6,9007) atommag(i),posmag(:,i),coordmag(:,i)
          end do
          write(6,*)
       end if


       !
       !***********************************************************************  
       !
       !----- Recentrage des atomes en pos0
       !---   Verification qu'il n'y a pas de double comptage
       !

       do i= 1, nch
          pos(:,i)= pos(:,i) - pos0(:) !+ 0.5d0
       end do
       !--- Coord. at. entre -0.5 et 0.5
       do i = 1,nch
          do idim = 1,3
             if (pos(idim,i)+0.5d0 < 0.0d0) then 
                nn = floor(pos(idim,i)+0.5d0)
                if ( iprint>= 1 )   write(6,9021) i,atom(i),pos(:,i),idim,pos(idim,i)-dble(nn)
                pos(idim,i) = pos(idim,i) - dble(nn)
             end if
             if (pos(idim,i)+0.5d0 >= 1.0) then 
                nn = floor(pos(idim,i)+0.5d0)
                if ( iprint>= 1 )   write(6,9021) i,atom(i),pos(:,i),idim,pos(idim,i)-dble(nn)
                pos(idim,i) = pos(idim,i) - dble(nn)
             end if
          end do
       end do
       if (ydebug) write(6,*) 
       if (ydebug) write(6,*) ' Positions recentrees'
       if (ydebug) then 
          do i = 1,nch
             write(6,9001) atom(i), pos(1,i), pos(2,i), pos(3,i)
          end do
       end if

       !---      Verification qu'il n'y a pas de double comptage

       write(6,*)
       write(6,"(' Verification de double comptage dans la cellule unite')")
       nn = nch
       i= 1
       do while ( i<=nn )
          j= i+1
          do while ( j<=nn )
             tpol(:) =  abs(pos(:,i)-pos(:,j))
             if (tpol(1)*a <= lim2 .and. tpol(2)*b <= lim2 .and. tpol(3)*c <= lim2) then 
                write(6,9023) i,atom(i),pos(:,i),j,atom(j),pos(:,j) 
                atom(j:nn-1)  = atom(j+1:nn)
                az(j:nn-1)  = az(j+1:nn)
                pos(:,j:nn-1) = pos(:,j+1:nn)
                q(j:nn-1) = q(j+1:nn)
                nn = nn - 1
             else
                j= j+1
             end if
          end do
          i=i+1
       end do
       nch = nn 

       !--- coord at magn recentres en pos0 
       !--- +verif qu'ils soient dans une bonne cell
       write(6,*)
       if (ydebug) write(6,*) ' Positions centres magnetiques recentrees'
       do i= 1, nmag
          posmag(:,i) = posmag(:,i) - pos0(:) ! + 0.5d0
          if (ydebug) write(6,9001) atommag(i), posmag(:,i)
          do j= 1, 3
             if ( posmag(j,i)<-ncel .or. posmag(j,i)>ncel+1 ) then
                write(6,"(' ***** atome magnetique ',i2,' sort des cellules centrales *****')") i
                write(6,"(' ***** dans la direction ',i2,' : pos= ',F10.5,' *****')") j, posmag(j,i)
             end if
          end do
       end do
       !--- coord reelles
       coordmag= matmul(R,posmag)
       if (iprint >=1) then
          write(6,*)
          write(6,*) ' Coordonnees reelles des centres magnetiques'
          do i=1,nmag 
             write(6,9001) atommag(i),coordmag(:,i)
          end do
       end if

       !--- centre
       pos0(:) = 0.d0

       !
       !***********************************************************************  
       !
       !---  Verification de la neutralite de la maille
       !
       write(6,"(' pour la maille elementaire :')")
       tmp = 0.d0
       tmp = sum(q(1:nch))
       write(6,*)
       write(6,9024) tmp
       tpol(:) = 0.D0
       do i=1,3
          tpol(i) = sum((pos(i,1:nch)+0.5d0)*q(1:nch)) 
       end do
       write(6,9025) tpol(:)

!!$
!!$***********************************************************************  
!!$
!!$----- Sortie donnees
!!$
!!$       

       write(6,"(/x,70('='))")
       write(6,"(x,'=== donnees generales')")
       write(6,*)
       write(6,*) 'Maille elementaire :'
       write(6,9011) nch
       write(6,9012) a,b,c
       write(6,9013) alpha,beta,gamma
       write(6,*)
       write(6,"(' axes de reference pour la base de sortie :')")
       write(6,"('       ',7x,'a',14x,'b',14x,'c')")
       write(6,"('  a2 : ',3(4x,ES12.5))") a2
       write(6,"('  b2 : ',3(4x,ES12.5))") b2
       write(6,"('  c2 : ',3(4x,ES12.5))") c2
       write(6,*)
       write(6,"(x,'Centre et ',i4,' centres magnetiques :')") nmag 
       write(6,*) '              Atomes    coord. fractionnaires'
       write(6,9016) '   ', pos0
       do i=1,nmag
          write(6,9017) atommag(i), (posmag(j,i),j=1,3) 
       end do
       if (iprint.ge.1) then
          write(6,*) ' Atomes    coord. fractionnaires        charge'
          do i = 1,nch
             write(6,9001) atom(i), (pos(j,i),j=1,3), q(i)
          end do
          write(6,*)
       end if
       if ( SINGLE ) then
          write(6,*)
          write(6,"('Moments multipolaires � annuler :',4x,i2)") ndip
       end if
       write(6,*)
       write(6,9015) ncel
       write(6,*)
!!$        write(6,"(' Systeme, pseudos, environnement ')")
       write(6,9018) dsys
       write(6,9019) dpseud
       write(6,*)
       write(6,"(' Options generales, debug= ',a5)") ydebug
       write(6,"(' Options generales, impression= ',i4)") iprint


       !***********************************************************************  
       !
       !----- ouverture des fichiers de sortie
       !
       open(fi_sys,file=nfi_sys) 
       if (nmag.le.5) then 
          write(fi_sys,'("Atom",2x, &
               &" x (a.u.) ",1x," y (a.u.) ",1x," z (a.u.) ",1x, &
               &" x (A)   ",1x,"  y (A)   ",1x,"  z (A)   ",14x, &
               &" charge   ",1x," distances (A) ")')
       else
          write(fi_sys,'("Atom",2x, &
               &" x (a.u.) ",1x," y (a.u.) ",1x," z (a.u.) ",1x, &
               &" x (A)   ",1x,"  y (A)   ",1x,"  z (A)   ",14x, &
               &" charge   ",1x," distance min (A) ")')
       end if

       open(fi_ps,file=nfi_ps)
       if (nmag.le.5) then 
          write(fi_ps,'("Atom",2x, & 
               &" x (a.u.) ",1x," y (a.u.) ",1x," z (a.u.) ",1x, &
               &" x (A)   ",1x,"  y (A)   ",1x,"  z (A)   ",14x, &
               &" charge   ",1x," distances (A) ")')
       else
          write(fi_ps,'("Atom",2x, & 
               &" x (a.u.) ",1x," y (a.u.) ",1x," z (a.u.) ",1x, &
               &" x (A)   ",1x,"  y (A)   ",1x,"  z (A)   ",14x, &
               &" charge   ",1x," distance min (A) ")')
       end if

       open(fi_env,file=nfi_env)
       write(fi_env,'(3x," x (a.u.) ",3x, " y (a.u.) ",3x," z (a.u.) ",3x, &
            &" charge  ",3x, &
            &"  xdip    ",3x,"  ydip    ",3x,"  zdip    ",3x)')

       open(fi_bid,file=nfi_bid)

       if (sortie.eq.'orca') then
          open(fi_xyz,file=nfi_xyz)
          open(fi_pc,file=nfi_pc)
       else
          open(fi_sew,file=nfi_sew)
       end if

!!$        open (unit=n,file=trim(nf),form="formatted",status="replace")
       open(unit=fi_fra,file=trim(nfi_fra),form='formatted',status='replace')
       write(fi_fra,"(' [Molden Format]'/' [5D]'/' [Atoms] (AU)')")

       open(unit=fi_fra2,file=trim(nfi_fra2),form='formatted',status='replace')
       write(fi_fra2,"(' [Molden Format]'/' [5D]'/' [Atoms] (AU)')")

       !***********************************************************************  
       !
       !----- atomes systemes de ref
       !
       write(6,"(/x,70('='))")
       write(6,"(x,'=== Atomes du systeme (utilises pour la convergence)')")

       allocate ( dist(nmag+1) )
       natsys= nch*(1+2*ncel)*(1+2*ncel)*(1+2*ncel)
       allocate ( coord_sys(3,natsys), atom_sys(natsys) )

       natsys= 0

       write(6,"(/' Atomes correspondant aux centres magnetiques')")
       write(6,"(' et autres atomes du systeme :')") 
       write(6,"('  nom        x         y         z          q      dist')")
       !---  cellules non renormalisees 
       if (ncel/=0) then
          do i = -ncel, ncel
             do j = -ncel, ncel
                do k = -ncel, ncel
                   do l= 1, nch
                      pos_1(1) = pos(1,l) + dble(i)
                      pos_1(2) = pos(2,l) + dble(j)
                      pos_1(3) = pos(3,l) + dble(k)
                      coord_1= matmul(R,pos_1)
                      dist(nmag+1) = dsqrt( &
                           coord_1(1)*coord_1(1) +  &
                           coord_1(2)*coord_1(2) +  &
                           coord_1(3)*coord_1(3))
                      do m= 1, nmag
                         dist(m) =  dsqrt( &
                              (coord_1(1)-coordmag(1,m))*(coord_1(1)-coordmag(1,m)) +&
                              (coord_1(2)-coordmag(2,m))*(coord_1(2)-coordmag(2,m)) +&
                              (coord_1(3)-coordmag(3,m))*(coord_1(3)-coordmag(3,m)))
                      end do
                      if ( minval(dist)<dsys .and. minval(dist)>lim2 ) then
                         natsys= natsys +1
                         coord_sys(:,natsys)= coord_1
                         atom_sys(natsys)= atom(l)
                         write(6,"(3x,a3,3x,3(f8.4,2x),2x,f6.3,4x,f6.2)") atom(l), pos_1(:), q(l), minval(dist)
                      end if
                      if ( minval(dist)<lim2 ) then
                         write(6,"(3x,a3,3x,3(f8.4,2x),2x,f6.3,3x,'(magn)')") atom(l), pos_1(:), q(l)
                      end if
                   end do
                end do
             end do
          end do
       end if

!!$        write(6,"(/' Autres atomes du systeme  : ',i4)")  natsys
!!$        write(6,"(100(3x,a3,3x,3(f8.4,2x)/))") (atom_sys(i),coord_sys(:,i),i=1,natsys) 

!!$***********************************************************************  
!!$
!!$----- boucle sur ndip
!!$

       write(6,"(x,70('='))")
       write(6,"(x,'=== iteration du calcul')")

       allocate ( pot(nmag+1+natsys,mindip:maxdip) )
       allocate ( dpot(nmag+1+natsys,mindip+1:maxdip) )
       allocate ( vpot(nmag+natsys,mindip:maxdip) )
       allocate ( dvpot(nmag+natsys,mindip+1:maxdip) )

       write(6,"(/' ndip   nbr atm')",advance='no')
       do i=1, nmag
          write(6,"('     VP',i2)",advance='no') i
       end do
       write(6,"('     VP C ')",advance='no')
       if ( natsys>0 )   write(6,"('    mVPsys')",advance='no')
       do i=1, nmag
          write(6,"('    VDP',i2)",advance='no') i
       end do
       if ( natsys>0 )   write(6,"('   mVDPsys')",advance='no')
       write(6,*)


       do ndip= mindip, maxdip !!-Debut boucle ndip------------------------

          !  dealloc
          if ( ndip > mindip ) then
             deallocate (lambda, mu)
             deallocate (pos_env,q_env)
             deallocate (coord_env,atom_env,az_env)
          end if

          natom = nch*(1+2*ncel+2*ndip)*(1+2*ncel+2*ndip)*(1+2*ncel+2*ndip)
          !
          !----- Calcul des facteurs de renormalisation des charges
          !
          allocate (lambda(nch,3,ndip+1), mu(nch,3))

          call fctmu(nch,ndip+1,pos+0.5d0,lambda)

          mu(:,:) = 0.0d0
          do p=1,ndip+1
             mu(1:nch,1:3) = mu(1:nch,1:3) + lambda(1:nch,1:3,p)
          end do
          do idim = 1,3
             do i = 1,nch
                if (abs(mu(i,idim)-1.0) >= 1.d-6) then
                   write(6,*)
                   write(6,*) ' erreur dans le calcul de mu'
                   write(6,*) i,idim,mu(i,idim), abs(mu(i,idim)-1.0d0)
                   write(6,*)
                   STOP
                end if
             end do
          end do

          if (ydebug) then 
             write(6,*)
             write(6,*) ' lambda(nch,3,ndip+1) '
             do k=1,ndip+1
                write(6,*) '   cellule ', k 
                do i=1,nch
                   write(6,*) lambda(i,:,k)
                end do
             end do
          end if

          !
          !----- Calcul de la renormalisation des charges et recentrage en (0,0,0)
          !
          allocate (pos_env(3,natom),q_env(natom))
          allocate (coord_env(3,natom),atom_env(natom),az_env(natom))

          !--- cellule centrale
          atom_env(1:nch)  = atom(1:nch)
          az_env(1:nch)  = az(1:nch)
          pos_env(:,1:nch) = pos(:,1:nch)
          q_env(1:nch) = q(1:nch)
          nn = nch
          if (ydebug) then 
             write(6,*)
             write(6,*) ' Maille centrale apres renormalisation des charges'
             do i= 1,nch
                write(6,9001) atom_env(i), pos_env(:,i), q_env(i)
             end do
          end if

          !--- cellules non renormalisees 
          if (ncel/=0) then
             do i = -ncel, ncel
                do j = -ncel, ncel
                   do k = -ncel, ncel
                      if (i/=0.or.j/=0.or.k/=0) then
                         atom_env(nn+1:nn+nch)  = atom(1:nch)
                         az_env(nn+1:nn+nch)  = az(1:nch)
                         pos_env(1,nn+1:nn+nch) = pos(1,1:nch) + dble(i)
                         pos_env(2,nn+1:nn+nch) = pos(2,1:nch) + dble(j)
                         pos_env(3,nn+1:nn+nch) = pos(3,1:nch) + dble(k)
                         q_env(nn+1:nn+nch) = q(1:nch)
                         if (ydebug) then 
                            write(6,*)
                            write(6,*) ' Maille non touchee apres renormalisation des charges ',i,j,k
                            do p = 1,nch
                               write(6,9001) atom_env(p+nn), pos_env(:,p+nn), q_env(p+nn)
                            end do
                         end if
                         nn = nn + nch
                      end if
                   end do
                end do
             end do
          end if


          !--- cellules renormalisees 
          if (ndip/=0) then
             do i = -ncel-ndip, ncel+ndip
                pmin(1) = max(1, -ncel+1+i)
                pmax(1) = min(ndip+1, ncel+ndip+1+i)
                do j = -ncel-ndip, ncel+ndip
                   pmin(2) = max(1, -ncel+1+j)
                   pmax(2) = min(ndip+1, ncel+ndip+1+j)
                   do k = -ncel-ndip, ncel+ndip
                      pmin(3) = max(1, -ncel+1+k)
                      pmax(3) = min(ndip+1, ncel+ndip+1+k)
                      if (abs(i)>ncel.or.abs(j)>ncel.or.abs(k)>ncel) then
                         atom_env(nn+1:nn+nch)  = atom(1:nch)
                         az_env(nn+1:nn+nch)  = az(1:nch)
                         pos_env(1,nn+1:nn+nch) = pos(1,1:nch) + dble(i) 
                         pos_env(2,nn+1:nn+nch) = pos(2,1:nch) + dble(j)
                         pos_env(3,nn+1:nn+nch) = pos(3,1:nch) + dble(k)

                         mu(:,:) = 0.0d0
                         do idim=1,3
                            do p = pmin(idim), pmax(idim)
                               mu(:,idim) = mu(:,idim) + lambda(:,idim,p)
                            end do
                         end do
                         if (ydebug) then
                            write(6,*) ' mu ', i,j,k
                            write(6,*) mu(:,1),mu(:,2),mu(:,3)
                         end if

                         q_env(nn+1:nn+nch) = q(1:nch)*mu(1:nch,1)*mu(1:nch,2)*mu(1:nch,3) 
                         if (ydebug) then 
                            write(6,*)
                            write(6,9022) i,j,k
                            do p = 1,nch
                               write(6,9001) atom_env(p+nn), pos_env(:,p+nn), q_env(p+nn)
                            end do
                         end if
                         nn = nn + nch
                      end if
                   end do
                end do
             end do
             if (nn /= natom) write(6,*) '****** ERREUR ****** nch=',&
                  nn,'/= natom=',natom
          end if

          !
          !----- mise en place de la sym�trie d'ordre 3
          !
          if (yaxe3) then
             write(6,*) 
             write(6,*) " Mise en place de la symetrie d'ordre 3"
             write(6,*) " Nbre d'atomes avant : natom =", natom

             allocate(atom_axe3(3*natom), az_axe3(3*natom), pos_axe3(3,3*natom), q_axe3(3*natom))

             call axe_sym3(natom,atom_env,az_env,pos_env,q_env,&
                  lim2,&
                  natom_axis3,atom_axe3,az_axe3,pos_axe3,q_axe3)

             deallocate (pos_env,q_env)
             deallocate (coord_env,atom_env,az_env)
             allocate (pos_env(3,natom_axis3), q_env(natom_axis3))
             allocate (coord_env(3,natom_axis3), atom_env(natom_axis3), az_env(natom_axis3))

             natom = natom_axis3
             atom_env(1:natom)    = atom_axe3(1:natom)
             az_env(1:natom)      = az_axe3(1:natom)
             pos_env(1:3,1:natom) = pos_axe3(1:3,1:natom)
             q_env(1:natom)       = q_axe3(1:natom)
             deallocate(atom_axe3, az_axe3, pos_axe3, q_axe3)
             write(6,*) " Nbre d'atomes apres : natom =", natom
          end if

          !
          !----- Test de la symetrie d'ordre 3
          !

          coord_env= matmul(R,pos_env)
          if (yaxe3) call  test_sym3(natom,atom_env,coord_env,q_env,lim2,lim4)

          !
          !----- Elimination des atomes de charge nulle 
          !----- Sauf pour systeme et pseudos
          !

          nn= natom
          iatom= 1
          do while ( iatom<= nn )
             !--- distance aux at. magn.
             dist(:)  = 0.d0
             dist_min = 0.d0
             dist(nmag+1) = dsqrt(coord_env(1,iatom)*coord_env(1,iatom) +  &
                  coord_env(2,iatom)*coord_env(2,iatom) +  &
                  coord_env(3,iatom)*coord_env(3,iatom))
             do i = 1,nmag
                dist(i) =  dsqrt( &
                     (coord_env(1,iatom)-coordmag(1,i))*(coord_env(1,iatom)-coordmag(1,i)) +&
                     (coord_env(2,iatom)-coordmag(2,i))*(coord_env(2,iatom)-coordmag(2,i)) +&
                     (coord_env(3,iatom)-coordmag(3,i))*(coord_env(3,iatom)-coordmag(3,i)))
             end do
             dist_min= minval(dist)
             !--- syst.
             if (dist_min <= dsys) then 
                if (ydebug) write(6,'(1x,a3,1x,3(f10.4,1x),x,f10.6,1x,a5)') &
                     atom_env(iatom),coord_env(:,iatom),q_env(iatom), "  sys"
                iatom= iatom + 1
                goto 1001
             end if
             !--- pseud.
             if (dist_min <= dsys + dpseud) then
                if (ydebug)  write(6,'(1x,a3,1x,3(f10.4,1x),x,f10.6,1x,a5)') &
                     atom_env(iatom),coord_env(:,iatom),q_env(iatom), "pseud"
                iatom= iatom + 1
                goto 1001
             end if
             !--- autre elimination des atomes de charge nulle
             if (abs(q_env(iatom)) < lim3) then
                nn = nn -1
                if (iprint >= 2) write(6,9041) iatom, atom_env(iatom),&
                     pos_env(:,iatom),q_env(iatom)
                atom_env(iatom:natom-1) = atom_env(iatom+1:natom)
                az_env(iatom:natom-1) = az_env(iatom+1:natom)
                pos_env(:,iatom:natom-1) = pos_env(:,iatom+1:natom)
                coord_env(:,iatom:natom-1) = coord_env(:,iatom+1:natom)
                q_env(iatom:natom-1) = q_env(iatom+1:natom)
             else
                iatom= iatom + 1
             end if
1001         continue
          end do
          natom = nn

          !
          !----- Calcul des positions dans l'espace reel
          !
          coord_env= matmul(R,pos_env)

          !
          !----- Calcul du potentiel electro
          !      des atommag et du centre 

          pot(:,ndip) = 0.d0

          do iatom = 1,natom
             !--- centre
             dist_1 = dsqrt(coord_env(1,iatom)*coord_env(1,iatom) +  &
                  coord_env(2,iatom)*coord_env(2,iatom) +  &
                  coord_env(3,iatom)*coord_env(3,iatom))
             if (dist_1>=lim5) then
                dist_1 = dist_1/ua
                pot(nmag+1,ndip) = pot(nmag+1,ndip) + q_env(iatom)/dist_1
             end if
             !--- atomes magnetiques
             do i = 1,nmag
                dist_1 =  sqrt( &
                     (coord_env(1,iatom)-coordmag(1,i))*(coord_env(1,iatom)-coordmag(1,i)) +&
                     (coord_env(2,iatom)-coordmag(2,i))*(coord_env(2,iatom)-coordmag(2,i)) +&
                     (coord_env(3,iatom)-coordmag(3,i))*(coord_env(3,iatom)-coordmag(3,i)))
                if (dist_1>=lim5) then
                   dist_1 = dist_1/ua
                   pot(i,ndip) = pot(i,ndip) + q_env(iatom)/dist_1
                end if
             end do
             !--- autres atomes du systeme
             do i= 1, natsys
                dist_1 =  sqrt( &
                     (coord_env(1,iatom)-coord_sys(1,i))*(coord_env(1,iatom)-coord_sys(1,i)) +&
                     (coord_env(2,iatom)-coord_sys(2,i))*(coord_env(2,iatom)-coord_sys(2,i)) +&
                     (coord_env(3,iatom)-coord_sys(3,i))*(coord_env(3,iatom)-coord_sys(3,i)))
                if (dist_1>=lim5) then
                   dist_1 = dist_1/ua
                   pot(nmag+1+i,ndip) = pot(nmag+1+i,ndip) + q_env(iatom)/dist_1
                end if
             end do
          end do


          !  en mev
          pot(:,ndip) = pot(:,ndip)*mev
          !---  difference de potentiel
          ! avec magnetiques
          do i= 1, nmag
             vpot(i,ndip)= pot(i,ndip)-pot(nmag+1,ndip)
          end do
          ! et systeme
          do i= nmag+1, nmag+natsys
             vpot(i,ndip)= pot(i+1,ndip)-pot(nmag+1,ndip)
          end do


          !***********************************************************************  
          !
          !----- sortie de la boucle
          !
          ! precision
          if ( ndip>mindip ) then
             minpot= 0.D0
             do i= 1, nmag+1+natsys
                dpot(i,ndip)= pot(i,ndip)-pot(i,ndip-1)
                if ( abs(dpot(i,ndip))>minpot )   minpot= abs(dpot(i,ndip))
             end do
             minvpot= 0.D0
             do i= 1, nmag+natsys
                dvpot(i,ndip)= vpot(i,ndip)-vpot(i,ndip-1)
                if ( abs(dvpot(i,ndip))>minvpot )   minvpot= abs(dvpot(i,ndip))
             end do

             if ( natsys>0 )   then
                maxsysdpot= maxval(abs(dpot(nmag+2:nmag+natsys+1,ndip)))             
                maxsysdvpot= maxval(abs(dvpot(nmag+1:nmag+natsys,ndip)))             
                write(6,"(x,i2,3x,ES8.1,2x,100(2x,ES7.0))") &
                     ndip, dble(natom),dpot(1:nmag+1,ndip),maxsysdpot,dvpot(1:nmag,ndip),maxsysdvpot
             else
                write(6,"(x,i2,3x,ES8.1,2x,100(2x,ES7.0))") ndip, dble(natom),dpot(:,ndip),dvpot(:,ndip)
             end if
             if ( minvpot<prec ) exit
          end if

          ! nombre d''atomes
          iatom = nch*(1+2*ncel+2*(ndip+1))*(1+2*ncel+2*(ndip+1))*(1+2*ncel+2*(ndip+1))
          if ( iatom>maxatom ) then
             write(6,"(/60('-')/'----- nombre d''atomes trop grand pour la prochaine iteration')")
             write(6,"('----- nombre d''atomes : ',i10,4x,'maximum : ',i10/60('-')/)") iatom, maxatom
             exit
          end if

       end do !!-Fin boucle ndip-----------------------------------------------

       ! si ca n'a pas converge :
       if ( ndip==maxdip+1 )  ndip= maxdip


       write(6,"(/x,70('='))")
       write(6,"(x,'=== valeures finales'/' ===  ')")
       write(6,"(' ===    moment dipolaires annules  : ', i4/' ===  ')") ndip
       write(6,"(' ===    nombre de charges          : ', i8/' ===  ')") natom
       if ( .not. SINGLE )  write(6,"(' ===    precision sur le potentiel  : '&
            &, ES8.1,'  meV'/' ===  ')") minpot
       if ( .not. SINGLE )  write(6,"(' ===    precision sur la difference : '&
            &, ES8.1,'  meV'/' ===  ')") minvpot

       !
       !***********************************************************************  
       !
       !----- Calcul des moments multipolaires
       !

       write(6,"(x,70('='))")
       write(6,"(x,'=== moments multipolaires'/)")
       do adip= 1,ndip+adipm
          write(6,"(' moment d''ordre :  ',i2,2x)",advance='no') adip
          if ( iprint > 0 )  write(6,*)
          mpolm= 0.D0
          do k=0,adip
             do j=0,adip
                do i=0,adip
                   if (i+j+k == adip) then
                      mpol = 0.0d0
                      do iatom = 1,natom
                         tmp = 1.d0
                         do  kk=1,k
                            tmp = tmp * pos_env(3,iatom)
                         end do
                         do jj = 1,j
                            tmp = tmp * pos_env(2,iatom)
                         end do
                         do ii = 1,i
                            tmp = tmp * pos_env(1,iatom)
                         end do
                         mpol = mpol + q_env(iatom) * tmp 
                      end do

                      if ( dabs(mpol) > mpolm )   mpolm= dabs(mpol)
                      if ( iprint > 0 ) then
                         if ( dabs(mpol) > 1.D-13 ) then
                            write(6,9003) i,j,k,mpol
                         else 
                            write(6,"(' ',3(i2,2x),'<    10^-13')") i,j,k
                         end if
                      end if
                   end if
                end do
             end do
          end do
          write(6,"('   valeur absolue maximale : ',ES10.3)") mpolm
       end do

!!$
!!$***********************************************************************  
!!$
!!$----- Calcul du potentiel electrostatique 
!!$
!!$       allocate(dist(nmag+1))
!!$       pot(:) = 0.d0
!!$       
!!$       do iatom = 1,natom
!!$          !  centre
!!$          dist(nmag+1) = dsqrt(coord_env(1,iatom)*coord_env(1,iatom) +  &
!!$               coord_env(2,iatom)*coord_env(2,iatom) +  &
!!$               coord_env(3,iatom)*coord_env(3,iatom))
!!$          if (dist(nmag+1)>=lim5) then
!!$             dist(nmag+1) = dist(nmag+1)/ua
!!$             pot(nmag+1) = pot(nmag+1) + q_env(iatom)/dist(nmag+1)
!!$          end if
!!$          !  atomes magnetiques
!!$          do i = 1,nmag
!!$             dist(i) =  sqrt( &
!!$                  (coord_env(1,iatom)-coordmag(1,i))*(coord_env(1,iatom)-coordmag(1,i)) +&
!!$                  (coord_env(2,iatom)-coordmag(2,i))*(coord_env(2,iatom)-coordmag(2,i)) +&
!!$                  (coord_env(3,iatom)-coordmag(3,i))*(coord_env(3,iatom)-coordmag(3,i)))
!!$             if (dabs(dist(i))>=lim5) then
!!$                dist(i) = dist(i)/ua
!!$                pot(i) = pot(i) + q_env(iatom)/dist(i)
!!$             end if
!!$          end do
!!$       end do
!!$          pot(:) = pot(:)*mev


       write(6,"(/x,70('='))")
       write(6,"(x,'=== potentiel electrostatique')")
       write(6,'(/" Potentiel electrostatique ",/," Atomes magnetiques :")') 
       do  i = 1,nmag
          write(6,9050) atommag(i), coordmag(:,i), pot(i,ndip)
       end do
       write(6,'(" Centre :")') 
       write(6,9050)  'ctr', pos0(:), pot(nmag+1,ndip)
       write(6,'(" Reste du systeme :")') 
       do  i = 1,natsys
          write(6,9050) atom_sys(i), coord_sys(:,i), pot(nmag+1+i,ndip)
       end do
       write(6,*)
       write(6,'(" Differences de potentiel :")')
       do i= 1, nmag
          write(6,9051)  atommag(i), coordmag(:,i), 'ctr', pos0(:), &
               pot(i,ndip)-pot(nmag+1,ndip)
       end do
       do i = 1, natsys
          write(6,9051)  atom_sys(i), coord_sys(:,i), 'ctr', pos0(:), &
               pot(nmag+1+i,ndip)-pot(nmag+1,ndip)
       end do
       do i = 1,nmag
          do j =  i+1, nmag
             write(6,9051) atommag(j), coordmag(:,j), atommag(i), coordmag(:,i), &
                  pot(j,ndip)-pot(i,ndip)
          end do
          do j =  1, natsys
             write(6,9051) atom_sys(j), coord_sys(:,j), atommag(i), coordmag(:,i), &
                  pot(nmag+1+j,ndip)-pot(i,ndip)
          end do
       end do
       do i = 1,natsys
          do j= i+1, natsys
             write(6,9051) atom_sys(j), coord_sys(:,j), atom_sys(i), coord_sys(:,i), &
                  pot(nmag+1+j,ndip)-pot(nmag+1+i,ndip)
          end do
       end do

       !
       !***********************************************************************  
       !
       !----- Separation sys-pseud-env
       !
       write(6,"(/x,70('='))")
       write(6,"(x,'=== ecritures atomes sans symetrie')")
       natfra=0
       natfra2=0
       natenv=0
       allocate(typat(natom))  ! type d'atomes pour ecriture seward (1 syst, 2 psd, 3 env)
       typat(:)=0

       do iatom = 1,natom
          !--- distance aux at. magn.
          dist(nmag+1) = dsqrt(coord_env(1,iatom)*coord_env(1,iatom) +  &
               coord_env(2,iatom)*coord_env(2,iatom) +  &
               coord_env(3,iatom)*coord_env(3,iatom))
          do i = 1,nmag
             dist(i) =  dsqrt( &
                  (coord_env(1,iatom)-coordmag(1,i))*(coord_env(1,iatom)-coordmag(1,i)) +&
                  (coord_env(2,iatom)-coordmag(2,i))*(coord_env(2,iatom)-coordmag(2,i)) +&
                  (coord_env(3,iatom)-coordmag(3,i))*(coord_env(3,iatom)-coordmag(3,i)))
          end do
          write (fi_bid,9004) atom_env(iatom),coord_env(:,iatom), q_env(iatom), dist(:)
          dist_min= minval(dist)
          !--- distance aux at. du syst. non renormalises
!!$           dist_min_sys= 1.D99
          dist_min_sys= dist_min
          do i= 1, natsys
             dist_1=  dsqrt( &
                  (coord_env(1,iatom)-coord_sys(1,i))*(coord_env(1,iatom)-coord_sys(1,i)) +&
                  (coord_env(2,iatom)-coord_sys(2,i))*(coord_env(2,iatom)-coord_sys(2,i)) +&
                  (coord_env(3,iatom)-coord_sys(3,i))*(coord_env(3,iatom)-coord_sys(3,i)))
             if ( dist_1 < dist_min_sys )    dist_min_sys= dist_1
          end do
          !--- syst.
          if (dist_min <= dsys) then
             typat(iatom) = 1
             if (nmag.le.5) then 
                write(fi_sys,9005) atom_env(iatom), &
                     coord_env(:,iatom)/ua, coord_env(:,iatom),&
                     q_env(iatom),dist(:)
             else
                write(fi_sys,9005) atom_env(iatom), &
                     coord_env(:,iatom)/ua, coord_env(:,iatom),&
                     q_env(iatom),dist_min
             end if
             natfra= natfra+1
             write(fi_fra,"(a3,4x,i3,4x,i3,3(4x,f9.5))") atom_env(iatom), natfra, az_env(iatom), &
                  coord_env(:,iatom)/ua
             !--- pseudo
          else if (dist_min_sys <= dpseud) then 
             typat(iatom) = 2
             if (nmag.le.5) then 
                write (fi_ps,9005) atom_env(iatom), &
                     coord_env(:,iatom)/ua, coord_env(:,iatom), &
                     q_env(iatom), dist(:), dist_min_sys
             else
                write (fi_ps,9005) atom_env(iatom), &
                     coord_env(:,iatom)/ua, coord_env(:,iatom), &
                     q_env(iatom),  dist_min_sys
             end if
             natfra2= natfra2+1
             write(fi_fra2,"(a3,4x,i3,4x,i3,3(4x,f9.5))") atom_env(iatom), natfra2, az_env(iatom), &
                  coord_env(:,iatom)/ua
!!$              OK= .false.
             !--- ou env
          else
             natenv= natenv+1
             typat(iatom) = 3
             write (fi_env,9006) coord_env(:,iatom)/ua, q_env(iatom),zero,zero,zero
          end if

       end do

       !--- fin fichier molekel
       write(fi_fra,"(' [GTO] (AU)')")
       do j= 1, natfra
          write(fi_fra,"(x,i3,'   0'/'   s   0   1.00'/' ')") j
       end do
       write(fi_fra,"(' [MO]'/x/x)")

       close(fi_fra)

       write(fi_fra2,"(' [GTO] (AU)')")
       do j= 1, natfra2
          write(fi_fra2,"(x,i3,'   0'/'   s   0   1.00'/' ')") j
       end do
       write(fi_fra2,"(' [MO]'/x/x)")

       close(fi_fra2)

       !***********************************************************************  
       !
       !----- Symetrie
       !

       if ( ysym ) then
          write (6,*) "call symetrie",  natom, nmag, maxnmag, natsys
          flush(6)
          deallocate (typat)
          call symetrie(natom, nmag, maxnmag, natsys, &
            atom_env, pos_env, coord_env, q_env, &
            coordmag, coord_sys, & 
            DELSYM,ndelsym,opdelsym, &
            dsys, dpseud,&
            fi_sys, fi_ps, fi_env, fi_sew, &
            natomirr, nirrsys, nirrps, nirrenv, &
            sortie, iprint, lim1,lim2, lim4)
          write (6,*) "apres symetrie"
          flush(6)
       else
          natomirr = natom
          if (natfra .ne. natsys+nmag) write(6,*) ">>> pb comptage <<<", &
               natfra, " different de ", natsys+nmag
          nirrsys  = natfra
          nirrps   = natfra2 
          nirrenv  = natenv
       end if


       !--- ecriture du fichier pour seward si pas de symetrie
       if (sortie.eq.'molcas'.and. .not.ysym ) then 
          write(fi_sew,"(' Fragment :')")
          i=0
          do iatom = 1,natom
             if ( typat(iatom)==1 ) then
                i= i+1
!!$                 write (fi_sew,"(2x,a,i2.2,x,3(4x,f17.12),8x,'Bohr')") &
!!$                      trim(atom_env(iatom)), i, coord_env(:,iatom)/ua
                if (atom_env(iatom)(2:2)==" ")  atom_env(iatom)(2:2)='_'
                write (fi_sew,"(2x,a2,i2.2,3(4x,f17.12),8x,'Bohr')") &
                     atom_env(iatom)(1:2), i, coord_env(:,iatom)/ua
             end if
          end do
          write(fi_sew,"(' Pseudos :')")
!!$           i=0
          do iatom = 1,natom
             if ( typat(iatom)==2 ) then
                i= i+1
!!$                 write (fi_sew,"(2x,a,i3.3,3(4x,f17.12),8x,'Bohr')") &
!!$                      trim(atom_env(iatom)), i, coord_env(:,iatom)/ua
                if (atom_env(iatom)(2:2)==" ")  atom_env(iatom)(2:2)='_'
                write (fi_sew,"(2x,a2,i2.2,3(4x,f17.12),8x,'Bohr')") &
                     atom_env(iatom)(1:2), i, coord_env(:,iatom)/ua
             end if
          end do
          write(fi_sew,"('XFIEld')")
          write(fi_sew,"(i8)") natenv
          i=0
          do iatom = 1,natom
             if ( typat(iatom)==3 ) then
                i= i+1
                write (fi_sew,"(2x,3(f17.12,4x),2x,f16.12,4x,'0.0   0.0   0.0')") &
                     coord_env(:,iatom)/ua, q_env(iatom)
             end if
          end do
       end if

       !--- ecriture du fichier pour orca
       if (sortie.eq.'orca') then 
          write(fi_xyz,"(i7)") nirrsys + nirrps 
          write(fi_xyz,*)  "Quantum fragment + Embedding potential  ", &
               trim(prefix)," Angstroms"
          i=0
          do iatom = 1,natom
             if ( typat(iatom)==1 ) then
                i= i+1
                write (fi_xyz,"(2x,a2,3(4x,f17.12))") &
                      lat(az_env(iatom)), coord_env(:,iatom)
             end if
          end do
!!$           i=0
          do iatom = 1,natom
             if ( typat(iatom)==2 ) then
                i= i+1
                write (fi_xyz,"(2x,a3,4x,f17.12, 4x,3(2x,f17.12))") &
                     trim(lat(az_env(iatom)))//'>', q_env(iatom), coord_env(:,iatom)
             end if
          end do
          !
          write(fi_pc,"(i8)") natenv
          i=0
          do iatom = 1,natom
             if ( typat(iatom)==3 ) then
                i= i+1
                write (fi_pc,"(2x,f16.12,4x,3(f17.12,4x))") &
                     q_env(iatom), coord_env(:,iatom)
             end if
          end do
       end if

       !--- fin de la partie sur la sym
!!$        end if



       write(6,*)
       write(6,*) 'Nbre total d''atomes = ',natom
       write(6,9042) natomirr,nirrsys,nirrps,nirrenv
       write(6,*)
       write(6,*) ' Ouf! Tout c''est probablement bien fini. '
       write(6,*)
       flush(6)

!!$***********************************************************************  
!!$
!!$----- Formats
!!$
9001   format (2x,a3,3x,3(f8.4,2x),ES12.5)
9002   format ('                      magnetiques=',1x,a3,3x,3(f8.4,2x))
9003   format (' Moment mutipolaire ',3(i2,2x),'= ',ES12.5)
9004   format (1x,a3,3x,3(f10.4,2x),x,f10.6,4x,10(f10.6,x))
9005   format (1x,a3,1x,6(f10.6,1x),1x,'Angstrom',2x,10(f10.6,x))
9006   format (3x,7(f15.6,3x))
9007   format (2x,a3,3x,6(f8.4,2x),f15.10)
9008   format (2x,i8,2x,a3,3x,4(f10.6,2x),i3)


9011   format('  Nbre d''atomes=',i5)
9012   format('  a, b, c            =',3(f10.4,2x)) 
9013   format('  alpha, beta, gamma =',3(f10.4,2x)) 
9014   format('  Nbre a annuler      =',  i3)
9015   format(' Nbre cell. intactes =',  i3)
9016   format('  Centre     :', 1x,a3,3x,3(f8.4,2x))
9017   format('  Magnetique :', 1x,a3,3x,3(f8.4,2x))
9018   format(' Systeme a ',f8.4,' du centre ou des atomes magnetiques')
9019   format(' Pseudos a ',f8.4,' du centre ou des atomes magnetiques')

9021   format(' Atome ',i3,1x,a3,' repositionne en ',3(f8.4,x), &
            ' selon la direction ', i3,' soit ', f8.4) 
9022   format(' Maille renormalisee apres renormalisation des charges ',3(i5,x))
9023   format(' Atomes ',i5,1x,a3,1x,3(f8.4,1x),' et ',i5,1x,a3,1x,3(f8.4,1x),&
            ' identiques, second supprime')
9024   format(' charge           =',f10.6)
9025   format(' moment dipolaire =',3(f10.6,2x))

9031   format(' Operation de symetrie     valide = ',3(I3,x))
9032   format('  Nbre maximum d''peration de symetrie = ', I3)
9033   format(' Operation de symetrie non valide = ',3(I3,x))

9041   format(' Atome elimine :',i9,2x,a3,x,3(f8.4,x),x,f10.6)
9042   format(' Nbre d''atomes irreductibles = ',I11,' (',2(i7,','),i7,')')

9050   format (2x,a3,1x,3(f8.4,2x),3x,f15.6)
9051   format (2x,2(a3,1x,3(f8.4,2x),2x),3x,f20.6)


!!$
!!$***********************************************************************  
!!$
End program Env


