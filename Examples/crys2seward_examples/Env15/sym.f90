  subroutine symetrie(natom, nmag, maxnmag, natomsys, &
       atom_env, pos_env, coord_env, q_env, &
       coordmag, coord_sys, & 
       DELSYM,ndelsym,opdelsym, &
       dsys, dpseud,&
       fi_sys, fi_ps, fi_env, fi_sew, &
       natomirr, nirrsys, nirrps, nirrenv, &
       sortie, iprint, lim1,lim2, lim4)
      Implicit none
!!$***************************************************************************************
!!$
!!$    Gestion de la symmetrie
!!$
!!$  
!!$***************************************************************************************
       real*8,  parameter :: pi=3.14159265358979d0, ua=0.52917720859d0, &
            mev=27211.39628d0, zero=0.0d0

      !--- Variables en entree
      integer, intent(in) :: natom, nmag, maxnmag, natomsys            ! nbre total d'atomes
      real*8, intent(in) :: lim1, lim2, lim4
      integer, intent(in) :: iprint
      integer, intent(in) :: fi_sys, fi_ps, fi_env, fi_sew
      character(len=6), intent(in) :: sortie
      real*8, intent(in) :: dsys, dpseud

      !--- Variables en entree pour virer des symetries
      logical, intent(inout) :: DELSYM
      integer, intent(inout) :: ndelsym
      Integer, dimension(3,8), intent(inout) :: opdelsym

      character*3, dimension(natom), intent(inout) :: atom_env ! nom des atomes
      real*8, dimension(3,natom), intent(in) :: pos_env     ! positions fractionnaires 
      real*8, dimension(3,natom), intent(in) :: coord_env   ! positions r�elles 
      real*8, dimension(natom), intent(in)   :: q_env       ! charge 
      real*8, dimension(3,maxnmag) :: coordmag  ! coord.  atomes mag.
      real*8, dimension(3,natomsys+nmag) ::  coord_sys

      !--- Variables en sortie
      integer, intent(out) :: natomirr, nirrsys, nirrps, nirrenv


      !--- Variables locales
      Integer, dimension(natom) :: sym,  typat
      Logical, dimension(natom) :: iw, ietoile
      Integer :: nopsym, maxsym, nopgpe
      Integer, dimension(3,8) :: opsym, opgpe
      Integer, dimension(8,natom) :: isym
      Real*8, dimension(3) :: coordsym
      Real*8, dimension(:,:), allocatable :: coordetoile
      logical, dimension(3) :: psym, asym      ! psym : plans de sym, asym : axes de sym
      logical :: ptsym, oksym, yaxe3=.false.
      real*8, dimension(nmag+1) :: dist
      integer :: natsys


      !--- Variables muettes
      integer :: i, j, k, ii, nn, is, p
      integer :: iatom, jatom
      logical :: OK
      logical, dimension(:), allocatable :: OKatomsym
      real*8 :: dist_1, dist_min, dist_min_sys

!!$
!!$
!!$***********************************************************************  
!!$

      natsys = natomsys + nmag
      write (6,*) "entree symetrie",  natom, nmag, maxnmag, natsys
      flush(6)

      typat(:)  = 0  ! type d'atomes pour ecriture seward (1 syst, 2 psd, 3 env)
      nopsym= 1
      opsym(:,1)= (/1,1,1/)
      oksym= .false.

      sym(:)    = 0  ! nbre de symetriques      (sym(i)== nbr copain de i)
      isym(:,:) = 0  ! numero des symetriques   (isym(j,i)== jeme copain de i)
      write(6,"(/x,70('='))")
      write(6,"(x,'=== detection  atomes symetriques')")

      !--- Initialisations 
          
      !--- Calcul du nbre de symetriques par atome
      maxsym = 1     ! nbre maximum d'operations de symetrie
      ii= 0
      do i=1,natom
         do j=i+1,natom
            if ( sym(i)==7 )  exit ! tous les symetriques possibles sont deja trouves
            if ( atom_env(i)==atom_env(j) ) then
               coordsym(:) = abs(coord_env(:,j)) - abs(coord_env(:,i))
               if ((abs(coordsym(1)) < lim2) .and. (abs(coordsym(2)) < lim2) &
                    .and.(abs(coordsym(3)) < lim2)) then
                  if (abs(q_env(i)-q_env(j)) < lim1) then  
                     sym(i) = sym(i) +1
                     sym(j) = sym(j) +1
                     if (maxsym-1 < sym(i))  maxsym= maxsym + 1
                     if (maxsym-1 < sym(j))  maxsym= maxsym + 1
                     isym(sym(i),i) = j
                     isym(sym(j),j) = i
                  else if (ii<10) then
                     write(6,*) ' Atomes ',i,' et ',j,' sont symetriques avec des charges differentes'
                     write(6,9007) atom_env(i), pos_env(:,i), coord_env(:,i), q_env(i)
                     write(6,9007) atom_env(j), pos_env(:,j), coord_env(:,j), q_env(j)
                     ii= ii+1
                     if (ii==10) write(6,"(/'...'/)")
                  end if
               end if
            end if
         end do
      end do

      write(6,"(/' nbre maximal d''operations de symetries : ',i2)") maxsym

      !--- Identification des operations de symetrie valides
      write(6,"(/x,70('='))")
      write(6,"(x,'=== test des symetries'/)")
!!$           call flush(6)
      if (iprint >= 1) write(6,9032) maxsym


      nopsym = 0
      do i = 1,-1,-2      ! boucle sur les symmetries (i,j,k)
      do j = 1,-1,-2
      do k = 1,-1,-2   
         OK= .true.
         !--- Test sur tous les atomes pour voir si l'operation (i,j,k) est valide 
         do iatom = 1,natom
            coordsym(1) =  dble(i)*coord_env(1,iatom)
            coordsym(2) =  dble(j)*coord_env(2,iatom)
            coordsym(3) =  dble(k)*coord_env(3,iatom)
            is = 0
            !---  l'atome est-il son propre symetrique pour la symmetrie (i,j,k)           
            if ( (abs(coord_env(1,iatom)-coordsym(1)) < lim4).and.&
                 (abs(coord_env(2,iatom)-coordsym(2)) < lim4).and.&
                 (abs(coord_env(3,iatom)-coordsym(3)) < lim4)) then 
               ! oui il est son symmetrique par (i,j,k)
               is = is + 1
            end if
            !---  l'atome a-il un symetrique different de lui pour la symmetrie (i,j,k)
            do  p = 1,sym(iatom)
               if ( (abs(coord_env(1,isym(p,iatom))-coordsym(1)) < lim2).and.&
                    (abs(coord_env(2,isym(p,iatom))-coordsym(2)) < lim2).and.&
                    (abs(coord_env(3,isym(p,iatom))-coordsym(3)) < lim2)) then 
                  ! oui il a un symmetrique par (i,j,k)
                  is = is + 1
               end if
            end do
            !--- si c'est pas le cas
            if (is/=1) then 
               if (is==0) then ! operation de symmetrie non valide pour le systeme
                  write(6,9033) i,j,k
                  write(6,9008) iatom,atom_env(iatom),pos_env(:,iatom),q_env(iatom)
                  write(6,9034) is
                  write(6,"('    autres atomes potentiellement symetriques de cet atome :')")
                  do  p = 1,sym(iatom)
                     write(6,9008) isym(p,iatom),atom_env(isym(p,iatom)),pos_env(:,isym(p,iatom)),&
                          q_env(isym(p,iatom)),is
                  end do
                  OK= .false.
                  exit
               end if
               write(6,*)
               write(6,*) ' ERREUR - ERREUR - ERREUR - ERREUR - ERREUR' 
               write(6,*) ' Operateur de symetrie',i,j,k
               write(6,9008) iatom,atom_env(iatom),coord_env(:,iatom),q_env(iatom),is
               OK= .false.
               exit
            end if
            ! end if
         end do
         if(OK) then  ! si (i,j,k) est une operation valide pour le systeme
            write(6,9031) i,j,k
            nopsym = nopsym + 1
            opsym(1,nopsym) = i
            opsym(2,nopsym) = j
            opsym(3,nopsym) = k
         end if
      end do
      end do
      end do

      write(6,"(/' nbre d''operations de symetries : ',i2)") nopsym

      !--- Suppression des sym
      if (DELSYM) then
         write(6,"(/' Suppression des symetries : ')") 
         do i = 1, nopsym
            !--- Faut-il enlever cette operation de symetrie 
            OK=.false.
            do j=1, ndelsym
               if ( opsym(1,i).eq.opdelsym(1,j) .and. opsym(2,i).eq.opdelsym(2,j) &
                    .and. opsym(3,i).eq.opdelsym(3,j) ) then
                  OK=.true.
                  exit
               end if
            end do
            if (OK) then  !--- On enleve l'operation de symetrie
               write(6,"('  op : ',3(x,i3))") opsym(:,i)
               do j=i+1, nopsym
                  opsym(:,j-1)=opsym(:,j)
               end do
               nopsym= nopsym -1
            end if
         end do
      end if
      write(6,"(/' nbre final d''operations de symetries : ',i2)") nopsym

      write(6,"(/x,70('='))")
      write(6,"(x,'=== symetries')")

      !--- pour les sym
      oksym= .false.

      !--- Quels sont les plans de symetrie
      psym(:)= .false.
      write(6,"(/' Plans de symetrie :')") 
      do i= 1, nopsym
         if (  sum(opsym(:,i)) == 1 ) then
            do j= 1, 3
               if ( opsym(j,i)==-1 ) then
                  oksym= .true.
                  psym(j)= .true.
                  write(6,"('  plan orth. a l''axe ',i1)") j
               end if
            end do
         end if
      end do
          

      !--- S'il n'y a pas de plan de symetrie 
      ptsym  = .false.
      asym(:)= .false.
      write(6,"(/' Symetries supplementaires :')")
      ii = 0 ! compte les axes de symetrie
      if ( (.not. psym(1)) .or. (.not. psym(2)) .or. (.not. psym(3)) ) then
         !--- point de sym ? 
         do i=1, nopsym
            if ( sum(opsym(:,i))==-3 ) then
               oksym= .true.
               ptsym= .true.
               write(6,"('  point de symetrie  ')")
            end if
         end do
         !--- axes d'ordre 2 ?
         do i= 1,  nopsym
            OK= .true.
            if ( sum(opsym(:,i))==-1 ) then
               !--- genere par des plans de sym !?
               do j= 1, 3
                  if ( opsym(j,i)==-1 .and. .not. psym(j) )  then
                     OK= .false.
                     exit
                  end if
               end do
               !--- ou le point de sym et un plan !?
               if ( .not. OK .and. ptsym ) then
                  do j= 1, 3
                     if ( opsym(j,i)==1 .and. psym(j) )  then
                        OK= .false.
                        exit
                     end if
                  end do
               end if
            end if
            !--- sinon ...
            if ( .not. OK ) then ! vrai axe de symetrie du systeme
               do j= 1, 3
                  if ( opsym(j,i)== 1 ) then
                     oksym= .true.
                     asym(j)= .true.
                     ii = ii + 1
                     write(6,"('  axe d''ordre 2 : ',i1)") j
                  end if
               end do
            end if
         end do
      end if

      if (ii.eq.3) then 
         write(6,*) " Attention il y a trois axes de symetrie il faut en enlever un"
         write(6,*) " On enleve l'axe z"
         DELSYM = .true. 
         ndelsym = ndelsym + 1
         opdelsym(1, ndelsym) = -1
         opdelsym(2, ndelsym) = -1
         opdelsym(3, ndelsym) =  1
         do i = 1, nopsym
            !--- numero de cette operation de symetrie 
            OK=.false.
            if ( opsym(1,i).eq.opdelsym(1,ndelsym) .and. opsym(2,i).eq.opdelsym(2,ndelsym) &
                 .and. opsym(3,i).eq.opdelsym(3,ndelsym) ) then
               OK=.true.
               exit
            end if
         end do
         if (OK) then  !--- On enleve l'operation de symetrie
            write(6,"('  op : ',3(x,i3))") opsym(:,i)
            do j=i+1, nopsym
               opsym(:,j-1)=opsym(:,j)
            end do
            nopsym= nopsym -1
            i= i -1
         end if
         write(6,"(/' nbre final d''operations de symetries : ',i2)") nopsym
      end if


      !    
      !--- Symetries finales   (levons nous et demain ...)
      !
      write(6,"(/' Symetries finales : ',i2)") nopsym
      do i= 1, nopsym
         write(6,"('  op ',i1,' : ',3(x,i3))") i, opsym(:,i)
      end do

      !    
      !--- ne garder pour chaque atome que les operations de symetrie finales
      !    pour chaque atome : rafraichissement des tableaux 
      !       sym(iatom)=n et isym(in) = jatom
      !
      sym(:) = 1
      do iatom = 1,natom
         isym(1,iatom) = iatom
      end do

      ! reconstruite toutes les operations du gpe et pas seulement les irreductibles
      nopgpe =  nopsym
      opgpe(1:3,1:nopsym) = opsym(1:3,1:nopsym)
      do i = 1,nopsym
         do j = 1, nopsym
            opgpe(:,nopgpe+1) = opsym(:,i)*opsym(:,j)
            do k = 1,nopgpe
               if ( opgpe(1,nopgpe+1).eq.opgpe(1,k) .and.  &
                    opgpe(2,nopgpe+1).eq.opgpe(2,k) .and. &
                    opgpe(3,nopgpe+1).eq.opgpe(3,k) ) goto 11
            end do
            nopgpe = nopgpe + 1
11          continue
         end do
      end do

      allocate(coordetoile(3,nopgpe), OKatomsym(nopgpe))

      iw(:) = .false.      
      do iatom = 1,natom
         ! si l'etoile a deja ete construite je saute
         if ( iw(iatom) ) cycle
         ! construction de l'etoile 
         do i = 1,nopgpe
            coordetoile(:,i) = opgpe(:,i)*coord_env(:,iatom)
            OKatomsym(i) = .true.
         end do
         ! On enleve les doubles pour cet atome
         ii = nopgpe
         do i = 1,nopgpe
            if (.not.OKatomsym(i)) cycle
            do j=i+1,nopgpe
               if ( (abs(coordetoile(1,i)-coordetoile(1,j)) < lim2).and.&
                    (abs(coordetoile(2,i)-coordetoile(2,j)) < lim2).and.&
                    (abs(coordetoile(3,i)-coordetoile(3,j)) < lim2)) then 
                  ! symetries i et j donnent le mm atome
                  ii = ii - 1
                  OKatomsym(j)= .false.
               end if
            end do
         end do
         ! ii : nombre de membre de l'etoile pour iatom 
         ! je cherche les symetriques de iatom
         do i = 1,nopgpe
            if (opgpe(1,i).eq.1 .and. opgpe(2,i).eq.1 .and. opgpe(3,i).eq.1) cycle
            if (.not.OKatomsym(i)) cycle
            ! je cherche le symetrique de iatom
            coordsym(:) = coordetoile(:,i)
            do jatom = 1,natom
               if ( (abs(coord_env(1,jatom)-coordsym(1)) < lim4).and.&
                    (abs(coord_env(2,jatom)-coordsym(2)) < lim4).and.&
                    (abs(coord_env(3,jatom)-coordsym(3)) < lim4)) then 
                  sym(iatom) = sym(iatom) + 1 
                  isym(sym(iatom),iatom) = jatom
                  goto 10
               end if
            end do
            write(6,*) " pas trouve le symmetrique ", iatom, opgpe(:,i)
            stop 
10          continue
         end do
         if (sym(iatom).ne.ii) stop "erreur nombre de symetriques"
         ! j'atribue aux autre membres de l'etoile les symetriques
         do i = 1,ii
            jatom = isym(i,iatom)
            if (jatom.ne.iatom.and.sym(jatom).ne.1) &
                 write(6,*) iatom, jatom, sym(iatom), sym(jatom) 
         end do
         do  i = 1,ii
            p = isym(i,iatom) 
            iw(p) = .true.
            do j = 2,ii ! de l'etoile de iatom
               jatom = isym(j,iatom) 
               if (i.ne.j) then 
                  sym(jatom) = sym(jatom) + 1
                  isym(sym(jatom),jatom) = p
               end if
            end do
         end do
      end do


      !--- Ecriture des atomes irreductibles 
      write(6,"(/x,70('='))")
      write(6,"(x,'=== ecritures atomes avec symetrie')")

      write(fi_sys,*) 
      write(fi_sys,*) 
      write(fi_sys,*) ' En symetrie'
      write(fi_ps,*) 
      write(fi_ps,*) 
      write(fi_ps,*) ' En symetrie'
      write(fi_env,*) 
      write(fi_env,*) 
      write(fi_env,*) ' En symetrie'

!       do iatom = 1,natom
!          write(6,*)  iatom, sym(iatom), (isym(p,iatom),p=1,sym(iatom))
!       end do


       !--- Je cherche quel atome de l'etoile j'ecris
       iw(:) = .false.    ! vrai = j'ecrit cet atome dans la classe d'equiv. de sym (etoile)
       ietoile(:) = .false. ! vrai = j'ai trouve l'atome a ecrire dans la classe d'equiv. de sym (etoile)
       do iatom = 1,natom 
          if (.not.ietoile(iatom)) then  ! si je n'ai pas encore choisi d'atome a ecrire dans son etoile 
             if (sym(iatom).eq.1) then !--- l'atome n'a pas de symmetrique autre que lui-meme je le mets
                iw(iatom) = .true.
                ietoile(iatom) = .true. 
                cycle
             else 
                !---  est ce qu'on met cet atome ou un de ses symetriques
                is = 0
                coordsym(:) =  coord_env(:,iatom)
                do p = 2, sym(iatom)  ! boucle sur tous les symetrieques de l'etoile pour les tester
                   ! je choisi l'atome le plus proche de l'origine
                   if (  sum(coord_env(:,isym(p,iatom))) < sum(coordsym(:)) ) is = p
                   coordsym(:) =  coord_env(:,isym(p,iatom))
                end do
                if (is.eq.0) then 
                   iw(iatom) = .true.
                   ietoile(iatom) = .true.
                else  
                   iw(isym(is,iatom)) = .true.
                   ietoile(iatom) = .true.
                end if
                do  p = 1,sym(iatom) 
                   ietoile(isym(p,iatom)) = .true.
                end do
             end if
          end if
       end do


       natomirr = 0
       nirrsys  = 0
       nirrps   = 0
       nirrenv  = 0

       do iatom = 1,natom 
          if (iw(iatom)) then
             !--- distance aux at. magn.
             dist(nmag+1) = sqrt(coord_env(1,iatom)*coord_env(1,iatom) +  &
                  coord_env(2,iatom)*coord_env(2,iatom) +  &
                  coord_env(3,iatom)*coord_env(3,iatom))
             dist(1:nmag) =  sqrt( &
                  (coord_env(1,iatom)-coordmag(1,1:nmag))*(coord_env(1,iatom)-coordmag(1,1:nmag)) +&
                  (coord_env(2,iatom)-coordmag(2,1:nmag))*(coord_env(2,iatom)-coordmag(2,1:nmag)) +&
                  (coord_env(3,iatom)-coordmag(3,1:nmag))*(coord_env(3,iatom)-coordmag(3,1:nmag)))
             dist_min= minval(dist)
             !--- distance aux at. du syst. non renormalises
!!$                    dist_min_sys= 1.D99
             dist_min_sys= dist_min
             do i= 1, natsys
                dist_1=  dsqrt( &
                     (coord_env(1,iatom)-coord_sys(1,i))*(coord_env(1,iatom)-coord_sys(1,i)) +&
                     (coord_env(2,iatom)-coord_sys(2,i))*(coord_env(2,iatom)-coord_sys(2,i)) +&
                     (coord_env(3,iatom)-coord_sys(3,i))*(coord_env(3,iatom)-coord_sys(3,i)))
                if ( dist_1 < dist_min_sys )    dist_min_sys= dist_1
             end do
             !--- syst
             if ( dist_min <= dsys) then 
                natomirr = natomirr + 1
                nirrsys  = nirrsys  + 1
                write(fi_sys,9005) atom_env(iatom), &
                     coord_env(:,iatom)/ua, coord_env(:,iatom),&
                     q_env(iatom), dist(:)
                write(6,*) " etoile de ", iatom, ':',(isym(p,iatom),p=1,sym(iatom))
                typat(iatom)= 1
                !--- pseudo
             else if ( dist_min_sys <= dpseud) then
                natomirr = natomirr + 1
                nirrps   = nirrps   + 1
                write (fi_ps,9005) atom_env(iatom), &
                     coord_env(:,iatom)/ua, coord_env(:,iatom), &
                     q_env(iatom),dist(:),dist_min_sys
                typat(iatom)= 2
                !--- ou env
             else
                natomirr = natomirr + 1
                nirrenv  = nirrenv  + 1
                write (fi_env,9006) coord_env(:,iatom)/ua, q_env(iatom),zero,zero,zero
                typat(iatom)= 3
             end if
          end if
       end do

       !--- ecriture du fichier pour seward
       if (sortie.eq.'molcas') then 
          if ( oksym ) then
             write(fi_sew,"('SYMMETRY')")
             if ( ptsym )    write(fi_sew,"(x,a)", advance='no') "xyz"
             if ( asym(1) )  write(fi_sew,"(x,a)", advance='no') "yz"
             if ( asym(2) )  write(fi_sew,"(x,a)", advance='no') "xz"
             if ( asym(3) )  write(fi_sew,"(x,a)", advance='no') "xy"
             if ( psym(1) )  write(fi_sew,"(x,a)", advance='no') "x"
             if ( psym(2) )  write(fi_sew,"(x,a)", advance='no') "y"
             if ( psym(3) )  write(fi_sew,"(x,a)", advance='no') "z"
             write(fi_sew,*)
          end if
          write(fi_sew,"(' Fragment :')")
          i=0
          do iatom = 1,natom
             if ( typat(iatom)==1 .and. iw(iatom).eqv. .true. ) then  ! j'ecris les atomes du systeme 
                i= i+1
                if (atom_env(iatom)(2:2)==" ")  atom_env(iatom)(2:2)='_'
                write (fi_sew,"(2x,a2,i2.2,3(4x,f17.12),8x,'Bohr')") &
                     atom_env(iatom)(1:2), i, coord_env(:,iatom)/ua
             end if
          end do
          write(fi_sew,"(' Pseudos :')")
          do iatom = 1,natom
             if (typat(iatom)==2 .and. iw(iatom).eqv..true.) then  ! j'ecris les atomes du pseudo
                i= i+1
                if (atom_env(iatom)(2:2)==" ")  atom_env(iatom)(2:2)='_'
                write (fi_sew,"(2x,a2,i2.2,3(4x,f17.12),8x,'Bohr')") &
                     atom_env(iatom)(1:2), i, coord_env(:,iatom)/ua
             end if
          end do
          write(fi_sew,"('XFIEld')")
          write(fi_sew,"(i8)") nirrenv
          i=0
          do iatom = 1,natom
             if ( typat(iatom)==3  .and. iw(iatom).eqv..true.) then ! j'ecris les atomes de l'env
                i= i+1
                write (fi_sew,"(2x,3(f17.12,4x),2x,f16.12,4x,'0.0   0.0   0.0')") &
                     coord_env(:,iatom)/ua, q_env(iatom)
             end if
          end do
       end if

       !--- formats
9005   format (1x,a3,1x,6(f10.6,1x),1x,'Angstrom',2x,10(f10.6,x))
9006   format (3x,7(f15.6,3x))
9007   format (2x,a3,3x,6(f8.4,2x),f15.10)
9008   format (2x,i8,2x,a3,3x,4(f10.6,2x),3x,i4)
9031   format(' Operation de symetrie     valide = ',3(I3,x))
9032   format('  Nbre maximum d''peration de symetrie = ', I3)
9033   format(' Operation de symetrie non valide = ',3(I3,x), /, '    car l''atome :')
9034   format('    a ',i3,' symetrique par cette operation')
       
       return
     end subroutine symetrie
