!-----Programme permettant de calculer les nouvelles positions--------!
!-----atomiques à partir de l'application d'un champ électrique-------!


      PROGRAM Position
      
      real*8, dimension(:,:), allocatable :: Hessien,Hessien2
      real*8, dimension(:,:), allocatable :: Tmp
      integer,dimension(:)  , allocatable :: Bid2,Tmp2
      real*8, dimension(:,:), allocatable :: Inverse,Hessienbidon
      real*8, dimension(:,:), allocatable :: Verification
      real*8, dimension(:,:), allocatable :: Born
      real*8, dimension(:,:), allocatable :: Position_atom
      real*8, dimension(:,:), allocatable :: Electric_field
      real*8, dimension(:),   allocatable :: Force
      real*8, dimension(:),   allocatable :: Deplacement
      real*8, dimension(:,:), allocatable :: Deplacement2
      real*8, dimension(:,:), allocatable :: Position_atom2
      real*8, dimension(:),   allocatable :: Maille

      INTEGER, parameter:: dim=90
      INTEGER, parameter:: dim2=2025
      INTEGER, parameter:: dim3=4
      INTEGER, parameter:: dim4=30


      integer l,m,o,p
      real*8 bid,a
      real*8 a0
      real*8 ha
      a0=0.52917720859e-10
      ha=27.2113838668
      allocate(Maille(3))
      allocate(Bid2(dim))
      allocate(Bid2(dim*dim))
      allocate(Deplacement(dim))
      allocate(Tmp(dim2,dim3))
      allocate(Hessien(dim,dim))
      allocate(Hessien2(dim,dim))
      allocate(Deplacement2(dim,dim))
      allocate(Inverse(dim,dim))
      allocate(Verification(dim,dim))
      allocate(Born(dim,3))
      allocate(Electric_field(3,1))
      allocate(Position_atom(dim4,3))
      allocate(Position_atom2(dim4,3))
      allocate(Force(dim))

      PRINT*, 'boucle de contrôle'


      do i=1,dim
         Bid2(i)=0
         do j=1,dim
            Hessien(i,j)=0
            Inverse(i,j)=0
            Verification(i,j)=0
         enddo
      enddo
      
!-----------Lecture du Hessien------------------------------------!
      OPEN(8,FILE='HESSIEN.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
      DO i=1,dim2
         READ(8,*) (Tmp(i,j),j=1,dim3,1)
      ENDDO
      l=1
      m=1
      Do i=1,dim2
         Do j=1,dim3
            o=m
            p=l
            Hessien(l,m)=Tmp(i,j)
            Hessien2(l,m)=Tmp(i,j)
            m=o+1
            if (m.gt.dim)   then 
               l=p+1
               m=1
            else
            endif
         enddo
      enddo
!-----------fin-de-lecture du Hessien-------------------------------!

      
          
      OPEN(8,FILE='BORN.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
      DO i=1,dim
         READ(8,*) (BORN(i,j),j=1,3,1)
      ENDDO
      
      OPEN(8,FILE='POSITION.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
      DO i=1,30
         READ(8,*) (Position_atom(i,j),j=1,3,1)
      ENDDO

      OPEN(8,FILE='MAILLE.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
      READ(8,*) (Maille(i),i=1,3,1)
      
      OPEN(8,FILE='ELECTRIC.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
      READ(8,*) (Electric_field(i,1),i=1,3,1)

      PRINT*, ' fin de lecture des différents fichiers'
!--------------------------------------------------------------------!
!--------------------------------------------------------------------!
!--------------------------------------------------------------------!
      PRINT*, ' Inversion de la matrice hessienne'

      CALL DCOPY(dim*dim,Hessien,1,Inverse,1)
      CALL DGETRF(dim,dim,Inverse,dim,Bid2,INFO)
      CALL DGETRI(dim,Inverse,dim,Bid2,Tmp2,dim*dim,INFO)


      bid=0
      a=0
      DO i=1,dim
         DO j=1,dim
            DO n=1,dim
               bid=Hessien2(i,n)*Inverse(n,j)
               Verification(i,j)=a+bid
               a=Verification(i,j)
            bid=0
            ENDDO
            a=0
         ENDDO
      ENDDO

      Do i=1,dim
         do j=1,dim
            if(Verification(i,j).lt.10E-5) then
               Verification(i,j)=0
            else
               if(Verification(i,j).gt.10E-5) then 
                  print*, 'valueur non nulle', Verification(i,j), 
     &                    'indice ', i, j
               endif
            endif
         enddo
      enddo

c      do i=1,dim
c         write(*,*) (Verification(i,j),j=1,dim,1)
c         print*, ' '
c      enddo

      bid=0
      Do i=1,dim
         Do j=1,3
            Force(i)=bid+Born(i,j)*Electric_field(j,1)
            bid=Force(i)
         enddo
         bid=0
      enddo
                  

!--------------------------------------------------------------------!
!---------------calculs des différents déplacements------------------!
!--------------------------------------------------------------------!
      bid=0
      Do i=1,dim
         Do j=1,dim
         Deplacement(i)=bid+Inverse(i,j)*Force(j)
         bid=Deplacement(i)
         enddo
         bid=0
      enddo


      j=1
      k=1
      o=1
      p=1
      Do i=1,dim
         o=j
         p=k
         Deplacement2(k,j)=Deplacement(i)*a0*0.52917720859/27.2113838668
         j=o+1
         if(j.gt.3) then
            k=p+1
            j=1
         else
         endif
      enddo
     
  
      Do i=1,dim4
         do j=1,3
            Position_atom2(i,j)=Position_atom(i,j)
     &      +(Deplacement2(i,j)/(Maille(j)))
         enddo
      Enddo


      OPEN(8,FILE='SORTIE.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
      Do i=1,dim4
            WRITE(8,*), (Position_atom2(i,j), j=1,3,1)
      Enddo
      
       
      END PROGRAM
