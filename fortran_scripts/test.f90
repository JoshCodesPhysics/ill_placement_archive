Program test
  Implicit none
  ! E (u.a.) = alpha^2 m c^2 / (e a0) E[V/m]
  ! alpha =  7.297 352 5693 x 10-3
  ! c     =  299 792 458 m s-1
  ! m     =  9.109 383 7015 x 10-31 kg
  ! e     =  1.602 176 634  x 10-19 C 
  ! -e/m  = -1.758 820 010 76 x 1011 C kg-1
  ! a0    =  5.291 772 109 03 x 10-11 m 
  Real*8, parameter :: &
       alpha =  7.2973525693D0, &
       c     =  2.99792458D0, &
       m     =  9.1093837015D0, &
       e     =  1.602176634D0, &
       a0    =  5.29177210903D0
  Real*8 :: kVm
  character*80 ligne
  real*8 :: tmp
  Integer :: i, iatom, j, n, sh, jj, err
  real*8, dimension(30) :: Efield
  real*8, dimension(3,3,10) :: bch
  real*8, dimension(30,30)  :: H , Q, V
  real*8, dimension(30*31/2):: HL
  real*8, dimension(30,1)   :: QE,d
  integer, dimension(30)    :: ipiv
  real*8, dimension(90)     :: work 
  real*8, dimension(30)     :: Ener

  ! write(6,*) alpha
  kVm = alpha*alpha * m * c*c / (e*a0) * 1.d6
  Write(6,'("E(u.a.) = ",E20.12," x E(kV/m)")') kVm

  open(1,file="/home/joshhorswill10/Documents/git_new/joshua_3/Examples/disp_solve_examples/ht.frequence.B1PW_PtBs.loto.out",&
          form="formatted")
  
!!$ Lecture du champ electrique
  print *, "Ex, Ey, Ez en kV/m?"
  Efield(:) = 0.d0
  read(5,*) Efield(1:3)
  write (6,*) " Champ electrique en kV/m= ",  Efield(1:3)
  Efield(1:3) = Efield(1:3) / kVm
  write (6,*) " Champ electrique en ua  = ",  Efield(1:3)

  do i=2,10
     sh = (i-1)*3
      Efield(1+sh:3+sh) = Efield(1:3)
  end do

  write(6,*)
  write(6,*) " Champ electrique en u.a."
  call wrtv(Efield,30)

!!$ Lecture des charges de Born
  bch(:,:,:) = 0.d0
  do i=1,1243
     read(1,*) 
  end do

  do iatom = 1, 10
     read(1,'(80a)') ligne
     write(6,*) ligne
     read(1,*) 
     read(1,*)
     read(1,9001) j, bch(1,1:3,iatom)
     read(1,9001) j, bch(2,1:3,iatom)
     read(1,9001) j, bch(3,1:3,iatom)
     read(1,*) 
     ! write(6,9001) 1,bch(1,1:3,iatom)
     ! write(6,9001) 2,bch(2,1:3,iatom)
     ! write(6,9001) 3,bch(3,1:3,iatom)
     ! write(6,*)
  end do

!!$ mise à 0 des 0
  do iatom = 1,10
     do j=1,3
        do i=1,3
           if (abs(bch(i,j,iatom)) .lt. 10.d-12) bch(i,j,iatom) = 0.d0
        end do
     end do
  end do

!!$ Q
  Q(:,:) = 0.d0
  do iatom = 1,10
     sh = (iatom-1)*3
     do j=1,3
        do i=1,3
            Q(i+sh,j+sh) = bch(i,j,iatom) 
        end do
     end do
  end do

  ! Write(6,*)
  ! Write(6,*) " Charges de Born"
  ! call wrtm(q,30,30)
  
!!$ QE
  QE(:,:) = 0.d0
  do j=1,30
     QE(:,j) = Q(:,j)*Efield(j)
  end do

  Write(6,*)
  Write(6,*) " Q.E"
  call wrtm(QE,30,30)

  
!!$ Lecture du Hessien
  H(:,:) = 0.d0
  read(1,*) 
  read(1,'(80a)') ligne
  Write(6,*)
  write(6,*) ligne

  do i=1,3
     read(1,*) 
  end do

  do n=1,6
     sh = (n-1)*5
     read (1,*)
     do j=1,30
        read (1,9002)  jj, H(j,1+sh:5+sh)
        if (jj.ne.j) stop " erreur "
     end do
     read (1,*)
  end do

  Write(6,*)
  Write(6,*) " Hessien "
  call wrtm(H,30,30)

  n=0
  do j=1,30
     do i= j,30
        n=n+1
        HL(n) = H(i,j)
     !   write(6,*) n,i,j,HL(n)
     end do
  end do
  
 ! Write(6,*)
 ! Write(6,*) " Hessien triangle bas"
  
  
!!$ -Hd + QE = 0 Hd = QE
  ! F = Q E
  ipiv(:) = 0
  work(:) = 0.d0
  V(:,:)  = H(:,:)
  Ener(:) = 0.d0
  ! call DSYSV('L',30,1,H,30,ipiv,QE,30,work,60, err)
  ! call DSPSV('L',30,1,HL,ipiv,QE,30, err)
  ! call DPOSV('L',30,1,H,30,QE,30, err)
  call DSYEV('V','L',30,V,30,Ener,work,90,err)

  if (err.ne.0) then
     write(6,*)
     write(6,*) " pb dans inversion systeme lineaire"
     write(6,*)
  end if

  write(6,*)
  write(6,*) " vecteurs pp "
  call wrtm(V,30,30)

  write(6,*)
  write(6,*) " valeurs pp "
  call wrtv(Ener,30)

  ! verif1
  Q(:,:) = 0.d0
  do i=1,30
     do j = 1,30
        do n=1,30
           Q(i,j) = Q(i,j) + V(i,n)  * V(j,n) 
        end do
     end do
  end do
  tmp = 0.d0
  do i=1,30
     if (abs(Q(i,i)-1).gt.tmp) tmp = abs(Q(i,i)-1)
  end do
  write(6,*) "erreur sur diag la + grande =", tmp
  tmp = 0.d0
  do i=1,30
     do j = 1,30
        if (j.eq.i) cycle
        if (abs(Q(i,j)).gt.tmp) tmp = abs(Q(i,j))
     end do
  end do
  write(6,*) "erreur sur extra-diag la + grande =", tmp

  ! verif2
  Q(:,:) = 0.d0
  do i=1,30
     do j = 1,30
        do n=1,30
           Q(i,j) = Q(i,j) + V(i,n)  * Ener(n) * V(j,n) 
        end do
     end do
  end do
  tmp = 0.d0
  do i=1,30
     do j = 1,30
        if (abs(H(i,j)-Q(i,j)).gt.tmp) tmp = abs(H(i,j)-Q(i,j))
     end do
  end do
  write(6,*)
  write(6,*) " Verif hessien "
  call wrtm(Q,30,30)
  write(6,*) "erreur sur diagonalisation =", tmp
  

  Q(:,:) = 0.d0
  do i=1,30
     do j = 1,30
        do n=1,30
           Q(i,j) = Q(i,j) + V(i,n)  * (1.d0/Ener(n)) * V(j,n) 
        end do
     end do
  end do

  write(6,*)
  write(6,*) " H^[-1] "
  call wrtm(Q,30,30)

  
  d(:,:) = 0.d0  
  do i=1,30
     do j = 1,30
        d(i,1) = d(i,1) + Q(i,j) *  QE(j,1)
     end do
  end do
  
  write(6,*)
  write(6,*) " d "
  call wrtm(d,30,1)
  
  write(6,*)
  write(6,*) " Deplacement"
  do iatom = 1,10
     sh = (iatom-1)*3
     write(6,9003) iatom, d(1+sh : 3+sh,1)
  end do

  
9001 format(3x,i1,3x,3ES12.4)
9002 format (1x,i3,3x,5(E13.6,1x))
9003 format(3x,i3,3x,3(E20.12,2x))

contains
  subroutine wrtm(a,n,p)
    implicit none
    integer, intent(in) :: n,p
    real*8, dimension(n,p), intent(in)  :: a
    integer :: q, r, sh, k, j

    q = p/5
    r = p-q*5

    if (q.eq.0) goto 10
    do k=1,q
     sh = (k-1)*5
     write (6,9001) (j, j=1+sh,5+sh)
     do j=1,n
        write (6,9002) j, a(j,1+sh:5+sh)
     end do
     write (6,*)
  end do

10 if (r.eq. 0) return
  sh = q*5
  write (6,*)  (j, j=1+sh,p)
  do j=1,n
     write (6,9002) j, a(j,1+sh:p)
  end do
  write (6,*)

9001 format (7x,5(5x,i4,5x))
9002 format (1x,i3,3x,5(E13.6,1x))
end subroutine wrtm

  subroutine wrtv(a,n)
    implicit none
    integer, intent(in) :: n
    real*8, dimension(n), intent(in)  :: a
    integer :: j

     do j=1,n
        write (6,9002) j, a(j)
     end do
     write (6,*)

9002 format (1x,i3,3x,5(E13.6,1x))
   end subroutine wrtv

end Program test
