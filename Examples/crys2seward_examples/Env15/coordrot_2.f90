!!$***********************************************************************
!!$
subroutine coord_rot_roti(R,MX,a2,b2,c2,a,b,c,alpha,beta,gamma,iprint)
!!$
!!$***********************************************************************
!!$
!!$          Routine de calcul de la rotation a faire
!!$          pour obtenir les coordonnées réelles (en A) des atomes
!!$          dans la base définie par (a2;b2;c2) re-orthonormalisée :
!!$
!!$            X= a2  ;  Y= b2-(X.b2)X  ;  Z= c2-(X.c2)X-(Y.c2)Y    +normalisation
!!$
!!$          où :  (.)= prod.scal. ,
!!$          ce qui donne (avec a2,b2,c2 renormalisés) :
!!$
!!$            X= a2
!!$            Y= ( b2 - cos(ab) a2 ) / sin(ab)
!!$            Z= ( sin²(ab) c2 +(cos(ab)cos(bc)-cos(ac)) a2 +(cos(ab)cos(ac)-cos(bc)) b2 ) / Z0
!!$
!!$          avec  Z0= sin(ab) sqrt(1 -cos²(ab) -cos²(ac) -cos²(bc) +2 cos(ab)cos(ac)cos(bc))
!!$          et ab= angle entre a2 et b2 
!!$
!!$          on obtient la matrice de passage :    MX
!!$
!!$          que l'on inverse                 :    R= (MX)-1
!!$
!!$          on peut verifier la rotation avec  (iprint>=1) :
!!$
!!$            x=  x2  +  y2  cos(ab)  +  z2  cos(ac)
!!$            y=         y2  sin(ab)  +  z2  (cos(bc)-cos(ab)cos(ac))/sin(ab)
!!$            z=                         z2  Z0 / sin(ab)²
!!$
!!$          et en calculant la matrice de passage de B->B2
!!$
!!$----- Déclarations  
!!$
       implicit none
!!$
!!$ Variables globales
       real*8, parameter :: pi=3.14159265358979d0 
!!$ variables de définition du cristal
       integer, intent(in) :: iprint
       real*8, dimension(3,3) :: R, MX
       real*8, dimension(3), intent(in) :: a2, b2, c2
       real*8, intent(in) :: a, b, c, alpha, beta, gamma  ! maille du cristal
       real*8 :: talpha, tbeta, tgamma
       logical :: yrot
       real*8 :: alpha3, beta3, gamma3, a3,b3,c3, ta,tb,tc, tmp
       real*8, dimension(3,3) :: Rot, RotI, Mtmp
       real*8 :: lim6=1.d-6 
       integer :: i
!!$
!!$---------------------------------------------------------------------------------------
!!$
!!$ Impressions 
!!$
      write(6,"(/' axes  a2,b2,c2 dans la base a,b,c  :')") 
      write(6,"(' a2=', 3(4x,ES12.5))") (a2(i),i=1,3)
      write(6,"(' b2=', 3(4x,ES12.5))") (b2(i),i=1,3)
      write(6,"(' c2=', 3(4x,ES12.5))") (c2(i),i=1,3)
!!$
!!$ Initialisations    
!!$
       R(:,:)  = 0.d0
       MX(:,:) = 0.d0
       yrot=.false.
       talpha = alpha /180.d0*pi
       tbeta  = beta  /180.d0*pi
       tgamma = gamma /180.d0*pi
!!$
!!$  Tests
!!$
       if       (a2(1)==1.d0.and.a2(2)==0.d0.and.a2(3)==0.d0 &
            .and.b2(1)==0.d0.and.b2(2)==1.d0.and.b2(3)==0.d0 &
            .and.c2(1)==0.d0.and.c2(2)==0.d0.and.c2(3)==1.d0) &
            then 
          a3     = a
          b3     = b
          c3     = c
          alpha3 = alpha
          beta3  = beta
          gamma3 = gamma
          ta  = (alpha3 -90.d0)/180.d0*pi
          tb  = (beta3 - 90.d0)/180.d0*pi
          tc  = (gamma3 -90.d0)/180.d0*pi
       else
          yrot=.true.
          rot(:,1) = a2(:)
          rot(:,2) = b2(:)
          rot(:,3) = c2(:)
          a3 = (a2(1)*a)*(a2(1)*a) + (a2(2)*b)*(a2(2)*b) + (a2(3)*c)*(a2(3)*c) &
              + 2.d0 * (a2(1)*a*a2(2)*b*dcos(tgamma) + &
                        a2(1)*a*a2(3)*c*dcos(tbeta) + &
                        a2(2)*b*a2(3)*c*dcos(talpha))
          a3 = dsqrt(a3)
          b3 = (b2(1)*a)*(b2(1)*a) + (b2(2)*b)*(b2(2)*b) + (b2(3)*c)*(b2(3)*c) &
              + 2.d0 * (b2(1)*a*b2(2)*b*dcos(tgamma) + &
                        b2(1)*a*b2(3)*c*dcos(tbeta) + &
                        b2(2)*b*b2(3)*c*dcos(talpha))
          b3 = dsqrt(b3)
          c3 = (c2(1)*a)*(c2(1)*a) + (c2(2)*b)*(c2(2)*b) + (c2(3)*c)*(c2(3)*c) &
              + 2.d0 * (c2(1)*a*c2(2)*b*dcos(tgamma) + &
                        c2(1)*a*c2(3)*c*dcos(tbeta) + &
                        c2(2)*b*c2(3)*c*dcos(talpha))
          c3 = dsqrt(c3)
          alpha3 =  b2(1)*a*c2(1)*a + b2(2)*b*c2(2)*b + b2(3)*c*c2(3)*c &
               + (b2(1)*c2(2) + c2(1)*b2(2))* a*b * dcos(tgamma) &
               + (b2(1)*c2(3) + c2(1)*b2(3))* a*c * dcos(tbeta)  &
               + (b2(2)*c2(3) + c2(2)*b2(3))* b*c * dcos(talpha)
          alpha3 =  alpha3 / (b3*c3)
          alpha3 =  dacos(alpha3)
          ta     =  alpha3 - pi/2.d0
          alpha3 =  alpha3*180.d0/pi
          beta3 =  a2(1)*a*c2(1)*a + a2(2)*b*c2(2)*b + a2(3)*c*c2(3)*c &
               + (a2(1)*c2(2) + c2(1)*a2(2))* a*b * dcos(tgamma) &
               + (a2(1)*c2(3) + c2(1)*a2(3))* a*c * dcos(tbeta)  &
               + (a2(2)*c2(3) + c2(2)*a2(3))* b*c * dcos(talpha)
          beta3 =  beta3 / (a3*c3)
          beta3 =  dacos(beta3)
          tb    =  beta3 - pi/2.d0
          beta3 =  beta3*180.d0/pi
          gamma3 =  a2(1)*a*b2(1)*a + a2(2)*b*b2(2)*b + a2(3)*c*b2(3)*c &
               + (a2(1)*b2(2) + b2(1)*a2(2))* a*b * dcos(tgamma) &
               + (a2(1)*b2(3) + b2(1)*a2(3))* a*c * dcos(tbeta)  &
               + (a2(2)*b2(3) + b2(2)*a2(3))* b*c * dcos(talpha)
          gamma3 =  gamma3 / (a3*b3)
          gamma3 =  dacos(gamma3)
          tc     =  gamma3 - pi/2.d0
          gamma3 =  gamma3*180.d0/pi
       end if

       write(6,*)
       write(6,"('a2, b2, c2, alpha2, beta2, gamma2')")
       write(6,"(3(4x,F10.5), 7x, 3(4x,F10.5))") a3,b3,c3, alpha3, beta3, gamma3

!!$ Construction des matrices de rotation
10     continue

      if (abs(alpha3-90.0d0).lt.lim6 .and. abs(beta3-90.0d0).lt.lim6 &
           .and. abs(gamma3-90.0d0).lt.lim6) then 
         ! base orthonormée
         R(1,1) = a3
         R(2,2) = b3
         R(3,3) = c3
         MX(1,1) = 1.d0/a3
         MX(2,2) = 1.d0/b3
         MX(3,3) = 1.d0/c3
!        do i= 1,3 
!           R(i,i)  = 1.d0
!           MX(i,i) = 1.d0
!        end do
      else if (abs(alpha3-90.0d0).lt.lim6 .and. abs(beta3-90.0d0).ge.lim6 &
           .and. abs(gamma3-90.0d0).lt.lim6) then
         ! beta3 /= 90
         R(1,1) =  a3
         R(1,3) = -c3*dsin(tb)
         R(2,2) =  b3
         R(3,3) =  c3*dcos(tb)
         MX(1,1) =  1.d0/a3
         MX(1,3) =  1d0/a3 * dtan(tb)
         MX(2,2) =  1.d0/b3
         MX(3,3) =  1.d0/(c3*dcos(tb))
      else if (abs(alpha3-90.0d0).lt.lim6 .and. abs(beta3-90.0d0).lt.lim6 &
           .and. abs(gamma3-90.0d0).ge.lim6) then
         ! gamma3 /= 90
         R(1,1) =  a3
         R(1,2) = -b3*dsin(tc)
         R(2,2) =  b3*dcos(tc)
         R(3,3) =  c3
         MX(1,1) =  1.d0/a3
         MX(1,2) =  1d0/a3 * dtan(tc)
         MX(2,2) =  1.d0/(b3*dcos(tc))
         MX(3,3) =  1.d0/c3
      else 
         ! cas general
         tmp = 1.d0 - dsin(ta)**2 - dsin(tb)**2 - dsin(tc)**2 &
              - 2.d0*dsin(ta)*dsin(tb)*dsin(tc)
         tmp = dsqrt(tmp)
         R(1,1) =  a3
         R(1,2) = -b3*dsin(tc)
         R(1,3) = -c3*dsin(tb)
         R(2,2) =  b3*dcos(tc)
         R(2,3) = -c3*(dsin(tb)*dsin(tc)+dsin(ta))/dcos(tc)
         R(3,3) =  c3*tmp/dcos(tc)
         MX(1,1) =  1.d0/a3
         MX(1,2) =  1.d0/a3 * dtan(tc)
         MX(1,3) =  1.d0/a3 * (dsin(ta)*dsin(tc)+dsin(tb)) / (dcos(tc)*tmp)
         MX(2,2) =  1.d0/(b3*dcos(tc))
         MX(2,3) =  1.d0/b3 * (dsin(tb)*dsin(tc)+dsin(ta)) / (dcos(tc)*tmp)
         MX(1,3) =  1.d0/c3 * dcos(tc)/tmp
      end if


      write(6,"(/' ---------------------- rotation de la base ---------------------------')")
      if (yrot)  then
         call inv_mat(RotI,Rot,3)
         Mtmp(:,:) = R(:,:)
         R(:,:)    = matmul(Mtmp,RotI)
         Mtmp(:,:) = MX(:,:)
         MX(:,:)   = matmul(Rot,Mtmp)
      end if
    
      write(6,"(/' axes  a,b,c dans la base X,Y,Z  :')") 
      write(6,"(3(4x,ES12.5))") (R(i,:),i=1,3)

      
end subroutine coord_rot_roti


!!$
!!$
!!$
!!$

!!$***********************************************************************
!!$
subroutine coord_rot(R,a2,b2,c2,a,b,c,alpha,beta,gamma,iprint)
!!$
!!$***********************************************************************
!!$
!!$          Routine de calcul de la rotation a faire
!!$          pour obtenir les coordonnées réelles (en A) des atomes
!!$          dans la base définie par (a2;b2;c2) re-orthonormalisée :
!!$
!!$            X= a2  ;  Y= b2-(X.b2)X  ;  Z= c2-(X.c2)X-(Y.c2)Y    +normalisation
!!$
!!$          où :  (.)= prod.scal. ,
!!$          ce qui donne (avec a2,b2,c2 renormailsés) :
!!$
!!$            X= a2
!!$            Y= ( b2 - cos(ab) a2 ) / sin(ab)
!!$            Z= ( sin²(ab) c2 +(cos(ab)cos(bc)-cos(ac)) a2 +(cos(ab)cos(ac)-cos(bc)) b2 ) / Z0
!!$
!!$          avec  Z0= sin(ab) sqrt(1 -cos²(ab) -cos²(ac) -cos²(bc) +2 cos(ab)cos(ac)cos(bc))
!!$          et ab= angle entre a2 et b2 
!!$
!!$          on obtient la matrice de passage :    MX
!!$
!!$          que l'on inverse                 :    R= (MX)-1
!!$
!!$          on peut verifier la rotation avec  (iprint>=1) :
!!$
!!$            x=  x2  +  y2  cos(ab)  +  z2  cos(ac)
!!$            y=         y2  sin(ab)  +  z2  (cos(bc)-cos(ab)cos(ac))/sin(ab)
!!$            z=                         z2  Z0 / sin(ab)²
!!$
!!$          et en calculant la matrice de passage de B->B2
!!$
!!$----- Déclarations  
!!$
       implicit none
!!$
!!$ Variables globales
       real*8, parameter :: pi=3.14159265358979d0 
!!$ variables de définition du cristal
       integer, intent(in) :: iprint
       real*8, dimension(3,3) :: R
       real*8, dimension(3), intent(in) :: a2, b2, c2
       real*8, intent(in) :: a, b, c, alpha, beta, gamma  ! maille du cristal
       real*8, dimension(3,3) :: ms, ms2, MB2, MX, msX, RX, R2, Di, Rv
       real*8 :: norm, sab, Z0
       integer :: i, j, k, l
!!$
!!$
!!$       write(6,"(/' ---------------------- rotation de la base ---------------------------')")
!!$       
!!$----- calc matrice de recouvrement
       ms(1,1)= a*a  ;  ms(2,2)= b*b  ;  ms(3,3)= c*c
       ms(1,2)= a*b* dcos(gamma*pi/180.D0)
       ms(1,3)= a*c* dcos(beta *pi/180.D0)
       ms(2,3)= b*c* dcos(alpha*pi/180.D0)
       ms(2,1)= ms(1,2)  ;  ms(3,1)= ms(1,3)  ;  ms(3,2)= ms(2,3)
       if ( iprint>=1 ) then
          write(6,"(/' matrice de recouvrement de a,b,c :')") 
          write(6,"(3(4x,ES12.5))") (ms(i,:),i=1,3)
       end if
!!$
!!$----- normalisation de a2, b2, c2; et recouvrement
!!$
       ! normalisation
       MB2(:,1)= a2(:)
       MB2(:,2)= b2(:)
       MB2(:,3)= c2(:)
       if ( iprint>=1 ) then
          write(6,"(/' axes a2,b2,c2 dans la base a,b,c :')") 
          write(6,"(3(4x,ES12.5))") (MB2(i,:),i=1,3)
       end if
       do i= 1, 3
          norm= dot_product(MB2(:,i),matmul(ms,MB2(:,i)))
          MB2(:,i)= MB2(:,i)/sqrt(norm)
       end do
       write(6,"(/' axes a2,b2,c2 renormes dans la base a,b,c  :')") 
       write(6,"(3(4x,ES12.5))")  (MB2(i,:),i=1,3)
       ! recouvrement 
       ms2= matmul(transpose(MB2),matmul(ms,MB2))
       if ( iprint>=1 ) then
          write(6,"(/' matrice de recouvrement de a2,b2,c2 :')") 
          write(6,"(3(4x,ES12.5))") (ms2(i,:),i=1,3)
       end if
       ! sinus (a2,b2)
       sab= sqrt(1.D0-ms2(1,2)*ms2(1,2))
!!$
!!$----- Z0
!!$
       Z0= sab*sqrt(1-ms2(1,2)*ms2(1,2)-ms2(1,3)*ms2(1,3)-ms2(2,3)*ms2(2,3) &
            +2*ms2(1,2)*ms2(1,3)*ms2(2,3))
!!$
!!$----- base X,Y,Z
!!$
       ! les vecteurs
       MX(:,1)= MB2(:,1)
       MX(:,2)= (MB2(:,2)-ms2(1,2)*MB2(:,1)) / sab
       MX(:,3)= ( sab*sab*MB2(:,3) +(ms2(1,2)*ms2(2,3)-ms2(1,3))*MB2(:,1) &
            +(ms2(1,2)*ms2(1,3)-ms2(2,3))*MB2(:,2) ) / Z0
       write(6,"(/' axes X,Y,Z dans la base a,b,c  :')") 
       write(6,"(3(4x,ES12.5))")  (MX(i,:),i=1,3)
       ! recouvrement
       if ( iprint>=1 ) then
          msX= matmul(transpose(MX),matmul(ms,MX))
          write(6,"(/' matrice de recouvrement de X,Y,Z :')") 
          write(6,"(3(4x,ES12.5))") (msX(i,:),i=1,3)
       end if
!!$
!!$----- Matrice de rotation
!!$
       call inv_mat(R,MX,3)
       write(6,"(/' axes  a,b,c dans la base X,Y,Z  :')") 
       write(6,"(3(4x,ES12.5))") (R(i,:),i=1,3)
!!$       
!!$----- verification de la rotation
!!$
       if ( iprint>=1 ) then
          ! matrice de passage B->B2
          call inv_mat(R2,MB2,3)
          write(6,"(/'  a,b,c dans a2,b2,c2  :')") 
          write(6,"(3(4x,ES12.5))")  (R2(i,:),i=1,3)
          ! matrice de passage B2->X
          RX(1,:)= (/ 1.D0, ms2(1,2), ms2(1,3) /)
          RX(2,:)= (/ 0.D0, sab, (ms2(2,3)-ms2(1,2)*ms2(1,3))/sab /)
          RX(3,:)= (/ 0.D0, 0.D0, Z0/sab/sab /)
          write(6,"(/'  a2,b2,c2 dans X,Y,Z  :')") 
          write(6,"(3(4x,ES12.5))")  (RX(i,:),i=1,3)
          ! matrice de passage tot.
          Rv= matmul(RX,R2)
          write(6,"(/' axes  a,b,c dans la base X,Y,Z (verif) :')") 
          write(6,"(3(4x,ES12.5))") (Rv(i,:),i=1,3)
          ! difference
          Di=R-Rv
          write(6,"(/' difference   :')") 
          write(6,"(3(4x,ES12.5))")  (Di(i,:),i=1,3)
       end if
!!$       
!!$----- fin
!!$
end subroutine coord_rot
!!$
!!$
!!$
!!$

