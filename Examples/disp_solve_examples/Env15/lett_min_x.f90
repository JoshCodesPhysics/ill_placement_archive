! version 1
!
!***********************************************************************
!
!      dans une chaine de caracteres : 
!
!      transforme les minuscules en majuscules
!
!***********************************************************************

subroutine lett_min(a)
       implicit none
       character(len=*) :: a
       integer :: l, i, ic

       l=len_trim(a)
       do i= 1, l
          ic= ichar(a(i:i))
          if ( ic>=64 .and. ic<=90 )   a(i:i)= char(ic+32)
       end do
end subroutine lett_min
