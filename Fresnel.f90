module Fresnel
contains
   subroutine FCS(C, S, t)
      implicit none

      real(8), intent(in) :: t
      real(8), intent(out) :: C, S

      call frenel(t, s, c)
   end subroutine FCS
   !subroutine FCS(C, S, t)
   !   implicit none
   !
   !   real(8), intent(in) :: t
   !   real(8) :: C, S
   !
   !   real(8) :: pi = dacos(-1.0D0), x
   !
   !   x = pi*t*t/2.0
   !
   !   call CS(C, S, x)
   !
   !   if (t < 0) then
   !      C = -C
   !      S = -S
   !   end if
   !end subroutine FCS
end module Fresnel
