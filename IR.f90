module IR_mod
contains
   complex(8) function IR(u, step, DT)
      implicit none

      integer(4), intent(in) :: step
      real(8), intent(in) :: DT
      complex(8), dimension(step), intent(in) :: u

      integer(4) j
      real(8) :: coeff_4_d_3_m_SQRDT, SQR2M2 = 2.828427124746190, SQR2D2 = 0.707106781186548

      coeff_4_d_3_m_SQRDT = 4.0D0/3.0D0*dsqrt(DT)

      if (step > 2) then
         IR = coeff_4_d_3_m_SQRDT*(u(1)*((step - 1.0D0)**(1.5) - (step - 1.5D0)*dsqrt(dble(step))) + u(step)*(SQR2M2 - 2.5D0))
         do j = 1, step - 2
          IR = IR + coeff_4_d_3_m_SQRDT*(u(j + 1)*((step - j - 1.0D0)**(1.5) - 2.0D0*(step - j)**(1.5) + (step - j + 1.0D0)**(1.5)))
         end do
      elseif (step == 2) then
         IR = coeff_4_d_3_m_SQRDT*(u(1)*(1 - SQR2D2) + u(2)*(SQR2M2 - 2.5D0))
      else
         IR = 0; 
      end if
   end function IR
end module IR_mod
