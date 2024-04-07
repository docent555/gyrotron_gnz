module integral_A4_Zend_mod
   use Fresnel
contains
   subroutine integral_A4_Zend(I1, I2, I3, I4, f0, f1, dfdz0, dfdz1, Zout, DT)
      implicit none

      real(8), intent(in) :: f0, f1, dfdz0, dfdz1, Zout, DT
      complex(8), intent(out) :: I1, I2, I3, I4

      real(8) :: pi = dacos(-1.0D0)
      real(8) :: Sqrt2PI
      real(8) :: SqrtDT
      real(8) :: SqrtPid2
      real(8) :: Zoutm2dpidh
      real(8) :: isqrZoutdh
      real(8) :: C_Zoutm2dpidh, S_Zoutm2dpidh, Cos_SqrZoutdDT, Sin_SqrZoutdDT

      Sqrt2PI = dsqrt(2*pi)
      SqrtDT = dsqrt(DT)
      SqrtPid2 = dsqrt(Pi/2.)
      Zoutm2dpidh = Zout*dsqrt(2/(pi*DT))
      isqrZoutdh = -(0, 1)*Zout*Zout/DT

      Cos_SqrZoutdDT = dcos(Zout**2/DT)
      Sin_SqrZoutdDT = dsin(Zout**2/DT)

      call FCS(C_Zoutm2dpidh, S_Zoutm2dpidh, Zoutm2dpidh)

      I1 = (SqrtPid2*(1 - 2*C_Zoutm2dpidh))/Zout + ((0, 1)*SqrtPid2*(-1 + 2*S_Zoutm2dpidh))/Zout

      I2 = (-1, -1)*Sqrt2PI*Zout + 2*SqrtDT*Cos_SqrZoutdDT + (0, 2)*Sqrt2PI*Zout*C_Zoutm2dpidh + 2*Sqrt2PI*Zout*S_Zoutm2dpidh - &
           (0, 2)*SqrtDT*Sin_SqrZoutdDT

   I3 = -0.5D0*(2*Pi*Zout - (1, -1)*SqrtDT*Sqrt2PI*Cos_SqrZoutdDT - (2, 2)*Pi*Zout*C_Zoutm2dpidh - (2, -2)*Pi*Zout*S_Zoutm2dpidh + &
                   (1, 1)*SqrtDT*Sqrt2PI*Sin_SqrZoutdDT)/Zout

     I4 = (DT*Pi + (0, 2)*Pi*Zout**2 - (1, 1)*SqrtDT*Sqrt2PI*Zout*Cos_SqrZoutdDT - (1, 1)*Pi*(DT + (0, 2)*Zout**2)*C_Zoutm2dpidh - &
            (1, -1)*Pi*(DT + (0, 2)*Zout**2)*S_Zoutm2dpidh - (1, -1)*SqrtDT*Sqrt2PI*Zout*Sin_SqrZoutdDT)/2.

   end subroutine integral_A4_Zend
end module integral_A4_Zend_mod
