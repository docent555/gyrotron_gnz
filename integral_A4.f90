module integral_A4_mod
   use Fresnel
contains
   subroutine integral_A4(Int, z, f0, f1, dfdz0, dfdz1, Zout, DT)
      implicit none

      real(8), intent(in) :: z, f0, f1, dfdz0, dfdz1, Zout, DT
      complex(8), intent(out) :: Int

      real(8) :: pi = dacos(-1.0D0)
      complex(8) :: I1, I2, I3, I4, C1, C2, C3, C4
      real(8) :: ZmZout
      real(8) :: ZpZout
      real(8) :: SqrZmZoutd4H
      real(8) :: SqrZpZoutd4H
      real(8) :: SqrZmZoutdH
      real(8) :: SqrZpZoutdH
      real(8) :: AbsZmZoutdSqrtH2PI
      real(8) :: AbsZpZoutdSqrtH2PI
      real(8) :: Sqrt2PI
      real(8) :: SqrtH
      real(8) :: SqrtPI
      complex(8) :: SqrtPII
      real(8) :: Sqrt2m2
      real(8) :: AbsZmZout
      real(8) :: AbsZpZout
      real(8) :: cos_SqrZmZoutd4H
      real(8) :: cos_SqrZpZoutd4H
      real(8) :: sin_SqrZmZoutd4H
      real(8) :: sin_SqrZpZoutd4H
      real(8) :: Sqrt2
      real(8) :: SqrtSqrZmZoutdH
      real(8) :: SqrtSqrZpZoutdH
      real(8) :: C_AbsZmZoutdSqrtH2PI
      real(8) :: S_AbsZmZoutdSqrtH2PI
      real(8) :: C_AbsZpZoutdSqrtH2PI
      real(8) :: S_AbsZpZoutdSqrtH2PI

      ZmZout = z - Zout
      ZpZout = z + Zout
      SqrZmZoutd4H = ZmZout*ZmZout/(4*DT)
      SqrZpZoutd4H = ZpZout*ZpZout/(4*DT)
      SqrZmZoutdH = ZmZout*ZmZout/DT
      SqrZpZoutdH = ZpZout*ZpZout/DT
      AbsZmZoutdSqrtH2PI = dabs(ZmZout)/dsqrt(DT*2.0D0*pi)
      AbsZpZoutdSqrtH2PI = dabs(ZpZout)/dsqrt(DT*2.0D0*pi)
      Sqrt2PI = dsqrt(2.0D0*pi)
      SqrtH = dsqrt(DT)
      SqrtPI = dsqrt(pi)
      SqrtPII = cdsqrt(pi*(0, 1))
      Sqrt2m2 = 2.0D0*dsqrt(2.0D0)
      AbsZmZout = dabs(ZmZout)
      AbsZpZout = dabs(ZpZout)
      Sqrt2 = dsqrt(2.0D0)
      SqrtSqrZmZoutdH = dsqrt(SqrZmZoutdH)
      SqrtSqrZpZoutdH = dsqrt(SqrZpZoutdH)

      cos_SqrZmZoutd4H = dcos(SqrZmZoutd4H)
      cos_SqrZpZoutd4H = dcos(SqrZpZoutd4H)

      sin_SqrZmZoutd4H = dsin(SqrZmZoutd4H)
      sin_SqrZpZoutd4H = dsin(SqrZpZoutd4H)

      call FCS(C_AbsZmZoutdSqrtH2PI, S_AbsZmZoutdSqrtH2PI, AbsZmZoutdSqrtH2PI)
      call FCS(C_AbsZpZoutdSqrtH2PI, S_AbsZpZoutdSqrtH2PI, AbsZpZoutdSqrtH2PI)

      I1 = (1, -1)*Sqrt2PI*((ZmZout*(-1 + (1, 1)*C_AbsZmZoutdSqrtH2PI + (1, -1)*S_AbsZmZoutdSqrtH2PI))/AbsZmZout + &
                            (ZpZout*(-1 + (1, 1)*C_AbsZpZoutdSqrtH2PI + (1, -1)*S_AbsZpZoutdSqrtH2PI))/AbsZpZout)

      I2 = (AbsZmZout*ZpZout*(-(SqrtPI*ZmZout**2) + 2*SqrtPI*ZmZout**2*S_AbsZmZoutdSqrtH2PI - (0,2)*DT*Sqrt2*SqrtSqrZmZoutdH*Sin_SqrZmZoutd4H) + &
         ZmZout*(SqrtH*ZpZout*((0,-1)*SqrtPI*SqrtSqrZmZoutdH*z - (0,1)*SqrtPI*SqrtSqrZpZoutdH*z + (0,1)*SqrtPI*SqrtSqrZmZoutdH*Zout - (0,1)*SqrtPI*SqrtSqrZpZoutdH*Zout + &
         Sqrt2m2*z*cos_SqrZmZoutd4H - 2*Sqrt2*Zout*Cos_SqrZmZoutd4H + Sqrt2m2*z*Cos_SqrZpZoutd4H + Sqrt2m2*Zout*cos_SqrZpZoutd4H + &
          (0, 2)*SqrtPI*SqrtSqrZmZoutdH*ZmZout*C_AbsZmZoutdSqrtH2PI + (0, 2)*SqrtPI*SqrtSqrZpZoutdH*ZpZout*C_AbsZpZoutdSqrtH2PI) + &
            AbsZpZout*(-(SqrtPI*ZpZout**2) + 2*SqrtPI*ZpZout**2*S_AbsZpZoutdSqrtH2PI - (0,2)*DT*Sqrt2*SqrtSqrZpZoutdH*sin_SqrZpZoutd4H)))/(dsqrt(2.0D0)*(-z + Zout)*ZpZout)

      I3 = (1, 1)*dsqrt(Pi/2.)*(-(Sqrt2PI*(ZmZout*C_AbsZmZoutdSqrtH2PI + ZpZout*C_AbsZpZoutdSqrtH2PI - &
                                           (0, 1)*((-1, -1)*z + ZmZout*S_AbsZmZoutdSqrtH2PI + ZpZout*S_AbsZpZoutdSqrtH2PI))) + &
         (2*SqrtH*ZmZout*((0,1)*cos_SqrZmZoutd4H + sin_SqrZmZoutd4H))/AbsZmZout + (2*SqrtH*ZpZout*((0,1)*cos_SqrZpZoutd4H + sin_SqrZpZoutd4H))/AbsZpZout)

      I4 = (0.25,0.25)*AbsZmZout*Sqrt2PI*SqrtH*ZmZout*cos_SqrZmZoutd4H + (0.25,0.25)*((0,1)*AbsZmZout**2 + 2*DT)*Pi*ZmZout*C_AbsZmZoutdSqrtH2PI + &
       (0.25,0.25)*((0,1)*AbsZpZout**2 + 2*DT)*Pi*ZpZout*C_AbsZpZoutdSqrtH2PI + (0.25,0.25)*Pi*(AbsZmZout**2*ZmZout - (0,2)*DT*ZmZout)*S_AbsZmZoutdSqrtH2PI + &
       (0.25,0.25)*Pi*(AbsZpZout**2*ZpZout - (0,2)*DT*ZpZout)*S_AbsZpZoutdSqrtH2PI + ((0.25,-0.25)*AbsZmZout**3*Sqrt2PI*SqrtH*sin_SqrZmZoutd4H)/ZmZout + &
       ((0,-2)*Pi*z*((0,-2)*DT + z**2 + 3*Zout**2)*ZpZout + (1,1)*AbsZpZout*Sqrt2PI*SqrtH*ZpZout**2*cos_SqrZpZoutd4H + (1,-1)*AbsZpZout**3*Sqrt2PI*SqrtH*sin_SqrZpZoutd4H)/(4.*ZpZout)

      C1 = SqrtPII*f1
      C2 = -SqrtPII/DT*f0 + SqrtPII/DT*f1
      C3 = -2.0D0/3.0D0*dfdz0 - 4.0D0/3.0D0*dfdz1
      C4 = -4.0D0/(3.0D0*DT)*dfdz0 + 4.0D0/(3.0D0*DT)*dfdz1
      Int = Zout/(2.0D0*pi)*(C1*I1 + C2*I2 + C3*I3 + C4*I4)
      print *, Int
   end subroutine integral_A4
end module integral_A4_mod
