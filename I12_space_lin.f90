module I12_space_mod
   use Fresnel
contains
   subroutine I12(I1, I2, N, DZ, DT)
      implicit none

      integer(4), intent(in) :: N
      real(8), intent(in) :: DZ, DT
      complex(8), dimension(:, :), intent(out) :: I1, I2

      integer(4) i, j

      do i = 1, N + 1
         do j = 1, N
            call calc_I12(I1(j, i), I2(j, i), j, i - 1, DZ, DT)
         end do
      end do
   end subroutine I12

   subroutine calc_I12(I1, I2, j, i, DZ, DT)
      implicit none
! Arguments declarations
      real(8), intent(in) :: DT
      real(8), intent(in) :: DZ
      integer(4), intent(in) :: j, i
      complex(8), intent(out) :: I1
      complex(8), intent(out) :: I2

! Variable declarations
      real(8) :: C_ZjmDZmZdSqrtDT2PI
      real(8) :: C_ZjmDZpZdSqrtDT2PI
      real(8) :: C_ZjmZdSqrtDT2PI
      real(8) :: C_ZjpZdSqrtDT2PI
      real(8) :: cos_SqrZjmDZmZd4DT
      real(8) :: cos_SqrZjmDZpZd4DT
      real(8) :: cos_SqrZjmZd4DT
      real(8) :: cos_SqrZjpZd4DT
      real(8) :: S_ZjmDZmZdSqrtDT2PI
      real(8) :: S_ZjmDZpZdSqrtDT2PI
      real(8) :: S_ZjmZdSqrtDT2PI
      real(8) :: S_ZjpZdSqrtDT2PI
      real(8) :: sin_SqrZjmDZmZd4DT
      real(8) :: sin_SqrZjmDZpZd4DT
      real(8) :: sin_SqrZjmZd4DT
      real(8) :: sin_SqrZjpZd4DT
      real(8) :: SqrZjmDZmZd4DT
      real(8) :: SqrZjmDZpZd4DT
      real(8) :: SqrZjmZd4DT
      real(8) :: SqrZjmZd8DT
      real(8) :: SqrZjpZd4DT
      real(8) :: ZjmDZmZdSqrtDT2PI
      real(8) :: ZjmDZpZdSqrtDT2PI
      real(8) :: ZjmZdSqrtDT2PI
      real(8) :: ZjpZdSqrtDT2PI
      real(8) :: SqrtDTSqrt2Pi
      real(8) :: pi = dacos(-1.0D0)

      ZjmZdSqrtDT2PI = (j - i)*DZ/dsqrt(DT*2*pi)
      ZjmDZmZdSqrtDT2PI = (j - 1 - i)*DZ/dsqrt(DT*2*pi)
      ZjpZdSqrtDT2PI = (j + i)*DZ/dsqrt(DT*2*pi)
      ZjmDZpZdSqrtDT2PI = (j - 1 + i)*DZ/dsqrt(DT*2*pi)
      SqrZjpZd4DT = (i + j)*(i + j)*DZ*DZ/(4*DT)
      SqrZjmDZpZd4DT = (j - 1 + i)*(j - 1 + i)*DZ*DZ/(4*DT)
      SqrZjmZd4DT = (j - i)*(j - i)*DZ*DZ/(4*DT)
      SqrZjmDZmZd4DT = (j - 1 - i)*(j - 1 - i)*DZ*DZ/(4*DT)

      call FCS(C_ZjpZdSqrtDT2PI, S_ZjpZdSqrtDT2PI, ZjpZdSqrtDT2PI)
      call FCS(C_ZjmZdSqrtDT2PI, S_ZjmZdSqrtDT2PI, ZjmZdSqrtDT2PI)
      call FCS(C_ZjmDZpZdSqrtDT2PI, S_ZjmDZpZdSqrtDT2PI, ZjmDZpZdSqrtDT2PI)
      call FCS(C_ZjmDZmZdSqrtDT2PI, S_ZjmDZmZdSqrtDT2PI, ZjmDZmZdSqrtDT2PI)

      cos_SqrZjmDZpZd4DT = dcos(SqrZjmDZpZd4DT)
      cos_SqrZjmDZmZd4DT = dcos(SqrZjmDZmZd4DT)
      sin_SqrZjmDZpZd4DT = dsin(SqrZjmDZpZd4DT)
      sin_SqrZjmDZmZd4DT = dsin(SqrZjmDZmZd4DT)
      sin_SqrZjmZd4DT = dsin(SqrZjmZd4DT)
      sin_SqrZjpZd4DT = dsin(SqrZjpZd4DT)
      cos_SqrZjmZd4DT = dcos(SqrZjmZd4DT)
      cos_SqrZjpZd4DT = dcos(SqrZjpZd4DT)

      SqrtDTSqrt2Pi = dsqrt(DT)*dsqrt(2*Pi)*DZ

      I1 = cdsqrt((0, 1)/2.0D0)*(-C_ZjmDZmZdSqrtDT2PI + C_ZjmDZpZdSqrtDT2PI + C_ZjmZdSqrtDT2PI - C_ZjpZdSqrtDT2PI - &
                                 (0, 1)*(-S_ZjmDZmZdSqrtDT2PI + S_ZjmDZpZdSqrtDT2PI + S_ZjmZdSqrtDT2PI - S_ZjpZdSqrtDT2PI))

      I2 = (cdsqrt((0, 1)/DT)*(-(SqrtDTSqrt2Pi*i*C_ZjmDZmZdSqrtDT2PI) - SqrtDTSqrt2Pi*i*C_ZjmDZpZdSqrtDT2PI + SqrtDTSqrt2Pi*i*C_ZjmZdSqrtDT2PI + &
                               SqrtDTSqrt2Pi*i*C_ZjpZdSqrtDT2PI + (0, 1)*(2*DT*(-Cos_SqrZjmDZmZd4DT + Cos_SqrZjmZd4DT) + &
                                                                       SqrtDTSqrt2Pi*i*(S_ZjmDZmZdSqrtDT2PI - S_ZjmZdSqrtDT2PI)) + &
                 (0, 1)*(2*DT*(Cos_SqrZjmDZpZd4DT - Cos_SqrZjpZd4DT) + SqrtDTSqrt2Pi*i*(S_ZjmDZpZdSqrtDT2PI - S_ZjpZdSqrtDT2PI)) + &
                      2*DT*(-Sin_SqrZjmDZmZd4DT + Sin_SqrZjmZd4DT) + 2*DT*(Sin_SqrZjmDZpZd4DT - Sin_SqrZjpZd4DT)))/(2.0D0*dsqrt(Pi))

   end subroutine calc_I12
end module I12_space_mod
