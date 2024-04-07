module I1234_space_mod
   use Fresnel
contains
   subroutine I1234(I1, I2, I3, I4, N, DZ, DT)
      implicit none

      integer(4), intent(in) :: N
      real(8), intent(in) :: DZ, DT
      complex(8), dimension(:, :), intent(out) :: I1, I2, I3, I4

      integer(4) i, j

      do i = 1, N + 1
         do j = 1, N
            call calc_I1234(I1(j, i), I2(j, i), I3(j, i), I4(j, i), j, i - 1, DZ, DT)
         end do
      end do
   end subroutine I1234

   real(8) function FresnelC(t)
      implicit none

      real(8), intent(in) :: t
      real(8) :: C, S

      real(8) :: pi = dacos(-1.0D0), x

      x = pi*t*t/2.0

      call CS(C, S, x)

      if (t > 0) then
         FresnelC = C
      else
         FresnelC = -C
      end if
   end function FresnelC

   real(8) function FresnelS(t)
      implicit none

      real(8), intent(in) :: t
      real(8) :: C, S

      real(8) :: pi = dacos(-1.0D0), x

      x = pi*t*t/2.0

      call CS(C, S, x)

      if (t > 0) then
         FresnelS = S
      else
         FresnelS = -S
      end if
   end function FresnelS

   subroutine calc_I1234(I1, I2, I3, I4, j, i, DZ, DT)
      implicit none
! Arguments declarations
      real(8), intent(in) :: DT  !<
      real(8), intent(in) :: DZ  !<
      integer(4), intent(in) :: j  !<
      integer(4), intent(in) :: i  !<
      complex(8), intent(out) :: I1  !< !m2f: check dim(:)!m
      complex(8), intent(out) :: I2  !< !m2f: check dim(:)!m
      complex(8), intent(out) :: I3  !< !m2f: check dim(:)!m
      complex(8), intent(out) :: I4  !< !m2f: check dim(:)!m

! Variable declarations
      !real(8) :: C_Zjm1mZdSqrtDT2PI  !<
      !real(8) :: C_Zjm1pZdSqrtDT2PI  !<
      !real(8) :: C_ZjmZdSqrtDT2PI  !<
      !real(8) :: C_ZjpZdSqrtDT2PI  !<
      !real(8) :: cos_SqrZjmDZmZd4DT  !<
      !real(8) :: cos_SqrZjmDZmZd8DT  !<
      !real(8) :: cos_SqrZjmDZpZd4DT  !<
      !real(8) :: cos_SqrZjmZd4DT  !<
      !real(8) :: cos_SqrZjmZd8DT  !<
      !real(8) :: cos_SqrZjpZd4DT  !<
      !real(8) :: S_Zjm1mZdSqrtDT2PI  !<
      !real(8) :: S_Zjm1pZdSqrtDT2PI  !<
      !real(8) :: S_ZjmZdSqrtDT2PI  !<
      !real(8) :: S_ZjpZdSqrtDT2PI  !<
      !real(8) :: sin_SqrZjmDZmZd4DT  !<
      !real(8) :: sin_SqrZjmDZmZd8DT  !<
      !real(8) :: sin_SqrZjmDZpZd4DT  !<
      !real(8) :: sin_SqrZjmZd4DT  !<
      !real(8) :: sin_SqrZjmZd8DT  !<
      !real(8) :: sin_SqrZjpZd4DT  !<
      !real(8) :: SqrZjmDZmZd4DT  !<
      !real(8) :: SqrZjmDZmZd8DT  !<
      !real(8) :: SqrZjmDZpZd4DT  !<
      !real(8) :: SqrZjmZd4DT  !<
      !real(8) :: SqrZjmZd8DT  !<
      !real(8) :: SqrZjpZd4DT  !<
      !real(8) :: Zjm1mZdSqrtDT2PI  !<
      !real(8) :: Zjm1pZdSqrtDT2PI  !<
      !real(8) :: ZjmZdSqrtDT2PI  !<
      !real(8) :: ZjpZdSqrtDT2PI  !<
      !real(8) :: SqrtDT2PI
      !real(8) :: Sqrt2
      real(8) :: pi = dacos(-1.0D0)
      !complex(8) :: SqrtI
      !complex(8) :: Im1 = (0, 1)
      !
      !ZjmZdSqrtDT2PI = (zj - Z)/dsqrt(DT*2*pi)
      !Zjm1mZdSqrtDT2PI = (zj - DZ - Z)/dsqrt(DT*2*pi)
      !ZjpZdSqrtDT2PI = (zj + Z)/dsqrt(DT*2*pi)
      !Zjm1pZdSqrtDT2PI = (zj - DZ + Z)/dsqrt(DT*2*pi)
      !SqrZjmZd8DT = (zj - Z)*(zj - Z)/(8*DT)
      !SqrZjmDZmZd8DT = (zj - DZ - Z)*(zj - DZ - Z)/(8*DT)
      !SqrZjpZd4DT = (Z + zj)*(Z + zj)/(4*DT)
      !SqrZjmDZpZd4DT = (zj - DZ + Z)*(zj - DZ + Z)/(4*DT)
      !SqrZjmZd4DT = (zj - Z)*(zj - Z)/(4*DT)
      !SqrZjmDZmZd4DT = (zj - DZ - Z)*(zj - DZ - Z)/(4*DT)

      !C_Zjm1mZdSqrtDT2PI = fresnelc(Zjm1mZdSqrtDT2PI)
      !S_Zjm1mZdSqrtDT2PI = fresnels(Zjm1mZdSqrtDT2PI)
      !call FCS(C_Zjm1mZdSqrtDT2PI,S_Zjm1mZdSqrtDT2PI,Zjm1mZdSqrtDT2PI)
      !
      !!C_Zjm1pZdSqrtDT2PI = fresnelc(Zjm1pZdSqrtDT2PI)
      !!S_Zjm1pZdSqrtDT2PI = fresnels(Zjm1pZdSqrtDT2PI)
      !call FCS(C_Zjm1pZdSqrtDT2PI,S_Zjm1pZdSqrtDT2PI,Zjm1pZdSqrtDT2PI)
      !
      !!C_ZjmZdSqrtDT2PI = fresnelc(ZjmZdSqrtDT2PI)
      !!S_ZjmZdSqrtDT2PI = fresnels(ZjmZdSqrtDT2PI)
      !call FCS(C_ZjmZdSqrtDT2PI,S_ZjmZdSqrtDT2PI,ZjmZdSqrtDT2PI)
      !
      !!C_ZjpZdSqrtDT2PI = fresnelc(ZjpZdSqrtDT2PI)
      !!S_ZjpZdSqrtDT2PI = fresnels(ZjpZdSqrtDT2PI)
      !call FCS(C_ZjpZdSqrtDT2PI,S_ZjpZdSqrtDT2PI,ZjpZdSqrtDT2PI)
      !
      !cos_SqrZjmDZmZd8DT = dcos(SqrZjmDZmZd8DT)
      !cos_SqrZjmDZpZd4DT = dcos(SqrZjmDZpZd4DT)
      !cos_SqrZjmZd8DT = dcos(SqrZjmZd8DT)
      !cos_SqrZjpZd4DT = dcos(SqrZjpZd4DT)
      !sin_SqrZjmDZmZd4DT = dsin(SqrZjmDZmZd4DT)
      !sin_SqrZjmDZmZd8DT = dsin(SqrZjmDZmZd8DT)
      !sin_SqrZjmDZpZd4DT = dsin(SqrZjmDZpZd4DT)
      !sin_SqrZjmZd4DT = dsin(SqrZjmZd4DT)
      !sin_SqrZjmZd8DT = dsin(SqrZjmZd8DT)
      !sin_SqrZjpZd4DT = dsin(SqrZjpZd4DT)
      !cos_SqrZjmDZmZd4DT = dcos(SqrZjmDZmZd4DT)
      !cos_SqrZjmZd4DT = dcos(SqrZjmZd4DT)
      !
      !SqrtDT2PI = dsqrt(DT*2*pi)
      !SqrtI = cdsqrt((0,1.0D0))
      !Sqrt2 = dsqrt(2.0D0)

      I1 = (Sqrt((0,1)/DT)*Sqrt(DT)*(-FresnelC((DZ*(-1 - i + j))/(Sqrt(DT)*Sqrt(2*Pi))) + FresnelC((DZ*(-i + j))/(Sqrt(DT)*Sqrt(2*Pi))) + &
                                FresnelC((DZ*(-1 + i + j))/(Sqrt(DT)*Sqrt(2*Pi))) - FresnelC((DZ*(i + j))/(Sqrt(DT)*Sqrt(2*Pi))) + &
           (0,1)*(FresnelS((DZ*(-1 - i + j))/(Sqrt(DT)*Sqrt(2*Pi))) - FresnelS((DZ*(-i + j))/(Sqrt(DT)*Sqrt(2*Pi))) - FresnelS((DZ*(-1 + i + j))/(Sqrt(DT)*Sqrt(2*Pi))) + &
                                            FresnelS((DZ*(i + j))/(Sqrt(DT)*Sqrt(2*Pi))))))/Sqrt(2.0D0)

      I2 = (Sqrt((0,1)/DT)*Sqrt(DT)*((0,2)*Sqrt(DT)*(cdexp(-((0,1)*DZ**2*(i - j)**2)/(4.*DT)) - cdexp(-((0,1)*DZ**2*(1 + i - j)**2)/(4.*DT)) + cdexp(-((0,1)*DZ**2*(-1 + i + j)**2)/(4.*DT)) - &
       cdexp(-((0, 1)*DZ**2*(i + j)**2)/(4.*DT))) + DZ*(-1 - i + j)*Sqrt(2*Pi)*FresnelC((DZ*(-1 - i + j))/(Sqrt(DT)*Sqrt(2*Pi))) + &
                                     DZ*(1 + i - j)*Sqrt(2*Pi)*FresnelC((DZ*(-i + j))/(Sqrt(DT)*Sqrt(2*Pi))) + &
           DZ*Sqrt(2*Pi)*(-((-1 + i + j)*FresnelC((DZ*(-1 + i + j))/(Sqrt(DT)*Sqrt(2*Pi)))) + (-1 + i + j)*FresnelC((DZ*(i + j))/(Sqrt(DT)*Sqrt(2*Pi))) + &
              (0,1)*((1 + i - j)*FresnelS((DZ*(-1 - i + j))/(Sqrt(DT)*Sqrt(2*Pi))) + (-1 - i + j)*FresnelS((DZ*(-i + j))/(Sqrt(DT)*Sqrt(2*Pi))) + &
   (-1 + i + j)*(FresnelS((DZ*(-1 + i + j))/(Sqrt(DT)*Sqrt(2*Pi))) - FresnelS((DZ*(i + j))/(Sqrt(DT)*Sqrt(2*Pi))))))))/(2.*Sqrt(Pi))

      I3 = (Sqrt((0,1)/DT)*Sqrt(DT)*(-(((0,-2)*DT + DZ**2*(1 + i - j)**2)*Sqrt(2*Pi)*FresnelC((DZ*(-1 - i + j))/(Sqrt(DT)*Sqrt(2*Pi)))) + &
                                    ((0, -2)*DT + DZ**2*(1 + i - j)**2)*Sqrt(2*Pi)*FresnelC((DZ*(-i + j))/(Sqrt(DT)*Sqrt(2*Pi))) - &
                                ((0, 2)*DT - DZ**2*(-1 + i + j)**2)*Sqrt(2*Pi)*FresnelC((DZ*(-1 + i + j))/(Sqrt(DT)*Sqrt(2*Pi))) + &
                                     ((0, 2)*DT - DZ**2*(-1 + i + j)**2)*Sqrt(2*Pi)*FresnelC((DZ*(i + j))/(Sqrt(DT)*Sqrt(2*Pi))) + &
                               (2*DT + (0, 1)*DZ**2*(1 + i - j)**2)*Sqrt(2*Pi)*FresnelS((DZ*(-1 - i + j))/(Sqrt(DT)*Sqrt(2*Pi))) - &
                                   (2*DT + (0, 1)*DZ**2*(1 + i - j)**2)*Sqrt(2*Pi)*FresnelS((DZ*(-i + j))/(Sqrt(DT)*Sqrt(2*Pi))) - &
                              (2*DT + (0, 1)*DZ**2*(-1 + i + j)**2)*Sqrt(2*Pi)*FresnelS((DZ*(-1 + i + j))/(Sqrt(DT)*Sqrt(2*Pi))) + &
                                   (2*DT + (0, 1)*DZ**2*(-1 + i + j)**2)*Sqrt(2*Pi)*FresnelS((DZ*(i + j))/(Sqrt(DT)*Sqrt(2*Pi))) - &
           (4*Sqrt(DT)*DZ*((0,1)*(0,1)*Cos((DZ**2*(0,1)*(-1 + j))/(2.*DT)) + (-1 + j)*Sin((DZ**2*(0,1)*(-1 + j))/(2.*DT))))/cdexp(((0,1)*DZ**2*((0,1)**2 + (-1 + j)**2))/(4.*DT)) + &
           (4*Sqrt(DT)*DZ*((0,1)*(0,1)*Cos((DZ**2*(0,1)*j)/(2.*DT)) + (-2 + j)*Sin((DZ**2*(0,1)*j)/(2.*DT))))/cdexp(((0,1)*DZ**2*((0,1)**2 + j**2))/(4.*DT))))/(2.*Sqrt(Pi))

 I4 = (Sqrt((0, 1)/DT)*Sqrt(DT)*((-2*Sqrt(DT)*(4*DT + (0, 1)*DZ**2*(1 + i - j)**2))/cdexp(((0, 1)*DZ**2*(1 + i - j)**2)/(4.*DT)) + &
                                (2*Sqrt(DT)*(4*DT + (0, 1)*DZ**2*(-1 + i + j)**2))/cdexp(((0, 1)*DZ**2*(-1 + i + j)**2)/(4.*DT)) + &
     (2*Sqrt(DT)*(4*DT + (0, 1)*DZ**2*(3 + (0, 1)**2 + (0, 1)*(3 - 2*j) + (-3 + j)*j)))/cdexp(((0, 1)*DZ**2*(i - j)**2)/(4.*DT)) - &
    (2*Sqrt(DT)*(4*DT + (0, 1)*DZ**2*(3 + (0, 1)**2 + (-3 + j)*j + (0, 1)*(-3 + 2*j))))/cdexp(((0, 1)*DZ**2*(i + j)**2)/(4.*DT)) - &
               DZ*(-1 + i + j)*((0, -6)*DT + DZ**2*(-1 + i + j)**2)*Sqrt(2*Pi)*FresnelC((DZ*(-1 + i + j))/(Sqrt(DT)*Sqrt(2*Pi))) + &
                    DZ*(-1 + i + j)*((0, -6)*DT + DZ**2*(-1 + i + j)**2)*Sqrt(2*Pi)*FresnelC((DZ*(i + j))/(Sqrt(DT)*Sqrt(2*Pi))) - &
                                      DZ*((0, -6)*DT + DZ**2*(1 + i - j)**2)*(1 + i - j)*Sqrt(2*Pi)* &
                  (FresnelC((DZ*(-1 - i + j))/(Sqrt(DT)*Sqrt(2*Pi))) - (0, 1)*FresnelS((DZ*(-1 - i + j))/(Sqrt(DT)*Sqrt(2*Pi)))) + &
           DZ*((0,-6)*DT + DZ**2*(1 + i - j)**2)*(1 + i - j)*Sqrt(2*Pi)*(FresnelC((DZ*(-i + j))/(Sqrt(DT)*Sqrt(2*Pi))) - (0,1)*FresnelS((DZ*(-i + j))/(Sqrt(DT)*Sqrt(2*Pi)))) + &
              DZ*(-1 + i + j)*(6*DT + (0, 1)*DZ**2*(-1 + i + j)**2)*Sqrt(2*Pi)*FresnelS((DZ*(-1 + i + j))/(Sqrt(DT)*Sqrt(2*Pi))) - &
       DZ*(-1 + i + j)*(6*DT + (0, 1)*DZ**2*(-1 + i + j)**2)*Sqrt(2*Pi)*FresnelS((DZ*(i + j))/(Sqrt(DT)*Sqrt(2*Pi)))))/(2.*Sqrt(Pi))
   end subroutine calc_I1234
end module I1234_space_mod
