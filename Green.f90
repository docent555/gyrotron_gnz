module Green
contains
   complex function G(x, z, t)
      implicit none

      real(8), intent(in) :: x, z, t

      real(8) :: pi = dacos(-1.0D0)
      complex(8) :: Im1 = (0, 1), argp, argm

      argp = -Im1*(x + z)*(x + z)/(4*t)
      argm = -Im1*(x - z)*(x - z)/(4*t)
      G = 1./2.*cdsqrt(Im1/(pi*t))*(cdexp(argm) - cdexp(argp))
   end function G
end module Green
