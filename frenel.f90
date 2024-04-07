        SUBROUTINE frenel(x, s, c)
           USE nrtype; USE nrutil, ONLY: nrerror
           IMPLICIT NONE
           REAL(DP), INTENT(IN) :: x
           REAL(DP), INTENT(OUT) :: s, c
           INTEGER(I4B), PARAMETER :: MAXIT = 100
           REAL(DP), PARAMETER :: EPS = epsilon(x), FPMIN = tiny(x), BIG = huge(x)*EPS, &
                                  XMIN = 1.5
           INTEGER(I4B) :: k, n
           REAL(DP) :: a, ax, fact, pix2, sign, sum, sumc, sums, term, test
           COMPLEX(DPC) :: b, cc, d, h, del, cs
           LOGICAL(LGT) :: odd
           ax = dabs(x)
           if (ax < dsqrt(FPMIN)) then
              s = 0.0
              c = ax
           else if (ax <= XMIN) then
              sum = 0.0
              sums = 0.0
              sumc = ax
              sign = 1.0
              fact = PIO2*ax*ax
              odd = .true.
              term = ax
              n = 3
              do k = 1, MAXIT
                 term = term*fact/k
                 sum = sum + sign*term/n
                 test = abs(sum)*EPS
                 if (odd) then
                    sign = -sign
                    sums = sum
                    sum = sumc
                 else
                    sumc = sum
                    sum = sums
                 end if
                 if (term < test) exit
                 odd = .not. odd
                 n = n + 2
              end do
              if (k > MAXIT) call nrerror('frenel: series failed')
              s = sums
              c = sumc
           else
              pix2 = PI*ax*ax
              b = cmplx(1.0_dp, -pix2, kind=DPC)
              cc = BIG
              d = 1.0_dp/b
              h = d
              n = -1
              do k = 2, MAXIT
                 n = n + 2
                 a = -n*(n + 1)
                 b = b + 4.0_dp
                 d = 1.0_dp/(a*d + b)
                 cc = b + a/cc
                 del = cc*d
                 h = h*del
                 if (absc(del - 1.0_dp) <= EPS) exit
              end do
              if (k > MAXIT) call nrerror('cf failed in frenel')
              h = h*cmplx(ax, -ax, kind=DPC)
              cs = cmplx(0.5_dp, 0.5_dp, kind=DPC)*(1.0_dp - &
                                                    cmplx(dcos(0.5_dp*pix2), dsin(0.5_dp*pix2), kind=DPC)*h)
              c = dreal(cs)
              s = dimag(cs)
           end if
           if (x < 0.0) then
              c = -c
              s = -s
           end if
        CONTAINS
!BL
           FUNCTION absc(z)
              IMPLICIT NONE
              COMPLEX(DPC), INTENT(IN) :: z
              REAL(DP) :: absc
              absc = dabs(dreal(z)) + dabs(dimag(z))
           END FUNCTION absc
        END SUBROUTINE frenel
