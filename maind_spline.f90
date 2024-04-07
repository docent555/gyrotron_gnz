program gyr
   use I1234_space_mod
   use integral_A4_Zend_mod
   use Green
   use cubic_spline
   use integral_A4_mod
   use ifport
   implicit none

   integer(4) hours, minutes, seconds, k
   real(8) start_time, stop_time, calc_time

   integer(4) N, i, j, err_alloc, NS, iz
   real(8) :: pi = dacos(-1.0D0), Zout, DZ, DT, Z, DZS
   complex(8) :: Is, Int, C1, C2, C3, C4

   complex(8), allocatable, target :: I1(:, :), I2(:, :), I3(:, :), I4(:, :), fs(:), fa(:), B(:), C(:), D(:), f0(:), f1(:)
   real(8), allocatable, target :: ZAxis(:), reb(:), rec(:), red(:), imb(:), imc(:), imd(:), x(:)

   N = 1000
   Zout = 47
   DZ = Zout/N
   DT = 5
   iz = 80
   Z = iz*DZ
   NS = 50000
   DZS = Zout/NS

   allocate (I1(N, N + 1), I2(N, N + 1), I3(N, N + 1), I4(N, N + 1), ZAxis(N + 1), fs(N + 1), B(N + 1), C(N + 1), D(N + 1), f0(N), f1(N), &
             reb(N + 1), rec(N + 1), red(N + 1), imb(N + 1), imc(N + 1), imd(N + 1), x(NS), fa(NS + 1), stat=err_alloc)
   if (err_alloc /= 0) then
      print *, "allocation error"
      pause
      stop
   end if

   do i = 1, N + 1
      ZAxis(i) = (i - 1)*DZ
      fs(i) = (dcos(2*pi*ZAxis(i)/Zout)**2 + (0, 1)*dsin(2*pi*ZAxis(i)/Zout)**2)
   end do

   call spline(N + 1, ZAxis, dreal(fs), reb, rec, red)
   call spline(N + 1, ZAxis, dimag(fs), imb, imc, imd)

   do i = 1, NS + 1
      fa(i) = seval_cmplx(N + 1, (i - 1)*DZS, ZAxis, dreal(fs), dimag(fs), reb, rec, red, imb, imc, imd)*G((i - 1)*DZS, Z, DT)
   end do

   B = dcmplx(reb, imb)
   C = dcmplx(rec, imc)
   D = dcmplx(red, imd)

   start_time = dclock()
   call I1234(I1, I2, I3, I4, N, DZ, DT)
   stop_time = dclock()
   calc_time = stop_time - start_time
   hours = calc_time/3600
   minutes = (calc_time - hours*3600)/60
   seconds = calc_time - hours*3600 - minutes*60
   print *, 'Execution time:', hours, 'h :', minutes, 'm :', seconds, 's'
   write (*, '(/)')

   Is = 0
   do j = 1, N
      !f0(j) = j*fs(j) - (j - 1)*fs(j + 1)
      !f1(j) = (fs(j + 1) - fs(j))/DZ
      Is = Is + fs(j)*I1(j, iz + 1) + B(j)*I2(j, iz + 1) + C(j)*I3(j, iz + 1) + D(j)*I4(j, iz + 1)
   end do

   write (*, '(a,2f18.10)') 'Is = ', Is
   print *, DZ

   open (1, file='test.dat')
   do i = 1, N
      write (1, '(1p2e17.8)') I4(i, iz + 1)
   end do
   close (1)

   !open (1, file='test.dat')
   !do i = 1, N + 1
   !   write (1, '(f12.6,1p2e17.8)') ZAxis(i), fs(i)*G(ZAxis(i), Z, DT)
   !end do
   !close (1)

   !open (1, file='testc.dat')
   !do i = 1, NS + 1
   !   write (1, '(f12.6,1p2e17.8)') (i - 1)*DZS, fa(i)
   !end do
   !close (1)

end program gyr
