program gyr
   use I12_space_mod
   use integral_A4_Zend_mod
   use Green
   use cubic_spline
   use integral_A4_mod
   use ifport
   implicit none

   integer(4) hours, minutes, seconds, k
   real(8) start_time, stop_time, calc_time

   integer(4) N, i, j, err_alloc, NS, iz
   real(8) :: pi = dacos(-1.0D0), Zout, DZ, DT, Z, DZS!, f0, f1, dfdz0, dfdz1
   complex(8) :: Is!, Int, C1, C2, C3, C4, I1, I2, I3, I4

   complex(8), allocatable, target :: I1s(:, :), I2s(:, :), fs(:), fa(:), B(:), C(:), D(:), f0(:), f1(:)
   real(8), allocatable, target :: ZAxis(:), reb(:), rec(:), red(:), imb(:), imc(:), imd(:), x(:)

   N = 1000
   Zout = 47
   DZ = Zout/N
   DT = 0.01
   iz = 500
   Z = iz*DZ
   NS = 50000
   DZS = Zout/NS

   allocate (I1s(N, N + 1), I2s(N, N + 1), ZAxis(N + 1), fs(N + 1), B(N + 1), C(N + 1), D(N + 1), f0(N), f1(N), &
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
   call I12(I1s, I2s, N, DZ, DT)
   stop_time = dclock()
   calc_time = stop_time - start_time
   hours = calc_time/3600
   minutes = (calc_time - hours*3600)/60
   seconds = calc_time - hours*3600 - minutes*60
   print *, 'Execution time:', hours, 'h :', minutes, 'm :', seconds, 's'
   write (*, '(/)')

   Is = 0
   do j = 1, N
      f0(j) = j*fs(j) - (j - 1)*fs(j + 1)
      f1(j) = (fs(j + 1) - fs(j))/DZ
      Is = Is + f0(j)*I1s(j, iz + 1) + f1(j)*I2s(j, iz + 1)
   end do

   write (*, '(a,2f18.10)') 'Is = ', Is
   print *, DZ

   !z = 47;
   !Zout = 47;
   !DT = 0.01;
   !f0 = 0.1;
   !f1 = 0.15;
   !dfdz0 = 0.1;
   !dfdz1 = 0.2;
   !call integral_A4(Int, z, f0, f1, dfdz0, dfdz1, Zout, DT)

   !Zout = 47.0D0;
   !DT = 0.1D0;
   !f0 = 0.1D0;
   !f1 = 0.15D0;
   !dfdz0 = 0.1D0;
   !dfdz1 = 0.2D0;
   !C1 = cdsqrt(pi*(0, 1))*f1;
   !C2 = -cdsqrt(pi*(0, 1))/DT*f0 + cdsqrt(pi*(0, 1))/DT*f1;
   !C3 = -2.0D0/3.0D0*dfdz0 - 4.0D0/3.0D0*dfdz1;
   !C4 = -4.0D0/(3.0D0*DT)*dfdz0 + 4.0D0/(3.0D0*DT)*dfdz1;
   !call integral_A4_Zend(I1, I2, I3, I4, f0, f1, dfdz0, dfdz1, Zout, DT)
   !
   !Int = Zout/(2.0D0*pi)*(C1*I1 + C2*I2 + C3*I3 + C4*I4)
   !
   !print *, Int

   open (1, file='test.dat')
   do i = 1, N
      write (1, '(1p2e17.8)') I1s(i, iz + 1)
   end do
   close (1)

   !open (1, file='test.dat')
   !do i = 1, N + 1
   !   write (1, '(f12.6,1p2e17.8)') ZAxis(i), fs(i)*G(ZAxis(i), Z, DT)
   !end do
   !close (1)
   !
   !open (1, file='testc.dat')
   !do i = 1, NS + 1
   !   write (1, '(f12.6,1p2e17.8)') (i - 1)*DZS, fa(i)
   !end do
   !close (1)

end program gyr
