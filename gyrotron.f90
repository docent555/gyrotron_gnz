program gyrotron
   use functions
#if __INTEL_COMPILER
   use ifport
#endif

   implicit none

   integer(4) hours, minutes, seconds
   real(8) start_time, stop_time, calc_time

   integer(4) i, j, Ne, Nz, Nt, INTT, INTZ, err_alloc, OUTNt, OUTNz
   real(8) Lz, Tend, Delta, Ic, dz, dt, ZBEG, ZEND, tol
   real(8) :: pi = dacos(-1.0d0)

   real(8), allocatable :: ZAxis(:), TAxis(:), OUTZAxis(:), OUTTAxis(:)
   complex(8), allocatable :: InitialField(:), OUTF(:, :), OUTJ(:, :)

   namelist /param/ Ne, Lz, Tend, Delta, Ic, dz, dt, tol

   open (unit=1, file='input_fortran.in', status='old', err=101)
   read (unit=1, nml=param, err=102)
   close (unit=1)

   write (*, nml=param)

   write (*, '(/)')
#if __INTEL_COMPILER
   start_time = dclock()
#endif

   Nz = Lz/dz + 1
   Nt = Tend/dt + 1

   allocate (ZAxis(Nz), TAxis(Nt), InitialField(Nz), stat=err_alloc)
   if (err_alloc /= 0) then
      print *, "allocation error"
      pause
      stop
   end if

   do i = 1, Nz
      ZAxis(i) = (i - 1)*dz
   end do

   do i = 1, Nt
      TAxis(i) = (i - 1)*dt
   end do

   ZBEG = 0
   ZEND = 0.5
   where ((ZAxis .GT. ZBEG) .AND. (ZAxis .LT. ZEND)) InitialField = sin(pi*(ZAxis - ZBEG)/(ZEND - ZBEG))**2
   !InitialField = dcmplx(10.0D0, 0.0D0)

   !infield = [real(InitialField) imag(InitialField)];
   !save('init_field.in', 'infield', '-ascii')

   INTT = Nt/500
   INTZ = Nz/500
   !INTT = 1;
   !INTZ = 1;
   if (INTT < 1) then
      print *, 'Too small "Tend"'
      pause
      stop
   end if

   if (INTZ < 1) then
      print *, 'Too small "Lz"'
      pause
      stop
   end if

   if (INTT > 1 .and. INTZ > 1) then
      OUTNt = (Nt - 1)/INTT + 1
      OUTNz = (Nz - 1)/INTZ + 1
   elseif (INTT == 1 .and. INTZ > 1) then
      OUTNt = Nt
      OUTNz = (Nz - 1)/INTZ + 1
   elseif (INTT > 1 .and. INTZ == 1) then
      OUTNt = (Nt - 1)/INTT + 1
      OUTNz = Nz
   else
      OUTNt = Nt
      OUTNz = Nz
   end if

   allocate (OUTJ(OUTNz, OUTNt), OUTF(OUTNz, OUTNt), OUTZAxis(OUTNz), OUTTAxis(OUTNt), stat=err_alloc)
   if (err_alloc /= 0) then
      print *, "allocation error"
      pause
      stop
   end if

   do i = 1, OUTNz
      OUTZAxis(i) = (i - 1)*INTZ*dz; 
   end do

   do i = 1, OUTNt
      OUTTAxis(i) = (i - 1)*INTT*dt; 
   end do

   call gyroscr(Nz, Nt, Ne, ZAxis, TAxis, Delta, Ic, dt, dz, tol, INTT, INTZ, InitialField, OUTNz, OUTNt, OUTF, OUTJ)

#if __INTEL_COMPILER
   stop_time = dclock()
   calc_time = stop_time - start_time
   hours = calc_time/3600
   minutes = (calc_time - hours*3600)/60
   seconds = calc_time - hours*3600 - minutes*60

   write (*, '(/)')
   print *, 'Calcualting took:', hours, 'h :', minutes, 'm :', seconds, 's'
#endif

   open (1, file='fvsz.dat')
   do i = 1, OUTNz
      write (1, fmt='(1x,f14.6)', advance="no") OUTZAxis(i)
      do j = 1, OUTNt
         write (1, fmt='(1x,1p2e17.7,a)', advance="no") OUTF(i, j)
      end do
      write (1, *) ! Assumes default "advance='yes'".
   end do
   close (1)

   open (1, file='ivsz.dat')
   do i = 1, OUTNz
      write (1, fmt='(1x,f14.6)', advance="no") OUTZAxis(i)
      do j = 1, OUTNt
         write (1, fmt='(1x,1p2e17.7,a)', advance="no") OUTJ(i, j)
      end do
      write (1, *) ! Assumes default "advance='yes'".
   end do
   close (1)

   stop
101 print *, "error of file open"; pause; stop
102 print *, 'error of file read'; pause; stop
103 print *, 'error of file write'; pause; stop
end
