!=====================================================================
!  pi_demo.f90
!
!  Demonstrates the equidistribution of the fractional parts of
!  n*π (n = 1,2,…,N) by building a histogram.
!
!  Compile with:
!      gfortran -O2 -Wall -fmax-errors=1 pi_equidistribution.f90 -o pi_eq
!  Run, e.g.:
!      ./pi_eq 1000000 20
!
!  The first command‑line argument is the number of points N,
!  the second argument is the number of histogram bins M (default 10).
!=====================================================================

program pi_equidistribution
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)   ! double precision
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp

    integer :: N          ! number of points to generate
    integer :: M          ! number of histogram bins
    integer :: i, ibin
    real(dp) :: xn        ! current fractional part
    real(dp) :: dx        ! bin width = 1/M
    real(dp), allocatable :: hist(:)   ! histogram counts (as real for easy printing)
    character(len=*), parameter :: fmt = '(I12,1X,F12.8)'
    real(dp) :: ideal, dev, maxdev
    !-----------------------------------------------------------------
    !  Read command line arguments (or use defaults)
    !-----------------------------------------------------------------
    N = 1000000000   ! default: one million points
    M = 10        ! default: ten bins

    if (N <= 0 .or. M <= 0) then
        write(*,*) 'Error: N and M must be positive integers.'
        stop 1
    end if

    allocate(hist(M))
    hist = 0.0_dp
    dx = 1.0_dp / real(M, dp)

    !-----------------------------------------------------------------
    !  Main loop: compute {n*π} and increment the appropriate bin
    !-----------------------------------------------------------------
    do concurrent (i=1:N) !reduce (+:hist)
        xn = mod( real(i, dp) * pi, 1.0_dp )   ! fractional part of n*π
        ibin = int( xn / dx ) + 1              ! bin index (1 … M)

        ! Guard against the very unlikely case xn == 1.0 (should never happen)
        if (ibin > M) ibin = M

        hist(ibin) = hist(ibin) + 1.0_dp
    end do

    !-----------------------------------------------------------------
    !  Output the histogram
    !-----------------------------------------------------------------
    write(*,*) '---------------------------------------------------'
    write(*,*) 'Equidistribution test for { n * π } (n = 1 ..', N, ')'
    write(*,*) 'Number of bins :', M
    write(*,*) 'Bin width       :', dx
    write(*,*) '---------------------------------------------------'
    write(*,'(A)') '  Bin   Interval               Count    Frequency'
    write(*,'(A)') '---------------------------------------------------'

    do ibin = 1, M
        write(*,'(I4,1X,"[",F6.4,1X,"-",F6.4,")",1X,I12,1X,F12.8)') &
            ibin, (real(ibin-1,dp)*dx), (real(ibin,dp)*dx), &
            nint(hist(ibin)), hist(ibin)/real(N,dp)
    end do

    !-----------------------------------------------------------------
    !  Compute a simple quality measure: max deviation from the ideal 1/M
    !-----------------------------------------------------------------

    ideal = 1.0_dp / real(M, dp)
    maxdev = 0.0_dp
    do ibin = 1, M
        dev = abs( hist(ibin)/real(N,dp) - ideal )
        if (dev > maxdev) maxdev = dev
    end do

    write(*,*) '---------------------------------------------------'
    write(*,'(A,F12.8)') 'Ideal frequency per bin      : ', ideal
    write(*,'(A,F12.8)') 'Maximum absolute deviation   : ', maxdev
    write(*,*) '---------------------------------------------------'
    write(*,*) 'The deviations shrink roughly like 1/sqrt(N).'
    write(*,*) 'Increasing N (or decreasing M) will make the histogram'
    write(*,*) 'appear more uniform, confirming equidistribution.'

    deallocate(hist)
end program pi_equidistribution
  
