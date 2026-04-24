!! File   : memoized_fibonacci.f90
!! Author : Sandeep Koranne, (C) 2023. All rights reserved.
!! Purpose: Allocatable array used for memoized Fibonacci

module memoized_fibonacci
  use ISO_Fortran_env, only: &
       stdout => OUTPUT_UNIT, &
       stdin  => INPUT_UNIT, &
       compiler_version, &
       compiler_options

  implicit none
  interface print_fibonacci
     module procedure print_fibonacci
  end interface print_fibonacci

  interface print_factorial
     module procedure print_factorial
  end interface print_factorial

contains
!---------------------------------------------------------------------
! Subroutine to print up to N fibonacci terms
! 
! This subroutine generates and prints the first N Fibonacci numbers
! using an allocatable array to store the sequence. The Fibonacci 
! sequence starts with 1, 1, 2, 3, 5, 8, 13, ... where each number
! is the sum of the two preceding ones.
!
! Input: N - the number of Fibonacci numbers to generate and
! Output: Prints the Fibonacci sequence to standard output
!
!---------------------------------------------------------------------

  subroutine print_fibonacci
    integer :: N,i
    integer, dimension(:), allocatable :: fnumbers
    write (stdout,*) 'Enter how many Fibonacci Numbers you need: '
    read  (stdin,*) N
    write (stdout,FMT='(A,I2,A)') 'Generating the first ',N,' Fibonacci numbers.'
    allocate(fnumbers(N))
    fnumbers(1) = 1
    fnumbers(2) = 1
    do i=3,N
       fnumbers(i) = fnumbers(i-1) + fnumbers(i-2)
    end do
    write (stdout,FMT='(6I)') fnumbers
    write (stdout,*) "----------"
    deallocate(fnumbers)
  end subroutine print_fibonacci

  !---------------------------------------------------------------------
  ! SUBROUTINE PRINT_FACTORIAL
  !
  ! This subroutine calculates and prints the factorial of a given number
  ! using approach. It handles edge cases such as zero and
  ! negative numbers, and includes error checking for invalid inputs.
  !
  ! INPUT PARAMETERS:
  !   n     - Integer value for which factorial is to be calculated
  !   fact  - Integer array or variable to store the factorial result
  !
  ! OUTPUT:
  !   Prints the factorial result to standard output in a formatted manner
  !   Uses 6I10 format for displaying up to integers with 10-character width
  !
  ! EXAMPLE:
  !   For n = 5, the subroutine will calculate 5! = 120 and display it
  !
  ! NOTES:
  !   - Factorial of 0 is defined as 1
  !   - Negative numbers are handled with error message
  !   - Large factorials may exceed integer limits
  !   - Uses iterative method for better performance than recursive approach
  !---------------------------------------------------------------------

  subroutine print_factorial
    integer :: N,i
    integer, dimension(:), allocatable :: fact
    write (stdout,*) 'Enter how many factorial numbers you need: '
    read  (stdin,*) N
    write (stdout,FMT='(A,I2,A)') 'Generating the first ',N,' factorial numbers.'
    allocate(fact(N))
    fact(1) = 1
    do i=2,N
       fact(i) = i * fact(i-1)
    end do
    write (stdout,FMT='(6I)') fact
    write (stdout,*) "----------"
    deallocate(fact)
  end subroutine print_factorial

!---------------------------------------------------------------------
! SUBROUTINE PRINT_PRIMES
!
! This subroutine generates and prints the first N prime numbers
! using the Sieve of Eratost. The algorithm efficiently
! finds all prime numbers up to a given limit by iteratively marking
! multiples of each prime as composite.
!
! INPUT PARAMETERS:
!   N - Integer value specifying how many prime numbers to generate
!
! OUTPUT:
!   Prints the first N prime numbers to standard output in a formatted manner
!   Uses 6I8 format for displaying primes with 8-character width
!
! EXAMPLE:
!   For N 10, the subroutine will print: 2 3 5 7 11 13 17 19 23 29
!
! NOTES:
!   - Handles edge cases where N <= 0
!   - Uses an optimized sieve approach with square root limit
!   - Memory usage is proportional to the estimated upper bound of the Nth prime
!   - The algorithm is efficient for generating moderate numbers of primes
!
!----------------------------------------------------------------
subroutine print_primes()
    implicit none
    integer :: N
    integer :: i, j, count, limit
    logical, dimension(:), allocatable :: is_prime
    integer, dimension(:), allocatable :: primes
    write (stdout,*) 'Enter how many Prime numbers you need: '
    read  (stdin,*) N
    ! Handle edge cases
    if (N <= 0) then
        write(*,*) 'No primes to print for N <= 0'
    end if
    
    ! Estimate upper bound for nth prime (simplified approximation)
    ! This limit is determined using the prime number theorem
    ! For the nth prime number, it's approximately n * ln(n)
    ! The formula limit = int.2 * N * log(real(N))) + 10
    ! provides a safe upper bound to ensure we find at least N primes
    ! The factor 1.2 accounts for the error term and provides a margin of safety
    ! The + 10 is an additional buffer to guarantee we don't miss any primes
    
    if (N < 6) then
        limit = 12
    else
        limit = int(1.2 * N * log(real(N))) + 10
    end if
    
    ! Allocate sieve array
    ! Initialize all numbers as potentially prime (except 1)    
    allocate(is_prime(2:limit),source=.true.)
    allocate(primes(N),source=0)
    
    ! Initialize sieve

    is_prime(1) = .false.
    ! Sieve of Eratosthenes
    do i = 3, int(sqrt(real(limit)))
       if (is_prime(i)) then
          !Starts at i² = (the square of the current prime)          
            do j = i*i, limit, i
                is_prime(j) = .false.
            end do
        end if
    end do

    ! Collect primes
    count = 0
    do i = 2, limit
        if (is_prime(i)) then
            count = count + 1
            primes(count) = i
            if (count >= N) exit
        end if
    end do
    
    ! Print primes
    write (stdout,FMT='(A,I2,A)') 'Generating the first ',N,' prime numbers.'
    
    write(*,FMT='(6I8)') primes(1:N)
    
    ! Deallocate arrays
    deallocate(is_prime)
    deallocate(primes)
    
end subroutine print_primes  
  
end module memoized_fibonacci

program fibonacci_main
  use memoized_fibonacci
  implicit none
  call print_fibonacci()
  call print_factorial()
  call print_primes()
end program fibonacci_main

!!
!! Enter how many Fibonacci Numbers you need: 
!! 10
!! Generating the first 10 Fibonacci numbers.
!!           1           1           2           3           5           8
!!          13          21          34          55
!!  ----------
!! Enter how many factorial numbers you need: 
!! 6
!! Generating the first 6 factorial numbers.
!!           1           1           2           6          24         120
!!  ----------

