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
end module memoized_fibonacci

program fibonacci_main
  use memoized_fibonacci
  implicit none
  call print_fibonacci()
  call print_factorial()
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

