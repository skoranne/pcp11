module shape_functor_mod
  implicit none
  private

  ! Expose the abstract type and shape constants
  public :: ShapeFunctor
  public :: SHAPE_RECTANGLE, SHAPE_CIRCLE, SHAPE_POLYGON, SHAPE_TRAPEZOID

  ! Shape type identifiers
  integer, parameter :: SHAPE_RECTANGLE = 1
  integer, parameter :: SHAPE_CIRCLE    = 2
  integer, parameter :: SHAPE_POLYGON   = 3
  integer, parameter :: SHAPE_TRAPEZOID = 4

  ! The Abstract Functor Base Class
  type, abstract :: ShapeFunctor
   contains
     ! Deferred methods act like pure virtual functions in C++
     procedure(check_shape_if), deferred, pass(this) :: CheckShape
     procedure(calc_rect_if),   deferred, pass(this) :: CalculateAreaRectangle
     procedure(calc_circ_if),   deferred, pass(this) :: CalculateAreaCircle
  end type ShapeFunctor

  ! Abstract interfaces define the exact signatures required for the deferred methods
  abstract interface
     function check_shape_if(this, i) result(shape_type)
       import :: ShapeFunctor
       class(ShapeFunctor), intent(in) :: this
       integer, intent(in)             :: i
       integer                         :: shape_type
     end function check_shape_if

     function calc_rect_if(this, i) result(area)
       import :: ShapeFunctor
       ! Passed as intent(inout) to allow the functor to accumulate state or cache data
       class(ShapeFunctor), intent(inout) :: this
       integer, intent(in)                :: i
       real(8)                            :: area
     end function calc_rect_if

     function calc_circ_if(this, i) result(area)
       import :: ShapeFunctor
       class(ShapeFunctor), intent(inout) :: this
       integer, intent(in)                :: i
       real(8)                            :: area
     end function calc_circ_if
  end interface

end module shape_functor_mod
module shape_processor_mod
  use shape_functor_mod
  implicit none
  private
  public :: process_shape_areas

contains

  subroutine process_shape_areas(indices, functor)
    integer, intent(in)                :: indices(:)
    class(ShapeFunctor), intent(inout) :: functor

    integer :: j, idx, current_shape_type
    real(8) :: current_area, total_area

    total_area = 0.0d0

    ! Iterate through the provided shape indices
    do j = 1, size(indices)
       idx = indices(j)

       ! 1. Call the functor to identify the shape
       current_shape_type = functor%CheckShape(idx)

       ! 2. Route to the correct calculation member based on the type
       select case (current_shape_type)
       case (SHAPE_RECTANGLE)
          current_area = functor%CalculateAreaRectangle(idx)
          total_area = total_area + current_area

       case (SHAPE_CIRCLE)
          current_area = functor%CalculateAreaCircle(idx)
          total_area = total_area + current_area

       case default
          ! Shape not supported or requires a different method
          continue 
       end select

    end do

    print *, "Total computed area for subset: ", total_area

  end subroutine process_shape_areas

end module shape_processor_mod
module colored_functor_mod
  use shape_functor_mod
  implicit none
  private
  public :: ColoredGeometryFunctor

  ! Concrete type extending the abstract interface
  type, extends(ShapeFunctor) :: ColoredGeometryFunctor
     ! --- THE NEW STATE ---
     ! Scalar state for the functor (e.g., the active layer we want to process)
     integer :: color 

     ! Backing arrays for the geometry properties
     integer, allocatable :: shape_ids(:)
     integer, allocatable :: shape_colors(:) ! Color tag for each specific shape
     real(8), allocatable :: rect_widths(:), rect_heights(:)
     real(8), allocatable :: circle_radii(:)
   contains
     procedure, pass(this) :: CheckShape             => impl_check_shape
     procedure, pass(this) :: CalculateAreaRectangle => impl_calc_rect
     procedure, pass(this) :: CalculateAreaCircle    => impl_calc_circ
  end type ColoredGeometryFunctor

contains

  function impl_check_shape(this, i) result(shape_type)
    class(ColoredGeometryFunctor), intent(in) :: this
    integer, intent(in)                       :: i
    integer                                   :: shape_type
    shape_type = this%shape_ids(i)
  end function impl_check_shape

  function impl_calc_rect(this, i) result(area)
    class(ColoredGeometryFunctor), intent(inout) :: this
    integer, intent(in)                          :: i
    real(8)                                      :: area

    ! Leverage the state: Only compute if the shape's color matches the functor's target color
    if (this%shape_colors(i) == this%color) then
       area = this%rect_widths(i) * this%rect_heights(i)
       print *, "  [Rect] Index", i, "matches active color", this%color, "- Area:", area
    else
       area = 0.0d0 ! Return 0 if it doesn't match the active color state
       print *, "  [Rect] Index", i, "skipped (shape color", this%shape_colors(i), "!= active color", this%color, ")"
    end if
  end function impl_calc_rect

  function impl_calc_circ(this, i) result(area)
    class(ColoredGeometryFunctor), intent(inout) :: this
    integer, intent(in)                          :: i
    real(8)                                      :: area
    real(8), parameter                           :: PI = 3.141592653589793d0

    if (this%shape_colors(i) == this%color) then
       area = PI * (this%circle_radii(i) ** 2)
       print *, "  [Circ] Index", i, "matches active color", this%color, "- Area:", area
    else
       area = 0.0d0
       print *, "  [Circ] Index", i, "skipped (shape color", this%shape_colors(i), "!= active color", this%color, ")"
    end if
  end function impl_calc_circ

end module colored_functor_mod
