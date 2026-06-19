program main
    use shape_functor_mod
    use shape_processor_mod
    use colored_functor_mod
    implicit none

    type(ColoredGeometryFunctor) :: my_functor
    integer, allocatable         :: process_subset(:)

    ! 1. Allocate backing store 
    allocate(my_functor%shape_ids(3))
    allocate(my_functor%shape_colors(3))
    allocate(my_functor%rect_widths(3), my_functor%rect_heights(3))
    allocate(my_functor%circle_radii(3))

    ! Shape 1: Rectangle (color 1 - e.g., Metal1)
    my_functor%shape_ids(1)    = SHAPE_RECTANGLE
    my_functor%shape_colors(1) = 1
    my_functor%rect_widths(1)  = 2.0d0
    my_functor%rect_heights(1) = 3.0d0

    ! Shape 2: Circle (color 2 - e.g., Via)
    my_functor%shape_ids(2)    = SHAPE_CIRCLE
    my_functor%shape_colors(2) = 2
    my_functor%circle_radii(2) = 5.0d0

    ! Shape 3: Rectangle (color 1 - e.g., Metal1)
    my_functor%shape_ids(3)    = SHAPE_RECTANGLE
    my_functor%shape_colors(3) = 1
    my_functor%rect_widths(3)  = 4.0d0
    my_functor%rect_heights(3) = 4.0d0

    ! Set up the indices to process
    allocate(process_subset(3))
    process_subset = [1, 2, 3]

    ! --- Pass 1: Set Functor State to Color 1 ---
    print *, "--- Processing subset with Functor State Color = 1 ---"
    my_functor%color = 1 
    call process_shape_areas(process_subset, my_functor)

    print *, ""

    ! --- Pass 2: Set Functor State to Color 2 ---
    print *, "--- Processing subset with Functor State Color = 2 ---"
    my_functor%color = 2 
    call process_shape_areas(process_subset, my_functor)

end program main
