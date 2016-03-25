program test

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT, &
        compiler_version, &
        compiler_options

    use lagrange_interpolator_library, only: &
        LagrangeInterpolator

    ! Explicit typing only
    implicit none

    !----------------------------------------------------------------------
    ! Dictionary
    !----------------------------------------------------------------------
    integer (ip)               :: i, j !! Counters
    integer (ip), parameter    :: NX=42
    integer (ip), parameter    :: NY=NX
    real (wp),    parameter    :: TOLERANCE=1.0e-16_wp
    real (wp)                  :: x(NX), y(NY)
    real (wp)                  :: func_1d(NX)
    real (wp)                  :: func_2d(NX,NY)
    real (wp)                  :: approx_value(2)
    real (wp)                  :: max_error(2)
    type(LagrangeInterpolator) :: interp
    !----------------------------------------------------------------------

    ! Calculate x-grid
    x = [ (real(i - 1, kind=wp)/(NX-1), i=1,NX) ]

    ! Calculate y-grid
    y = [(real(j - 1, kind=wp)/(NY-1), j=1,NY)]

    ! Calculate exact values
    do i=1,NX
        func_1d(i) = f1(x(i))
        do j=1,NY
            func_2d(i,j) = f2(x(i),y(j))
        end do
    end do


    ! compute max error at interpolation points
    max_error = 0.0_wp
    do i=1,NX
        call interp%perform_1d_interpolation( 0, x, func_1d, x(i), approx_value(1) )
        associate( exact_value => f1(x(i)) )
            associate( approximation_error => abs(exact_value-approx_value(1)) )
                max_error(1) = max(approximation_error,max_error(1))
            end associate
        end associate
        do j=1,NY
            associate( desired_point => [x(i), y(j)] )
                call interp%perform_2d_interpolation( 0, x, y, func_2d, desired_point, approx_value(2) )
            end associate
            associate( exact_value => f2(x(i),y(j)) )
                associate( approximation_error => abs(exact_value-approx_value(2)) )
                    max_error(2) = max(approximation_error,max_error(2))
                end associate
            end associate
        end do
    end do

    ! Check discretization error against tolerance
    do i = 1, size(max_error)
        write( stdout, * ) i,'D interpolation: max error =', max_error(i)
        if ( max_error(i) >= TOLERANCE ) then
            write( stdout, '(A)' )  ' ** Test failed ** '
        else
            write( stdout, '(A)' )  ' ** Test passed ** '
        end if
        write( stdout, '(A)' ) ''
    end do

    ! Print compiler info
    write( stdout, '(A)' ) ' '
    write( stdout, '(4A)' ) 'This result was compiled by ', &
        compiler_version(), ' using the options ', &
        compiler_options()
    write( stdout, '(A)' ) ' '


contains


    pure function f1(x) result (return_value) !! 1d test function
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real (wp), intent (in) :: x
        real (wp)              :: return_value
        !----------------------------------------------------------------------

        return_value = 0.5_wp * (x*exp(-x) + x*sin(x) )

    end function f1


    pure function f2(x,y) result (return_value) !! 2d test function
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real (wp), intent (in) :: x
        real (wp), intent (in) :: y
        real (wp)              :: return_value
        !----------------------------------------------------------------------

        associate( half_pi => acos(-1.0_wp)/2 )
            return_value = 0.5_wp * (y*exp(-x) + sin(half_pi * y) ) * x
        end associate

    end function f2

end program test
