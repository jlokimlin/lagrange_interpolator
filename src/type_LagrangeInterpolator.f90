module type_LagrangeInterpolator

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stderr => ERROR_UNIT

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: LagrangeInterpolator

    !--------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !--------------------------------------------------------------------------
    integer (ip), parameter :: EXTRAPOLATION_NONE = 0
    integer (ip), parameter :: EXTRAPOLATION_CONSTANT = 1
    integer (ip), parameter :: EXTRAPOLATION_LINEAR = 2
    !--------------------------------------------------------------------------

    ! Declare derived data type
    type, public :: LagrangeInterpolator
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        logical,                private  :: initialized = .false.
        logical,                public   :: success = .false.
        integer (ip),           private  :: EXTRAPOLATION_TYPE = 0
        real (wp), allocatable, private  :: input_data(:)
        real (wp), allocatable, private  :: range_data(:)
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure,         public  :: perform_1d_interpolation
        procedure,         public  :: perform_2d_interpolation
        procedure, nopass, public  :: perform_naive_interpolation
        procedure,         private :: create => create_lagrange_interpolator
        procedure,         private :: destroy => destroy_lagrange_interpolator
        final                      :: finalize_lagrange_interpolator
        !----------------------------------------------------------------------
    end type LagrangeInterpolator


contains


    subroutine create_lagrange_interpolator( this, x, y, extrapolation )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (LagrangeInterpolator), intent (in out) :: this
        real (wp),                    intent (in)     :: x(:)
        real (wp),                    intent (in)     :: y(:)
        integer (ip),                 intent (in)     :: extrapolation
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: i !! counter
        !----------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Set the extrapolation method
        select case (extrapolation)
            case(EXTRAPOLATION_CONSTANT)
                this%EXTRAPOLATION_TYPE = EXTRAPOLATION_CONSTANT
            case(EXTRAPOLATION_LINEAR)
                this%EXTRAPOLATION_TYPE = EXTRAPOLATION_LINEAR
            case default
                this%EXTRAPOLATION_TYPE = EXTRAPOLATION_NONE
        end select


        ! Check if sufficient data was supplied
        if ( size(x) < 2 .or. size(y) < size(x) ) then
            error stop  'TYPE (LagrangeInterpolator): '&
                //'Insufficient data in CREATE_LAGRANGE_INTERPOLATOR'&
                //'Either size(x) < 2 or size(y) < size(x)  '
        end if

        ! Sort data
        this%success = .true.
        associate( n => size(x) )
            do i = 2, n
                if ( x(i) < x(i-1) ) then
                    this%success = .false.
                    exit
                end if
            end do
            ! Return if data isn't sorted
            if ( .not. this%success ) then
                return
            end if

            ! Allocate memory
            allocate( this%input_data(n) )
            allocate( this%range_data(n) )

            ! Remark: we allow range data y to be larger than input x
            this%input_data(1:n) = x(1:n)
            this%range_data(1:n) = y(1:n)
        end associate

        ! Set initialization flag
        this%initialized = .true.

    end subroutine create_lagrange_interpolator


    subroutine destroy_lagrange_interpolator( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (LagrangeInterpolator), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check flag
        if ( this%initialized  .eqv. .false. ) return

        ! Release memory
        if (allocated(this%input_data)) deallocate( this%input_data )
        if (allocated(this%range_data)) deallocate( this%range_data )

        ! Reset constants
        this%EXTRAPOLATION_TYPE = 0

        ! Reset booleans
        this%initialized = .false.
        this%success = .false.

    end subroutine destroy_lagrange_interpolator


    subroutine perform_1d_interpolation( this, &
        extrapolation_type, x_data, fx_data, desired_point, estimate, success )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (LagrangeInterpolator), intent (in out) :: this
        integer (ip),                 intent (in)     :: extrapolation_type
        real (wp),                    intent (in)     :: x_data(:)
        real (wp),                    intent (in)     :: fx_data(:)
        real (wp),                    intent (in)     :: desired_point
        real (wp),                    intent (out)    :: estimate
        logical, optional,            intent (out)    :: success
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: i, idx
        !----------------------------------------------------------------------

        ! initialize estimate and flag
        estimate = 0.0_wp
        this%success = .false.

        ! initialize interpolation object
        call this%create( x_data, fx_data, extrapolation_type )

        ! Perform interpolation
        associate( &
            n => size(this%input_data), &
            x => this%input_data, &
            y => this%range_data, &
            xp => desired_point &
            )

            ! Check extrapolation
            if ( extrapolation_type == EXTRAPOLATION_NONE ) then
                if ( xp < x(1)  ) return
                if ( xp > x(n) ) return
            end if

            if ( extrapolation_type == EXTRAPOLATION_CONSTANT ) then
                if ( xp < x(1) ) then
                    estimate = x(1)
                    return
                end if
                if ( xp > x(n) ) then
                    estimate = x(n)
                    return
                end if
            end if

            ! Search for the interval that contains xp
            ! Linear extrapolation is taken care of automatically
            idx = n-1
            do i = 2, n-1
                if ( xp < x(i) ) then
                    idx = i-1
                    exit
                end if
            end do

            associate( dx => x(idx+1) - x(idx) )
                if ( dx /= 0.0_wp ) then
                    estimate = y(idx) + (xp - x(idx))  * (y(idx+1) - y(idx)) / dx
                else
                    estimate =  0.5_wp * (y(idx+1) + y(idx))
                end if
            end associate
        end associate

        ! Set interpolation status
        this%success = .true.

        ! Address optional argument
        if (present(success)) success = this%success

        ! Release memory
        call this%destroy()

    end subroutine perform_1d_interpolation



    subroutine perform_2d_interpolation( this, &
        extrapolation_type, x_data, y_data, fxy_data, desired_point, estimate )
        !
        ! Purpose:
        !
        !    Given arrays X_DATA(1:n) and Y_DATA(1:m) of independent variables, and
        !    an m-by-n array of function values FXY_DATA(1:n, 1:m), tabulated at the
        !    grid points defined by X_DATA and Y_DATA; and given values
        !    DESIRED_POINT = (x1, x2) of the independent variables;
        !    this routine returns an interpolated function value ESTIMATE
        !
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (LagrangeInterpolator), intent (in out) :: this
        integer (ip),                 intent (in)     :: extrapolation_type
        real (wp),                    intent (in)     :: x_data(:)
        real (wp),                    intent (in)     :: y_data(:)
        real (wp),                    intent (in)     :: fxy_data(:,:)
        real (wp),                    intent (in)     :: desired_point(:)
        real (wp),                    intent (out)    :: estimate
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)              :: i !! Counter
        real (wp), allocatable    :: x_estimate(:)
        real (wp), allocatable    :: y_estimate(:)
        !----------------------------------------------------------------------

        ! Associate dimension
        associate( n => size(x_data) )

            ! Allocate memory
            allocate( x_estimate(n) )
            allocate( y_estimate(n) )

            ! Perform interpolation
            do i = 1, n ! loop over rows
                ! copy the rows into temporary storage
                y_estimate = fxy_data(i,:)
                ! interpolate answer into temporary storage
                call this%perform_1d_interpolation( &
                    extrapolation_type, y_data, y_estimate, desired_point(2), x_estimate(i) )
                ! do the final interpolation
                call this%perform_1d_interpolation( &
                    extrapolation_type, x_data, x_estimate, desired_point(1), estimate )
            end do
        end associate

        ! Release memory
        deallocate( x_estimate )
        deallocate( y_estimate )

    end subroutine perform_2d_interpolation



    function perform_naive_interpolation( x, y, xp ) result ( return_value )
        !
        ! Purpose:
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real (wp), intent (in)  :: x(:)
        real (wp), intent (in)  :: y(:)
        real (wp), intent (in)  :: xp
        real (wp)               :: return_value
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)  :: i, idx
        !----------------------------------------------------------------------

        ! Search for the interval that contains xp
        idx = size(x) - 1

        do i = 2, size(x) - 1
            if ( xp < x(i) ) then
                idx = i - 1
                exit
            end if
        end do

        return_value = &
            y(idx) + (xp - x(idx)) * (y(idx+1) - y(idx)) / &
            (x(idx+1) -  x(idx))

    end function perform_naive_interpolation


    subroutine finalize_lagrange_interpolator( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (LagrangeInterpolator), intent (in out) :: this
        !----------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_lagrange_interpolator


end module type_LagrangeInterpolator
