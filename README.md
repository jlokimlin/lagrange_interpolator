# **lagrange\_interpolator**

This Fortran project implements 1 and 2-dimensional Lagrange interpolation.


-----------------------------------------------------------------------------

## Requirements
* The GNU Make tool https://www.gnu.org/software/make/
* The GNU gfortran compiler https://gcc.gnu.org/wiki/GFortran

-----------------------------------------------------------------------------


## To build the project

Type the following command line arguments
```
git clone https://github.com/jlokimlin/lagrange_interpolator.git

cd lagrange_interpolator; make all
```

-----------------------------------------------------------------------------

## Usage

```fortran

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use lagrange_interpolator_library, only: &
        LagrangeInterpolator

    ! Explicit typing only
    implicit none
    
    type (LagrangeInterpolator) :: interp
    real (wp)				    :: estimate
    real (wp)					:: desired_point, desired_point_2d(2)
    real (wp), allocatable      :: x(:), y(:), fx(:), fxy(:,:)
    integer (ip), parameter	    :: TYPE = 0 ! 0 = None, 1 = Constant, 2 = Linear
    
    !... generate some data
    
    call interp%perform_1d_interpolation( TYPE, x, fx, desired_point, estimate )
    
    call interp%perform_2d_interpolation( TYPE, x, y, fxy, desired_point_2d, estimate )

```