program test_utils_m
    use, intrinsic :: iso_fortran_env, only: dp => real64

    use fortuno_serial, only : execute_serial_cmd_app, is_equal, & 
       test => serial_case_item,&
       check => serial_check

    use utils_m

    implicit none

    ! Register tests by providing name and subroutine to run for each test.
    ! Note: this routine does not return but stops the program with the right exit code.
    call execute_serial_cmd_app(testitems=[&
            test("Real-space 1D grid", test_generate_real_space_grid_1D), &
            test("Real-space 2D grid", test_generate_real_space_grid_2D), &
            test("Real-space 3D grid", test_generate_real_space_grid_3D), &
            test("Linear spaced sampling", test_linspace) &
            ])
  
contains

    !> @brief Are \f$x\f$ and \f$y\f$ equal within a tolerance.
    elemental logical function close(x, y, rtol, atol)
        real(dp), intent(in) :: x, y
        real(dp), optional, intent(in) :: rtol
        real(dp), optional, intent(in) :: atol
        real(dp) :: atol_, rtol_

        if(present(rtol)) then
            rtol_ = rtol
        else
            rtol_ = 1.e-5_dp
        endif

        if(present(atol)) then
            atol_ = atol
        else
            atol_ = 1.e-8_dp
        endif

        close = abs(x - y) <= (atol_ + rtol_ * abs(y))
    end function close

    subroutine test_generate_real_space_grid_1D
        integer      :: sampling(1)
        real(dp) :: spacings(1, 1)
        real(dp), allocatable :: grid(:, :), ref_grid(:, :)
        
        sampling = 10
        spacings = 0.25_dp
        allocate(grid(1, product(sampling)), ref_grid(1, product(sampling)))

        ref_grid(1, :) = [0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50]
        call generate_real_space_grid(sampling, spacings, grid)

        call check(all(close(grid, ref_grid)), msg="1D grid disagrees with reference")

    end subroutine test_generate_real_space_grid_1D


    subroutine test_generate_real_space_grid_2D
        integer      :: sampling(2)
        real(dp) :: spacings(2, 2)
        real(dp), allocatable :: grid(:, :), ref_grid(:, :)
        integer :: n_total
        integer :: i, j, k, ir

        sampling = [3, 3]
        spacings = reshape([0.25_dp, 0.00_dp, &
                            0.00_dp, 0.25_dp], [2, 2])
        n_total = product(sampling)
        ref_grid = reshape([ &
            0.25,  0.25, &     
            0.50,  0.25, &     
            0.75,  0.25, &     
            0.25,  0.50, &     
            0.50,  0.50, &     
            0.75,  0.50, &     
            0.25,  0.75, &     
            0.50,  0.75, &     
            0.75,  0.75 ], [2, n_total])
    
        allocate(grid(2, n_total))
        call generate_real_space_grid(sampling, spacings, grid)
        call check(all(shape(grid) == [2, n_total]))
        call check(all(close(grid, ref_grid)), msg="2D grid disagrees with reference")

    end subroutine test_generate_real_space_grid_2D


    subroutine test_generate_real_space_grid_3D
        integer      :: sampling(3)
        real(dp) :: spacings(3, 3), origin(3)
        real(dp), allocatable :: grid(:, :), ref_grid(:, :)
        integer :: n_total, i_centre

        sampling = [3, 3, 3]
        spacings = reshape([1.00_dp, 0.00_dp, 0.00_dp,&
                            0.00_dp, 1.00_dp, 0.00_dp,&
                            0.00_dp, 0.00_dp, 1.00_dp], [3, 3])
        n_total = product(sampling)
        origin = [2._dp, 2._dp, 2._dp]
      
        ref_grid = reshape([ &
        -1.000000, -1.000000, -1.000000, &
         0.000000, -1.000000, -1.000000, &
         1.000000, -1.000000, -1.000000, &
        -1.000000,  0.000000, -1.000000, &
         0.000000,  0.000000, -1.000000, &
         1.000000,  0.000000, -1.000000, &
        -1.000000,  1.000000, -1.000000, &
         0.000000,  1.000000, -1.000000, &
         1.000000,  1.000000, -1.000000, &
        -1.000000, -1.000000,  0.000000, &
         0.000000, -1.000000,  0.000000, &
         1.000000, -1.000000,  0.000000, &
        -1.000000,  0.000000,  0.000000, &
         0.000000,  0.000000,  0.000000, &
         1.000000,  0.000000,  0.000000, &
        -1.000000,  1.000000,  0.000000, &
         0.000000,  1.000000,  0.000000, &
         1.000000,  1.000000,  0.000000, &
        -1.000000, -1.000000,  1.000000, &
         0.000000, -1.000000,  1.000000, &
         1.000000, -1.000000,  1.000000, &
        -1.000000,  0.000000,  1.000000, &
         0.000000,  0.000000,  1.000000, &
         1.000000,  0.000000,  1.000000, &
        -1.000000,  1.000000,  1.000000, &
         0.000000,  1.000000,  1.000000, &
         1.000000,  1.000000,  1.000000  &
         ], [3, n_total])
    
        allocate(grid(3, n_total))
        call generate_real_space_grid(sampling, spacings, grid, origin=origin)
        
        call check(all(shape(grid) == [3, n_total]))

        call check(all(close(grid, ref_grid)), msg="3D grid disagrees with reference")

        i_centre = 14
        call check(all(close(grid(:, i_centre), [0._dp, 0._dp, 0._dp])), msg="Grid centred on the origin")
   
    end subroutine test_generate_real_space_grid_3D


    subroutine test_linspace
        real(dp), allocatable :: x(:), y(:)

        allocate(x(10))
        call linspace(1._dp, 10.0_dp, size(x), x, end_point=.true.)
        call check(all(close(x, [1._dp, 2._dp, 3._dp, 4._dp, 5._dp, 6._dp, 7._dp, 8._dp, 9._dp, 10._dp])))

        allocate(y(9))
        call linspace(1._dp, 10.0_dp, size(y), y, end_point=.false.)
        call check(all(close(x, [1._dp, 2._dp, 3._dp, 4._dp, 5._dp, 6._dp, 7._dp, 8._dp, 9._dp])))


    end subroutine test_linspace



end program test_utils_m
