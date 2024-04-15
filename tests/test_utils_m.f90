program test_utils_m
    use, intrinsic :: iso_fortran_env

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
            test("Real-space 3D grid", test_generate_real_space_grid_3D) &
            ])
  
contains

    !> @brief Are \f$x\f$ and \f$y\f$ equal within a tolerance.
    elemental logical function close(x, y, rtol, atol)
        real(real64), intent(in) :: x, y
        real(real64), optional, intent(in) :: rtol
        real(real64), optional, intent(in) :: atol
        real(real64) :: atol_, rtol_

        if(present(rtol)) then
            rtol_ = rtol
        else
            rtol_ = 1.e-5_real64
        endif

        if(present(atol)) then
            atol_ = atol
        else
            atol_ = 1.e-8_real64
        endif

        close = abs(x - y) <= (atol_ + rtol_ * abs(y))
    end function close


    subroutine test_generate_real_space_grid_1D
        integer      :: sampling(1)
        real(real64) :: spacings(1, 1)
        real(real64), allocatable :: grid(:, :), ref_grid(:, :)
        
        sampling = 10
        spacings = 0.25_real64
        allocate(grid(1, product(sampling)), ref_grid(1, product(sampling)))

        ref_grid(1, :) = [0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50]
        call generate_real_space_grid(sampling, spacings, grid)

        call check(all(close(grid, ref_grid)), msg="1D grid disagrees with reference")

    end subroutine test_generate_real_space_grid_1D


    subroutine test_generate_real_space_grid_2D
        integer      :: sampling(2)
        real(real64) :: spacings(2, 2)
        real(real64), allocatable :: grid(:, :), ref_grid(:, :)
        integer :: n_total
        integer :: i, j, k, ir

        sampling = [3, 3]
        spacings = reshape([0.25_real64, 0.00_real64, &
                            0.00_real64, 0.25_real64], [2, 2])
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
        real(real64) :: spacings(3, 3), origin(3)
        real(real64), allocatable :: grid(:, :), ref_grid(:, :)
        integer :: n_total, i_centre

        sampling = [3, 3, 3]
        spacings = reshape([1.00_real64, 0.00_real64, 0.00_real64,&
                            0.00_real64, 1.00_real64, 0.00_real64,&
                            0.00_real64, 0.00_real64, 1.00_real64], [3, 3])
        n_total = product(sampling)
        origin = [2._real64, 2._real64, 2._real64]
      
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
        call check(all(close(grid(:, i_centre), [0._real64, 0._real64, 0._real64])), msg="Grid centred on the origin")
   
    end subroutine test_generate_real_space_grid_3D

end program test_utils_m
