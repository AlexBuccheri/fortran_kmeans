program test_grids_m
    use, intrinsic :: iso_fortran_env, only: dp => real64

    ! Get via wrapper
    ! use fortuno_serial, only : execute_serial_cmd_app, is_equal, & 
    !    test => serial_case_item,&
    !    check => serial_check

    use fortuno_interface_m, only: execute_cmd_app, test, check, is_equal
    use maths_m, only: all_close
    use grids_m, only: generate_real_space_grid, linspace, linspace_to_grid, generate_gaussian

    implicit none

    ! Register tests by providing name and subroutine to run for each test.
    ! Note: this routine does not return but stops the program with the right exit code.
    call execute_cmd_app(testitems=[&
            test("Real-space 1D grid", test_generate_real_space_grid_1D), &
            test("Real-space 2D grid", test_generate_real_space_grid_2D), &
            test("Real-space 3D grid", test_generate_real_space_grid_3D), &
            test("Linear spaced sampling", test_linspace), &
            test("2D linear grid", test_linspace_to_grid2d), & 
            test("2D Gaussian", test_generate_gaussian_2d) &
            ])
  
contains


    subroutine test_generate_real_space_grid_1D
        integer      :: sampling(1)
        real(dp) :: spacings(1, 1)
        real(dp), allocatable :: grid(:, :), ref_grid(:, :)
        
        sampling = 10
        spacings = 0.25_dp
        allocate(grid(1, product(sampling)), ref_grid(1, product(sampling)))

        ref_grid(1, :) = [0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50]
        call generate_real_space_grid(sampling, spacings, grid)

        call check(all_close(grid, ref_grid), msg="1D grid disagrees with reference")

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
        call check(all_close(grid, ref_grid), msg="2D grid disagrees with reference")

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

        call check(all_close(grid, ref_grid), msg="3D grid disagrees with reference")

        i_centre = 14
        call check(all_close(grid(:, i_centre), [0._dp, 0._dp, 0._dp]), msg="Grid centred on the origin")
   
    end subroutine test_generate_real_space_grid_3D


    subroutine test_linspace
        real(dp), allocatable :: x(:), y(:)

        allocate(x(10))
        call linspace(1._dp, 10.0_dp, size(x), x, end_point=.true.)
        call check(all_close(x, [1._dp, 2._dp, 3._dp, 4._dp, 5._dp, 6._dp, 7._dp, 8._dp, 9._dp, 10._dp]))

        allocate(y(9))
        call linspace(1._dp, 10.0_dp, size(y), y, end_point=.false.)
        call check(all_close(y, [1._dp, 2._dp, 3._dp, 4._dp, 5._dp, 6._dp, 7._dp, 8._dp, 9._dp]))

    end subroutine test_linspace


    subroutine test_linspace_to_grid2d
        real(dp) :: x(5), y(5)
        real(dp), allocatable :: grid(:, :), ref_grid(:, :)
        integer :: i

        allocate(grid(2, 25), ref_grid(2, 25))
        ref_grid = reshape([&
          1.0_dp,  1.00_dp, &    
          2.0_dp,  1.00_dp, &    
          3.0_dp,  1.00_dp, &    
          4.0_dp,  1.00_dp, &    
          5.0_dp,  1.00_dp, &    
          1.0_dp,  3.25_dp, &    
          2.0_dp,  3.25_dp, &    
          3.0_dp,  3.25_dp, &    
          4.0_dp,  3.25_dp, &    
          5.0_dp,  3.25_dp, &    
          1.0_dp,  5.50_dp, &    
          2.0_dp,  5.50_dp, &    
          3.0_dp,  5.50_dp, &    
          4.0_dp,  5.50_dp, &    
          5.0_dp,  5.50_dp, &    
          1.0_dp,  7.75_dp, &    
          2.0_dp,  7.75_dp, &    
          3.0_dp,  7.75_dp, &    
          4.0_dp,  7.75_dp, &    
          5.0_dp,  7.75_dp, &    
          1.0_dp,  10.0_dp, &    
          2.0_dp,  10.0_dp, &    
          3.0_dp,  10.0_dp, &    
          4.0_dp,  10.0_dp, &    
          5.0_dp,  10.0_dp], [2, 25])

        call linspace(1._dp,  5.0_dp, size(x), x)
        call linspace(1._dp, 10.0_dp, size(y), y)
        call linspace_to_grid(x, y, grid)

        call check(all_close(grid, ref_grid))

    end subroutine test_linspace_to_grid2d


    ! TODO(Alex) Note, this is a demonstration that can be plot, but it does not assert anything
    subroutine test_generate_gaussian_2d
        real(dp), allocatable :: x(:), y(:), grid(:, :)
        real(dp), allocatable :: mean(:, :), sigma(:)
        real(dp), allocatable :: gaussian(:), total_gaussian(:)
        integer :: nx, ny, ir, ix, iy, i

        nx = 50
        ny = 50
        allocate(x(nx), y(ny), grid(2, nx * ny))
        call linspace(1._dp, 5.0_dp, nx, x)
        call linspace(1._dp, 5.0_dp, ny, y)
        call linspace_to_grid(x, y, grid)

        ! Create a gaussian in each quadrant of the 2D grid
        allocate(gaussian(nx * ny), total_gaussian(nx * ny))
        sigma = [0.4_dp, 0.4_dp]
        mean = reshape([2._dp, 2._dp, &
                        2._dp, 4._dp, &
                        4._dp, 2._dp, &
                        4._dp, 4._dp], [2, 4])

        total_gaussian = 0._dp
        do i = 1, 4
            call generate_gaussian(mean(:, i), sigma, grid, gaussian)
            total_gaussian = total_gaussian + gaussian
        enddo

        ! Output for gnuplot
        ! ir = 0
        ! do iy = 1 , ny
        !     do ix = 1, nx
        !         ir = ir + 1
        !         write(100, *) grid(:, ir), total_gaussian(ir)
        !     enddo
        !     write(100, *)
        ! enddo

    end subroutine test_generate_gaussian_2d


    ! TODO(Alex) test_linspace_to_grid3d


end program test_grids_m
