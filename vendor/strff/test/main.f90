! Generated by cart. DO NOT EDIT
program main

#ifdef USE_CAFFEINE
   use caffeine_m, only : stop => caf_stop
#endif

    implicit none

    if (.not.run()) stop 1
contains
    function run() result(passed)
        use add_hanging_indentation_test, only: &
                add_hanging_indentation_add_hanging_indentation => &
                    test_add_hanging_indentation
        use format_hanging_indented_test, only: &
                format_hanging_indented_format_hanging_indented => &
                    test_format_hanging_indented
        use indent_test, only: &
                indent_indent => &
                    test_indent
        use join_test, only: &
                join_join => &
                    test_join
        use read_file_lines_test, only: &
                read_file_lines_read_file_lines => &
                    test_read_file_lines
        use read_file_test, only: &
                read_file_read_file => &
                    test_read_file
        use split_at_test, only: &
                split_at_split_at => &
                    test_split_at
        use starts_with_test, only: &
                starts_with_starts_with => &
                    test_starts_with
        use to_string_test, only: &
                to_string_to_string_for_doubles => &
                    test_to_string_for_doubles, &
                to_string_to_string_for_integers => &
                    test_to_string_for_integers
        use veggies, only: test_item_t, test_that, run_tests



        logical :: passed

        type(test_item_t) :: tests
        type(test_item_t) :: individual_tests(10)

        individual_tests(1) = add_hanging_indentation_add_hanging_indentation()
        individual_tests(2) = format_hanging_indented_format_hanging_indented()
        individual_tests(3) = indent_indent()
        individual_tests(4) = join_join()
        individual_tests(5) = read_file_lines_read_file_lines()
        individual_tests(6) = read_file_read_file()
        individual_tests(7) = split_at_split_at()
        individual_tests(8) = starts_with_starts_with()
        individual_tests(9) = to_string_to_string_for_doubles()
        individual_tests(10) = to_string_to_string_for_integers()
        tests = test_that(individual_tests)


        passed = run_tests(tests)

    end function
end program
